#include <array>
#include <cmath>

#include <Eigen/Core>

#include "common/Constants.hpp"

#if T2DS_DEBUG
#include "Utils/Logger.hxx"
#include "Utils/Utilities.hxx"
#endif
#include "KalmanFitter/KalmanFitterAlgebra.hxx"
#include "KalmanFitter/KalmanFitterTransport.hxx"
#include "Seeder/SeederHelixVertex.hxx"
#include "Seeder/SeederLineVertex.hxx"
#include "Seeder/SeederTypes.hxx"

#include "KalmanFitter/KalmanFitterProdVertexConstraint.hxx"

// == Production Vertex Constraint == //

namespace T2DS::KF::Internal {

namespace {

// Seed a particle onto an arbitrary point, picking the line or the helix route.
[[nodiscard]] Seeder::Result SeedToPoint(const Particle& part, const std::array<double, 3>& xyz, double bz) {

    const Seeder::State s0{.x = part.X(), .y = part.Y(), .z = part.Z(), .px = part.Px(), .py = part.Py(), .pz = part.Pz(), .charge = part.Charge()};

    if (s0.charge == 0 || std::abs(bz) < Common::AbsAlmostZero) return Seeder::LineVertex::FullPCA(s0, xyz);
    return Seeder::HelixVertex::FullPCA(s0, xyz, bz);
}

}  // namespace

// Add `vtx` as a point-like measurement to `part`.
// The returned copy has its parameters *at* the production vertex, and the eighth component of its state vector filled
// with the decay length to momentum ratio (s = l/p) back to `part`'s decay vertex, covariances included.
// The constraint's own contribution is `out.Chi2() - part.Chi2()`, with 2 d.o.f.
Particle AtProductionVertex(const Particle& part, const Vertex& vtx, double bz) {

    // remember where the particle decayed //

    const Eigen::Vector<double, 3> decay_xyz = part.fP.head<3>();
    const Eigen::Matrix<double, 3, 3> decay_cov = part.fC.topLeftCorner<3, 3>();

    // transport the particle to the vertex //

    const Seeder::Result s = SeedToPoint(part, vtx.GetXYZ(), bz);
    const TransportOp op = TransportOp::From(part, s, bz);

    Particle kf = part;
    SetStateAtPCA(kf, s.seed.pca);

    // -- cov(transported state, decay point) = F * C0[0:6, 0:3]. Kept for the var(S) cross term below; it has to
    //    be taken before `fC` is overwritten, and the vertex contributes nothing here (it is assumed independent
    //    of the candidate -- see the note on `fC(7,7)`).
    const Eigen::Matrix<double, 6, 3> mA = op.jacob * part.fC.block<6, 3>(0, 0);

    op.PropagateCov(kf.fC);

    // -- correlation between the transported particle and the vertex: D[i][j] = sum_k V[j,k] * F1[i,k]
    // -- fold in the vertex's own error, F1 * V * F1^T. `corr`'s columns 3-5 are zero (`ds_dr1` only has position
    //    components), so `V` only ever enters through its 3x3 block; and `corr` is rank-1, so both products
    //    collapse onto the same outer product `df3 * df3'`.
    //    NOTE: the cross term of `TransportToPCA` is deliberately not reused here: it adds the whole 6x6, whereas
    //          the original writes back `fC[0..5]` only, i.e. the position block. For a neutral mother the two are
    //          identical anyway -- `bq = 0` makes rows 3-5 of `corr` vanish as well.
    const Eigen::Vector<double, 3> df3 = op.df_ds.head<3>();
    const Eigen::Vector<double, 3> ds_dvtx = op.ds_dr1.head<3>();
    const Eigen::Vector<double, 3> vw = vtx.cov * ds_dvtx;  // `vtx.cov` is full symmetric, so this is also ds_dvtx' * V

    const Eigen::Matrix<double, 3, 3> mD = df3 * vw.transpose();
    kf.fC.topLeftCorner<3, 3>() += vw.dot(ds_dvtx) * df3 * df3.transpose();

    // update with the vertex as a position measurement //

    const Eigen::Matrix<double, 3, 3> mS = kf.fC.topLeftCorner<3, 3>() + vtx.cov;
    const Eigen::Matrix<double, 3, 3> mS_inv = Algebra::InvertSym3(mS);

    // -- residual = measured - estimated
    const Eigen::Vector<double, 3> res = vtx.xyz - kf.fP.head<3>();

    const Eigen::Matrix<double, 7, 3> mCHt = kf.fC.block<7, 3>(0, 0);

    // -- see the note in `VertexUpdate()`: S is symmetric, so both transposes drop out
    const Eigen::Matrix<double, 7, 3> mK = mCHt * mS_inv;

    // -- both `res` and `mS` are the pre-update ones, as in the original
    const double dchi2 = res.dot(mS_inv * res);

    // -- cov(updated state, decay point) = (1 - K * H) * mA, with H picking the position block
    const Eigen::Matrix<double, 6, 3> mG = mA - mK.topRows<6>() * mA.topRows<3>();

    kf.fP.head<7>().noalias() += mK * res;
    kf.fC.topLeftCorner<7, 7>().noalias() -= mK * mCHt.transpose();

    // correlation correction //

    // -- NOTE: `mD` goes in transposed w.r.t. the vertex updates' `cross`. Both mirror their own original, and the
    //          two `D`s carry opposite index orders: `AddDaughterWithEnergyFit` has `D[i][k] * K2[k][j]`,
    //          `SetProductionVertex` has `D[k][i] * K2[k][j]`.
    Algebra::ApplyCrossCorrection(kf.fC.topLeftCorner<3, 3>(), mD.transpose(), mK.topRows<3>());

    // decay length //

    // -- now that the particle sits at the production vertex, measure back to where it decayed
    const Seeder::Result back = SeedToPoint(kf, {decay_xyz(0), decay_xyz(1), decay_xyz(2)}, bz);
    const Eigen::Map<const Eigen::Vector<double, 6>> ds_dr(back.deriv.ds_dr.data());

    kf.fP(7) = back.seed.ds;

    // -- S depends on the state *and*, explicitly, on the decay point: d(S)/d(decay point) = -`ds_dr`'s position
    //    part, by translation invariance. So, with j = `ds_dr` and j3 = its position part,
    //        cov(S, state) = C * j - G * j3
    //        var(S)        = j' * C * j - 2 * j' * G * j3 + j3' * V_decay * j3
    //    NOTE: the original drops the middle term, treating the decay point as independent of the PV-constrained
    //          state. That is *exact before* the update -- the transport lands on the PCA, so `ds_dr`'s position
    //          part is parallel to p and the state's dependence on the decay point projects out -- but the update
    //          mixes in a K-weighted piece of the position measurement and breaks it. The dropped term vanishes
    //          again as K -> 0 and as K -> 1, so it is a correction, not a factor of two; it is restored here
    //          because `DecayLengthErr()` feeds a significance that gets cut on.
    //    IMPORTANT: this still assumes the vertex is independent of the candidate. If the granddaughter tracks
    //               entered the vertex fit upstream, that correlation is of the same order and nothing in here
    //               can recover it.
    const Eigen::Vector<double, 6> cds = kf.fC.topLeftCorner<6, 6>() * ds_dr - mG * ds_dr.head<3>();

    kf.fC.block<1, 6>(7, 0) = cds.transpose();
    kf.fC.block<6, 1>(0, 7) = cds;                                     // NOTE: `fC` is full symmetric, so both halves are written
    kf.fC(7, 7) = ds_dr.dot(cds)                                       //  j' * C * j - j' * G * j3
                  - ds_dr.head<3>().dot(mG.transpose() * ds_dr)        // -j3' * G' * j
                  + ds_dr.head<3>().dot(decay_cov * ds_dr.head<3>());  // +j3' * V_decay * j3

    // -- NOTE: `fC(7,6)`, the S-E covariance, is left alone. The original doesn't fill it either (its loop stops at
    //          `fC[33]`), so it stays at the 0 the transport produced.

    kf.fChi2 += dchi2;
    kf.fNDF += 2;

#if T2DS_DEBUG
    Logger::Debug(__FUNCTION__, "mS     = {}", mS);
    Logger::Debug(__FUNCTION__, "res    = {}", res);
    Logger::Debug(__FUNCTION__, "mK     = {}", mK);
    Logger::Debug(__FUNCTION__, "mD     = {}", mD);
    Logger::Debug(__FUNCTION__, "mG     = {}", mG);
    Logger::Debug(__FUNCTION__, "dchi2  = {:13.6e}", dchi2);
    Logger::Debug(__FUNCTION__, "out.fP = {}", kf.fP);
    Logger::Debug(__FUNCTION__, "out.fC = {}", kf.fC);
#endif

    return kf;
}

}  // namespace T2DS::KF::Internal
