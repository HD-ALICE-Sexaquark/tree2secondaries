#include <cmath>
#include <optional>

#include <Eigen/Eigen>

#include "common/Constants.hpp"

#include "Seeder/BaseSeeder.hxx"
#include "Seeder/SeederHelixVertex.hxx"
#include "Seeder/SeederLineVertex.hxx"
#if T2DS_DEBUG
#include "App/Logger.hxx"
#include "App/Utilities.hxx"
#endif

#include "KalmanFitter/BaseKalmanFitter.hxx"

namespace T2DS::KF {

// == Daughter Struct == //

// Calculate the transport Jacobian = d(fP new)/d(fP old) and the correlation matrix = d(fP new)/d(r1),
// given
// - this parametrization:
//     x' = x0 + sB * px0 + cB * py0
//     y' = y0 - cB * px0 + sB * py0
//     z' = z0 + dS * pz0
//     px' = c * px0 + s * py0
//     py' = -s * px0 + c * py0
//     pz' = pz0,
// - these cached `Seed` variables:
//     sB = sin(bq * ds)/bq
//     cB = (1 - cos(bq * ds))/bq
//     cos = cos(bq * ds)
//     sin = sin(bq * ds),
// - and the d(ds)/d(r0) derivatives from `Deriv` for (x0,y0,z0,px0,py0,z0), respectively.
void Daughter::PrepareJacobAndCorr(const Seeder::Result& s, double bz) {

    // get initial values //

    double bq = bz * static_cast<double>(Charge()) * Common::Kappa;
    double px0 = Px();
    double py0 = Py();
    double pz0 = Pz();

    // prepare Jacobian Matrix //

    jacob = Eigen::Matrix<double, 8, 8>::Identity();

    // for example, when (i,j)=(0,3)=(x',px0):
    // d(x')/d(px0) = d(sB * px0)/d(px0) + d(cB * py0)/d(px0)
    //              = sB + px0 * d(sB)/d(px0) + py0 * d(cB)/d(px0)
    //              = sB + px0 * d(sB)/d(ds) * d(ds)/d(px0) + py0 * d(cB)/d(ds) * d(ds)/d(px0)
    //              = sB + [ px0 * d(sB)/d(ds) + py0 * d(cB)/d(ds) ] * d(ds)/d(px0)
    //              = sB + d(x')/d(ds) * d(ds)/d(px0)
    //              = FC + DF_DS * DS_DR

    // -- first components (non-zero and non-unity)
    jacob(0, 3) = s.seed.sB;
    jacob(0, 4) = s.seed.cB;
    jacob(1, 3) = -s.seed.cB;
    jacob(1, 4) = s.seed.sB;
    jacob(2, 5) = s.seed.ds;
    jacob(3, 3) = s.seed.cos;
    jacob(3, 4) = s.seed.sin;
    jacob(4, 3) = -s.seed.sin;
    jacob(4, 4) = s.seed.cos;

    // -- d(r')/d(ds)
    Eigen::Vector<double, 6> df_ds;
    df_ds(0) = s.seed.cos * px0 + s.seed.sin * py0;
    df_ds(1) = -s.seed.sin * px0 + s.seed.cos * py0;
    df_ds(2) = pz0;
    df_ds(3) = -bq * s.seed.sin * px0 + bq * s.seed.cos * py0;
    df_ds(4) = -bq * s.seed.cos * px0 - bq * s.seed.sin * py0;
    df_ds(5) = 0.;

    Eigen::RowVector<double, 6> ds_dr = Eigen::Map<const Eigen::RowVector<double, 6>>(s.deriv.ds_dr.data());

    jacob.block<6, 6>(0, 0).noalias() += df_ds * ds_dr;

    // prepare Correlation Matrix //

    Eigen::RowVector<double, 6> ds_dr1 = Eigen::Map<const Eigen::RowVector<double, 6>>(s.deriv.ds_dr1.data());

    corr = df_ds * ds_dr1;

#if T2DS_DEBUG
    Logger::Debug(__FUNCTION__, "jacob  = {}", jacob);
    Logger::Debug(__FUNCTION__, "corr   = {}", corr);
#endif
}

// Transport particle.
// Input arguments: (packed in a single `Seeder::Seed` struct)
// - `ds`     -- transport parameter
// - `ds_dr`  -- partial derivatives of ds w.r.t. state parameters = d(ds1)/dr1
// - `ds_dr1` -- partial derivatives of current particle's ds w.r.t. other particle's state parameters = d(ds1)/dr2
// - `pca.xyz`, `pca.mom`, `theta`, `sin`, `cos`, `sB`, `cB`, `ds`
void Daughter::Transport(const Seeder::PCA& pca, const Eigen::Ref<const Eigen::Matrix<double, 8, 8>>& other_bt_cov) {

    // update state //

    fP(0) = pca.xyz[0];
    fP(1) = pca.xyz[1];
    fP(2) = pca.xyz[2];
    fP(3) = pca.mom[0];
    fP(4) = pca.mom[1];
    fP(5) = pca.mom[2];

    // update cov matrix //

    // -- with jacobian
    fC = jacob * fC.selfadjointView<Eigen::Lower>() * jacob.transpose();

    // -- with corr. matrix + other particle's cov matrix
    fC.block<6, 6>(0, 0).noalias() += corr * other_bt_cov.block<6, 6>(0, 0).selfadjointView<Eigen::Lower>() * corr.transpose();

#if T2DS_DEBUG
    Logger::Debug(__FUNCTION__, "fP = {}", fP);
    Logger::Debug(__FUNCTION__, "fC = {}", fC);
#endif
}

// == Mass Constraint == //

// Pin the state vector `p` (and propagate its covariance `c`) onto the mass shell defined by `mass`, by rescaling
//     p -> p/(1 - lambda)
//     E -> E/(1 + lambda)
// where lambda is chosen so that E'^2 - p'^2 = mass^2.
// Substituting gives the quartic:
//     f(lambda) = -m^2 * lambda^4 + a * lambda^2 + b * lambda + c
//     a = E^2 - p^2 + 2m^2
//     b = -2(E^2 + p^2)
//     c = E^2 - p^2 - m^2
// whose root lambda = 0 corresponds to a particle already on shell, since f(0) = c.
// The Jacobian d(p new)/d(p old) is stored in `j`, and the two rescaling factors are returned.
// If no usable lambda could be found, returns nullopt; leaving `p` and `c` untouched, and `j` as identity.
std::optional<MassScale> SetMassConstraint(Eigen::Vector<double, 8>& p, Eigen::Matrix<double, 8, 8>& c, Eigen::Matrix<double, 7, 7>& j, double mass) {

    j.setIdentity();

    // a negative energy cannot be rescaled onto a physical shell
    if (p(6) < 0.) return std::nullopt;

    const double energy2 = p(6) * p(6);
    const double mom2 = p.segment<3>(3).squaredNorm();
    const double mass2 = mass * mass;

    const double a = energy2 - mom2 + 2. * mass2;
    const double b = -2. * (energy2 + mom2);
    const double c0 = energy2 - mom2 - mass2;

    // solve for lambda //

    // -- seed with the smaller root of the quadratic part, i.e. dropping the -m^2*lambda^4 term.
    //    (energy2 + mom2) == -b/2 and d == (b^2 - 4ac)/4, so the textbook form is (energy2 + mom2 - sqrt(d))/a.
    //    That subtracts two nearly-equal numbers exactly in the near-on-shell case that dominates here, so use
    //    the conjugate instead: the roots multiply to c0/a, hence lambda_- == c0/((energy2 + mom2) + sqrt(d)).
    //    It also degrades gracefully to the linear root -c0/b as a -> 0, so no separate fallback is needed.
    const double d = 4. * energy2 * mom2 - mass2 * (energy2 - mom2 - 2. * mass2);
    const double q = energy2 + mom2 + (d > 0. ? std::sqrt(d) : 0.);

    double lambda = 0.;
    if (q > MassConstraint_MinDenom) lambda = c0 / q;

    // -- refine by Newton
    for (int i = 0; i < MassConstraint_MaxIter; ++i) {
        const double lambda2 = lambda * lambda;
        const double f = -mass2 * lambda2 * lambda2 + a * lambda2 + b * lambda + c0;
        const double df = -4. * mass2 * lambda2 * lambda + 2. * a * lambda + b;
        if (std::abs(df) < MassConstraint_MinDenom) break;
        const double delta = f / df;
        lambda -= delta;
        if (std::abs(delta) < MassConstraint_Tolerance) break;
    }

    // -- protection: lambda = +-1 would blow up the rescaling; leave caller's state untouched
    if (!std::isfinite(lambda) || std::abs(1. - lambda) < MassConstraint_MinDenom || std::abs(1. + lambda) < MassConstraint_MinDenom) {
        return std::nullopt;
    }

    const double lpi = 1. / (1. + lambda);
    const double lmi = 1. / (1. - lambda);
    const double lp2i = lpi * lpi;
    const double lm2i = lmi * lmi;

    // prepare Jacobian Matrix //

    // -- d(lambda)/d(px,py,pz,E), by implicit differentiation of f
    const double lambda2 = lambda * lambda;
    const double dfl = -4. * mass2 * lambda2 * lambda + 2. * a * lambda + b;

    // protection
    if (std::abs(dfl) < MassConstraint_MinDenom) return std::nullopt;

    Eigen::Vector<double, 4> dfx;
    dfx(0) = -2. * (1. + lambda) * (1. + lambda) * p(3);
    dfx(1) = -2. * (1. + lambda) * (1. + lambda) * p(4);
    dfx(2) = -2. * (1. + lambda) * (1. + lambda) * p(5);
    dfx(3) = 2. * (1. - lambda) * (1. - lambda) * p(6);

    Eigen::Vector<double, 4> dlx = -dfx / dfl;

    // -- d(px',py',pz',E')/d(lambda)
    Eigen::Vector<double, 4> dxx;
    dxx(0) = p(3) * lm2i;
    dxx(1) = p(4) * lm2i;
    dxx(2) = p(5) * lm2i;
    dxx(3) = -p(6) * lp2i;

    // -- position rows stay untouched, which is what keeps the position block of `c` (and hence the later
    //    D/K2/A/M correction in GetUpdatedMC) valid
    j.block<4, 4>(3, 3).noalias() = dxx * dlx.transpose();
    j.block<3, 3>(3, 3).diagonal().array() += lmi;
    j(6, 6) += lpi;

    // apply //

    // row/col 7 (S) is left alone, matching the original
    // NOTE: that leaves S's cross-covariances with the rescaled momenta stale. Harmless as long as nothing reads
    //       them -- `fP(7)` is never even filled -- but it is a real inconsistency, not just an omission.
    c.topLeftCorner<7, 7>() = (j * c.topLeftCorner<7, 7>() * j.transpose()).eval();

    p.segment<3>(3) *= lmi;
    p(6) *= lpi;

    return MassScale{lmi, lpi};
}

// Pin the fitted mother `part` onto the mass shell defined by `mass`, updating chi2, NDF (+= 1) and the mass bookkeeping.
// It always fires, when enabled. Returns the rescaling that was applied, so the caller can pass it on to the daughters.
std::optional<MassScale> SetNonlinearMassConstraint(Particle& part, double mass) {

    // h = d(m^2)/d(px,py,pz,E)
    Eigen::Vector<double, 4> h;
    h << -2. * part.Px(), -2. * part.Py(), -2. * part.Pz(), 2. * part.E();

    const double var_m2 = h.transpose() * part.fC.block<4, 4>(3, 3) * h;

    // -- the mass error is already 0, so the particle can't be constrained (the original doesn't guard this one,
    //    but its linearised counterpart does)
    if (var_m2 < MassConstraint_MinVariance) return std::nullopt;

    const double residual = part.E() * part.E() - part.SquaredMomentum() - mass * mass;

    // -- `j` is discarded: there's no cross-covariance left to rotate at this point
    Eigen::Matrix<double, 7, 7> j;
    const std::optional<MassScale> scale = SetMassConstraint(part.fP, part.fC, j, mass);
    if (!scale) return std::nullopt;

    // -- both bail-outs above happen before anything is mutated, so a failure leaves `part` untouched

    part.fChi2 += residual * residual / var_m2;
    part.fNDF += 1;
    part.fMassHypo = mass;
    part.fSumDaughterMass = mass;

    return scale;
}

// == Production Vertex == //

// Add `vtx` as a point-like measurement to `part`.
// The returned copy has its parameters *at* the production vertex, and the eighth component of its state vector filled
// with the decay length to momentum ratio (s = l/p) back to `part`'s decay vertex, covariances included.
// The constraint's own contribution is `out.Chi2() - part.Chi2()`, with 2 d.o.f.
KF::Particle SetProductionVertex(const KF::Particle& part, const KF::Vertex& vtx, double bz) {

    auto SeedToPoint = [](const KF::Particle& in_part, const std::array<double, 3>& xyz, double in_bz) -> Seeder::Result {
        if (in_part.fQ == 0 || std::abs(in_bz) < Common::AbsAlmostZero) {
            return Seeder::LineVertex::FullPCA(in_part.X(), in_part.Y(), in_part.Z(), in_part.Px(), in_part.Py(), in_part.Pz(), xyz);
        }
        return Seeder::HelixVertex::FullPCA(in_part.X(), in_part.Y(), in_part.Z(), in_part.Px(), in_part.Py(), in_part.Pz(), in_part.Charge(), xyz,
                                            in_bz);
    };

    // remember where the particle decayed //

    const Eigen::Vector<double, 3> decay_xyz = part.fP.head<3>();
    const Eigen::Matrix<double, 3, 3> decay_cov = part.fC.topLeftCorner<3, 3>();

    // transport the particle to the vertex //

    const Seeder::Result s = SeedToPoint(part, vtx.GetXYZ(), bz);

    KF::Daughter kf(part);
    kf.PrepareJacobAndCorr(s, bz);

    kf.fP.head<3>() = Eigen::Map<const Eigen::Vector<double, 3>>(s.seed.pca.xyz.data());
    kf.fP.segment<3>(3) = Eigen::Map<const Eigen::Vector<double, 3>>(s.seed.pca.mom.data());

    // -- cov(transported state, decay point) = F * C0[0:6, 0:3]. Kept for the var(S) cross term below; it has to
    //    be taken before `fC` is overwritten, and the vertex contributes nothing here (it is assumed independent
    //    of the candidate -- see the note on `fC(7,7)`).
    const Eigen::Matrix<double, 6, 3> mA = kf.jacob.topLeftCorner<6, 6>() * part.fC.block<6, 3>(0, 0);

    kf.fC = kf.jacob * kf.fC.selfadjointView<Eigen::Lower>() * kf.jacob.transpose();

    // -- correlation between the transported particle and the vertex: D[i][j] = sum_k V[j,k] * F1[i,k]
    // -- fold in the vertex's own error, F1 * V * F1^T. `corr`'s columns 3-5 are zero (`ds_dr1` only has position
    //    components), so `V` only ever enters through its 3x3 block.
    //    NOTE: `Daughter::Transport` is deliberately not reused here: it adds the whole 6x6, whereas the original
    //          writes back `fC[0..5]` only, i.e. the position block. For a neutral mother the two are identical
    //          anyway -- `bq = 0` makes rows 3-5 of `corr` vanish as well.
    const Eigen::Matrix<double, 3, 3> mD = kf.corr.topLeftCorner<3, 3>() * vtx.cov;
    kf.fC.topLeftCorner<3, 3>() += mD * kf.corr.topLeftCorner<3, 3>().transpose();

    // update with the vertex as a position measurement //

    Eigen::Matrix<double, 3, 3> mS = kf.fC.topLeftCorner<3, 3>() + vtx.cov;
    auto ldlt = mS.selfadjointView<Eigen::Lower>().ldlt();

    // -- residual = measured - estimated
    const Eigen::Vector<double, 3> res = vtx.xyz - kf.fP.head<3>();

    const Eigen::Matrix<double, 7, 3> mCHt = kf.fC.block<7, 3>(0, 0);
    const Eigen::Matrix<double, 7, 3> mK = ldlt.solve(mCHt.transpose()).transpose();

    // -- both `res` and `mS` are the pre-update ones, as in the original
    const double dchi2 = res.dot(ldlt.solve(res));

    // -- cov(updated state, decay point) = (1 - K * H) * mA, with H picking the position block
    const Eigen::Matrix<double, 6, 3> mG = mA - mK.topRows<6>() * mA.topRows<3>();

    kf.fP.head<7>().noalias() += mK * res;
    kf.fC.topLeftCorner<7, 7>().noalias() -= mK * mCHt.transpose();

    // correlation correction //

    const Eigen::Matrix<double, 3, 3> K = mK.topRows<3>();
    const Eigen::Matrix<double, 3, 3> K2 = Eigen::Matrix<double, 3, 3>::Identity() - K.transpose();
    // -- NOTE: transposed w.r.t. `GetUpdated`'s `A = D * K2`. Both mirror their own original, and the two `D`s carry
    //          opposite index orders: `AddDaughterWithEnergyFit` has `D[i][k] * K2[k][j]`, `SetProductionVertex` has
    //          `D[k][i] * K2[k][j]`.
    const Eigen::Matrix<double, 3, 3> A = mD.transpose() * K2;
    const Eigen::Matrix<double, 3, 3> M = K * A;

    kf.fC.topLeftCorner<3, 3>() += M + M.transpose();

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
    //    IMPORTANT: this still assumes the vertex is independent of the candidate. If the grand-daughter tracks
    //               entered the vertex fit upstream, that correlation is of the same order and nothing in here
    //               can recover it.
    const Eigen::Vector<double, 6> cds = kf.fC.topLeftCorner<6, 6>().selfadjointView<Eigen::Lower>() * ds_dr - mG * ds_dr.head<3>();

    kf.fC.block<1, 6>(7, 0) = cds.transpose();
    kf.fC.block<6, 1>(0, 7) = cds;                                     // NOTE: `fC` is full symmetric in this port, so both halves are written
    kf.fC(7, 7) = ds_dr.dot(cds)                                       //  j' * C * j - j' * G * j3
                  - ds_dr.head<3>().dot(mG.transpose() * ds_dr)        // -j3' * G' * j
                  + ds_dr.head<3>().dot(decay_cov * ds_dr.head<3>());  // +j3' * V_decay * j3

    // -- NOTE: `fC(7,6)`, the S-E covariance, is left alone. The original doesn't fill it either (its loop stops at
    //          `fC[33]`), so it stays at the 0 the transport produced.

    KF::Particle out = kf;

    out.fChi2 += dchi2;
    out.fNDF += 2;

#if T2DS_DEBUG
    Logger::Debug(__FUNCTION__, "mS     = {}", mS);
    Logger::Debug(__FUNCTION__, "res    = {}", res);
    Logger::Debug(__FUNCTION__, "mK     = {}", mK);
    Logger::Debug(__FUNCTION__, "mD     = {}", mD);
    Logger::Debug(__FUNCTION__, "mG     = {}", mG);
    Logger::Debug(__FUNCTION__, "M      = {}", M);
    Logger::Debug(__FUNCTION__, "dchi2  = {:13.6e}", dchi2);
    Logger::Debug(__FUNCTION__, "out.fP = {}", out.fP);
    Logger::Debug(__FUNCTION__, "out.fC = {}", out.fC);
#endif

    return out;
}

// == Main Fitting Methods == //

KF::FitResult GetUpdated(const KF::Daughter& kf_1, const KF::Daughter& kf_2) {

    KF::FitResult res;
    KF::Particle& out = res.mother;
    out.fC = kf_1.fC;
    out.fP = kf_1.fP;

    // sum of position covariances //

    Eigen::Matrix<double, 3, 3> mS = kf_1.fC.block<3, 3>(0, 0) + kf_2.fC.block<3, 3>(0, 0);
    auto ldlt = mS.selfadjointView<Eigen::Lower>().ldlt();

    // residual = measured - estimated //

    Eigen::Vector<double, 3> zeta = kf_2.fP.head<3>() - kf_1.fP.head<3>();

    // update chi2 //

    const Eigen::Vector<double, 3> mSz = ldlt.solve(zeta);
    out.fChi2 = zeta.dot(mSz);

    // correlation between state and position measurement //

    Eigen::Matrix<double, 7, 3> mCHt = kf_1.fC.block<7, 3>(0, 0);
    mCHt.block<4, 3>(3, 0) -= kf_2.fC.block<4, 3>(3, 0);

    // Kalman gain //

    Eigen::Matrix<double, 7, 3> mK = ldlt.solve(mCHt.transpose()).transpose();

    // add 4-momentum //

    out.fP.segment<4>(3).noalias() += kf_2.fP.segment<4>(3);
    out.fC.block<4, 4>(3, 3).noalias() += kf_2.fC.block<4, 4>(3, 3);

    // state update: P += K * zeta //

    out.fP.head<7>().noalias() += mK * zeta;

    // covariance update: C -= K * CHt' //

    out.fC.block<7, 7>(0, 0).noalias() -= mK * mCHt.transpose();

    // recover the individual daughter 4-momenta //

    res.dau_2.noalias() = kf_2.fP.segment<4>(3) - kf_2.fC.block<4, 3>(3, 0) * mSz;
    res.dau_1 = out.fP.segment<4>(3) - res.dau_2;

    // correlation correction //

    Eigen::Matrix<double, 3, 3> F3C1F1T =
        kf_2.corr.block<3, 6>(0, 0) * kf_1.bt_cov.block<6, 6>(0, 0).selfadjointView<Eigen::Lower>() * kf_1.jacob.block<3, 6>(0, 0).transpose();
    Eigen::Matrix<double, 3, 3> F4C2F2T =
        kf_2.jacob.block<3, 6>(0, 0) * kf_2.bt_cov.block<6, 6>(0, 0).selfadjointView<Eigen::Lower>() * kf_1.corr.block<3, 6>(0, 0).transpose();
    Eigen::Matrix<double, 3, 3> D = F3C1F1T + F4C2F2T;

    Eigen::Matrix<double, 3, 3> K = mK.topRows<3>();
    Eigen::Matrix<double, 3, 3> K2 = Eigen::Matrix<double, 3, 3>::Identity() - K.transpose();
    Eigen::Matrix<double, 3, 3> A = D * K2;
    Eigen::Matrix<double, 3, 3> M = K * A;

    out.fC.block<3, 3>(0, 0).noalias() += M + M.transpose();

    // update charge, NDF and mass bookkeeping //

    out.fQ = kf_1.Charge() + kf_2.Charge();
    out.fNDF += 2;
    out.fMassHypo = std::nullopt;
    out.fSumDaughterMass = kf_1.fSumDaughterMass + kf_2.fSumDaughterMass;

#if T2DS_DEBUG
    Logger::Debug(__FUNCTION__, "mS     = {}", mS);
    Logger::Debug(__FUNCTION__, "zeta   = {}", zeta);
    Logger::Debug(__FUNCTION__, "mCHt   = {}", mCHt);
    Logger::Debug(__FUNCTION__, "mK     = {}", mK);
    Logger::Debug(__FUNCTION__, "mSz    = {}", mSz);
    Logger::Debug(__FUNCTION__, "D      = {}", D);
    Logger::Debug(__FUNCTION__, "K      = {}", K);
    Logger::Debug(__FUNCTION__, "K2     = {}", K2);
    Logger::Debug(__FUNCTION__, "A      = {}", A);
    Logger::Debug(__FUNCTION__, "M      = {}", M);
    Logger::Debug(__FUNCTION__, "out.fP = {}", out.fP);
    Logger::Debug(__FUNCTION__, "out.fC = {}", out.fC);
#endif

    return res;
}

// Same vertex measurement as `GetUpdated()`, but the 4-momenta are combined only after each daughter has been pinned to
// its own mass shell; which guarantees mass(mother) >= sum(mass(daughters)).
// This forces a different assembly order:
// `GetUpdated()` folds the momentum sum into the Kalman update via the "CHt - D'" shortcut, whereas here
// the mother and the daughter measurement are updated separately, their cross-covariance `mDf` is tracked explicitly,
// and only then are the 4-momenta added.
KF::FitResult GetUpdatedMC(const KF::Daughter& kf_1, const KF::Daughter& kf_2) {

    KF::FitResult res;
    KF::Particle& out = res.mother;
    out.fC = kf_1.fC;
    out.fP = kf_1.fP;

    // the daughter measurement, mutated in place below //

    Eigen::Vector<double, 8> m = kf_2.fP;
    Eigen::Matrix<double, 8, 8> mV = kf_2.fC;

    // sum of position covariances //

    Eigen::Matrix<double, 3, 3> mS = kf_1.fC.block<3, 3>(0, 0) + kf_2.fC.block<3, 3>(0, 0);
    auto ldlt = mS.selfadjointView<Eigen::Lower>().ldlt();

    // residual = measured - estimated //

    Eigen::Vector<double, 3> zeta = m.head<3>() - kf_1.fP.head<3>();

    // update chi2 //

    const Eigen::Vector<double, 3> mSz = ldlt.solve(zeta);
    out.fChi2 = zeta.dot(mSz);

    // correlations with the position measurement //

    // -- different from `GetUpdated()`, the daughter is updated later
    Eigen::Matrix<double, 7, 3> mCHt = kf_1.fC.block<7, 3>(0, 0);
    Eigen::Matrix<double, 7, 3> mVHt = kf_2.fC.block<7, 3>(0, 0);

    // Kalman gains, one per particle //

    Eigen::Matrix<double, 7, 3> mK = ldlt.solve(mCHt.transpose()).transpose();
    Eigen::Matrix<double, 7, 3> mKm = ldlt.solve(mVHt.transpose()).transpose();

    // cross-covariance between the updated daughter and the updated mother //

    Eigen::Matrix<double, 7, 7> mDf = mKm * mCHt.transpose();

    // state updates: P += K * zeta, m -= Km * zeta //

    out.fP.head<7>().noalias() += mK * zeta;
    m.head<7>().noalias() -= mKm * zeta;

    // covariance updates //

    out.fC.topLeftCorner<7, 7>().noalias() -= mK * mCHt.transpose();
    mV.topLeftCorner<7, 7>().noalias() -= mKm * mVHt.transpose();

    // pin both particles to their mass shells //

    Eigen::Matrix<double, 7, 7> mJ1 = Eigen::Matrix<double, 7, 7>::Identity();
    Eigen::Matrix<double, 7, 7> mJ2 = Eigen::Matrix<double, 7, 7>::Identity();

    ConstrainToMassShell(out.fP, out.fC, mJ1, kf_1.fMassHypo, kf_1.fSumDaughterMass);
    ConstrainToMassShell(m, mV, mJ2, kf_2.fMassHypo, kf_2.fSumDaughterMass);

    mDf = (mJ2 * mDf * mJ1.transpose()).eval();

    // keep the two shell-pinned 4-momenta -- these are exactly what gets summed just below //

    res.dau_1 = out.fP.segment<4>(3);
    res.dau_2 = m.segment<4>(3);

    // add the daughter 4-momentum to the mother's //

    out.fP.segment<4>(3).noalias() += m.segment<4>(3);
    out.fC.block<4, 4>(3, 3).noalias() += mV.block<4, 4>(3, 3);

    // fold in cross-covariance //

    out.fC.block<4, 3>(3, 0).noalias() += mDf.block<4, 3>(3, 0);
    out.fC.block<4, 4>(3, 3).noalias() += mDf.block<4, 4>(3, 3) + mDf.block<4, 4>(3, 3).transpose();

    // correlation correction //

    Eigen::Matrix<double, 3, 3> F3C1F1T =
        kf_2.corr.block<3, 6>(0, 0) * kf_1.bt_cov.block<6, 6>(0, 0).selfadjointView<Eigen::Lower>() * kf_1.jacob.block<3, 6>(0, 0).transpose();
    Eigen::Matrix<double, 3, 3> F4C2F2T =
        kf_2.jacob.block<3, 6>(0, 0) * kf_2.bt_cov.block<6, 6>(0, 0).selfadjointView<Eigen::Lower>() * kf_1.corr.block<3, 6>(0, 0).transpose();
    Eigen::Matrix<double, 3, 3> D = F3C1F1T + F4C2F2T;

    Eigen::Matrix<double, 3, 3> K = mK.topRows<3>();
    Eigen::Matrix<double, 3, 3> K2 = Eigen::Matrix<double, 3, 3>::Identity() - K.transpose();
    Eigen::Matrix<double, 3, 3> A = D * K2;
    Eigen::Matrix<double, 3, 3> M = K * A;

    out.fC.block<3, 3>(0, 0).noalias() += M + M.transpose();

    // restore symmetry //

    // -- force `fC` stay a full symmetric matrix
    Eigen::Matrix<double, 7, 7> sym = out.fC.topLeftCorner<7, 7>();
    out.fC.topLeftCorner<7, 7>() = sym.selfadjointView<Eigen::Lower>();

    // update charge, NDF and mass bookkeeping //

    out.fQ = kf_1.Charge() + kf_2.Charge();
    out.fNDF += 2;
    out.fMassHypo = std::nullopt;
    out.fSumDaughterMass = kf_1.fSumDaughterMass + kf_2.fSumDaughterMass;

#if T2DS_DEBUG
    Logger::Debug(__FUNCTION__, "mS     = {}", mS);
    Logger::Debug(__FUNCTION__, "zeta   = {}", zeta);
    Logger::Debug(__FUNCTION__, "mCHt   = {}", mCHt);
    Logger::Debug(__FUNCTION__, "mVHt   = {}", mVHt);
    Logger::Debug(__FUNCTION__, "mK     = {}", mK);
    Logger::Debug(__FUNCTION__, "mSz    = {}", mSz);
    Logger::Debug(__FUNCTION__, "mKm    = {}", mKm);
    Logger::Debug(__FUNCTION__, "mJ1    = {}", mJ1);
    Logger::Debug(__FUNCTION__, "mJ2    = {}", mJ2);
    Logger::Debug(__FUNCTION__, "mDf    = {}", mDf);
    Logger::Debug(__FUNCTION__, "D      = {}", D);
    Logger::Debug(__FUNCTION__, "M      = {}", M);
    Logger::Debug(__FUNCTION__, "out.fP = {}", out.fP);
    Logger::Debug(__FUNCTION__, "out.fC = {}", out.fC);
#endif

    return res;
}

KF::FitResult FitVertex(const KF::Particle& part_1, const KF::Particle& part_2, const Seeder::Result& s_1, const Seeder::Result& s_2, double bz,
                        const FitPolicy& policy) {

    KF::Daughter kf_1(part_1);
    KF::Daughter kf_2(part_2);

    kf_1.PrepareJacobAndCorr(s_1, bz);
    kf_2.PrepareJacobAndCorr(s_2, bz);

    kf_1.Transport(s_1.seed.pca, kf_2.bt_cov);
    kf_2.Transport(s_2.seed.pca, kf_1.bt_cov);

    // constraint daughters' masses
    KF::FitResult res = policy.pin_daughters ? GetUpdatedMC(kf_1, kf_2) : GetUpdated(kf_1, kf_2);

    // constraint mother's mass
    // -- the pin rescales the mother's p by 1/(1-lambda) and its E by 1/(1+lambda); handing the daughters the
    //    very same two factors rescales their sum by exactly as much, so the tree stays closed
    if (policy.mother_mass) {
        if (const auto scale = SetNonlinearMassConstraint(res.mother, *policy.mother_mass)) res.RescaleDaughters(*scale);
    }

    // constraint production vertex
    if (policy.prod_vertex) res.at_pv = SetProductionVertex(res.mother, *policy.prod_vertex, bz);

    return res;
}

}  // namespace T2DS::KF
