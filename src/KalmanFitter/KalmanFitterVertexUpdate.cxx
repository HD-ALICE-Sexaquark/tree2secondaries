#include <optional>

#include <Eigen/Core>

#if T2DS_DEBUG
#include "App/Logger.hxx"
#include "App/Utilities.hxx"
#endif
#include "KalmanFitter/KalmanFitterAlgebra.hxx"
#include "KalmanFitter/KalmanFitterFitTypes.hxx"
#include "KalmanFitter/KalmanFitterMassConstraints.hxx"
#include "KalmanFitter/KalmanFitterTransport.hxx"

#include "KalmanFitter/KalmanFitterVertexUpdate.hxx"

// ## Vertex Update ## //

namespace T2DS::KF::Internal {

// NOTE: the mother starts out at `Constants::Initial_NDF` with a zero chi2, i.e. the first daughter's own chi2 is dropped.
//       That matches `original/AddDaughter.cxx`, whose first-daughter copy doesn't carry it over either.
FitResult VertexUpdate(const TransportedPair& pair) {

    const Particle& kf_1 = pair.first;
    const Particle& kf_2 = pair.second;

    FitResult res;
    Particle& out = res.mother;
    out.fC = kf_1.fC;
    out.fP = kf_1.fP;

    // sum of position covariances //

    const Eigen::Matrix<double, 3, 3> mS = kf_1.fC.topLeftCorner<3, 3>() + kf_2.fC.topLeftCorner<3, 3>();
    const Eigen::Matrix<double, 3, 3> mS_inv = Algebra::InvertSym3(mS);

    // residual = measured - estimated //

    const Eigen::Vector<double, 3> zeta = kf_2.fP.head<3>() - kf_1.fP.head<3>();

    // update chi2 //

    const Eigen::Vector<double, 3> mSz = mS_inv * zeta;
    out.fChi2 = zeta.dot(mSz);

    // correlation between state and position measurement //

    Eigen::Matrix<double, 7, 3> mCHt = kf_1.fC.block<7, 3>(0, 0);
    mCHt.block<4, 3>(3, 0) -= kf_2.fC.block<4, 3>(3, 0);

    // Kalman gain //

    // -- K = (S^-1 * CHt')' = CHt * S^-1, because S (and hence S^-1) is symmetric. Both transposes drop out and
    //    what's left is a plain gemm.
    const Eigen::Matrix<double, 7, 3> mK = mCHt * mS_inv;

    // add 4-momentum //

    out.fP.segment<4>(3).noalias() += kf_2.fP.segment<4>(3);
    out.fC.block<4, 4>(3, 3).noalias() += kf_2.fC.block<4, 4>(3, 3);

    // state update: P += K * zeta //

    out.fP.head<7>().noalias() += mK * zeta;

    // covariance update: C -= K * CHt' //

    out.fC.topLeftCorner<7, 7>().noalias() -= mK * mCHt.transpose();  // NOTE: not re-symmetrised on purpose

    // recover the individual daughter 4-momenta //

    res.dau_2.noalias() = kf_2.fP.segment<4>(3) - kf_2.fC.block<4, 3>(3, 0) * mSz;
    res.dau_1 = out.fP.segment<4>(3) - res.dau_2;

    // correlation correction //

    Algebra::ApplyCrossCorrection(out.fC.topLeftCorner<3, 3>(), pair.cross, mK.topRows<3>());

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
    Logger::Debug(__FUNCTION__, "out.fP = {}", out.fP);
    Logger::Debug(__FUNCTION__, "out.fC = {}", out.fC);
#endif

    return res;
}

// Same vertex measurement as `VertexUpdate()`, but the 4-momenta are combined only after each daughter has been pinned
// to its own mass shell; which guarantees mass(mother) >= sum(mass(daughters)).
// This forces a different assembly order:
// `VertexUpdate()` folds the momentum sum into the Kalman update via the "CHt - D'" shortcut, whereas here
// the mother and the daughter measurement are updated separately, their cross-covariance `mDf` is tracked explicitly,
// and only then are the 4-momenta added.
FitResult VertexUpdateMC(const TransportedPair& pair) {

    const Particle& kf_1 = pair.first;
    const Particle& kf_2 = pair.second;

    FitResult res;
    Particle& out = res.mother;
    out.fC = kf_1.fC;
    out.fP = kf_1.fP;

    // the daughter measurement, mutated in place below //

    Eigen::Vector<double, 8> m = kf_2.fP;
    Eigen::Matrix<double, 8, 8> mV = kf_2.fC;

    // sum of position covariances //

    const Eigen::Matrix<double, 3, 3> mS = kf_1.fC.topLeftCorner<3, 3>() + kf_2.fC.topLeftCorner<3, 3>();
    const Eigen::Matrix<double, 3, 3> mS_inv = Algebra::InvertSym3(mS);

    // residual = measured - estimated //

    const Eigen::Vector<double, 3> zeta = m.head<3>() - kf_1.fP.head<3>();

    // update chi2 //

    const Eigen::Vector<double, 3> mSz = mS_inv * zeta;
    out.fChi2 = zeta.dot(mSz);

    // correlations with the position measurement //

    // -- different from `VertexUpdate()`, the daughter is updated later
    const Eigen::Matrix<double, 7, 3> mCHt = kf_1.fC.block<7, 3>(0, 0);
    const Eigen::Matrix<double, 7, 3> mVHt = kf_2.fC.block<7, 3>(0, 0);

    // Kalman gains, one per particle //

    // -- see the note in `VertexUpdate()`: S is symmetric, so both transposes drop out
    const Eigen::Matrix<double, 7, 3> mK = mCHt * mS_inv;
    const Eigen::Matrix<double, 7, 3> mKm = mVHt * mS_inv;

    // cross-covariance between the updated daughter and the updated mother //

    // -- only its (px,py,pz,E) rows are ever folded back in, so the position rows are never formed
    Eigen::Matrix<double, 4, 7> mDf = mKm.bottomRows<4>() * mCHt.transpose();

    // state updates: P += K * zeta, m -= Km * zeta //

    out.fP.head<7>().noalias() += mK * zeta;
    m.head<7>().noalias() -= mKm * zeta;

    // covariance updates //

    out.fC.topLeftCorner<7, 7>().noalias() -= mK * mCHt.transpose();
    mV.topLeftCorner<7, 7>().noalias() -= mKm * mVHt.transpose();

    // pin both particles to their mass shells //

    const std::optional<Eigen::Matrix<double, 4, 4>> mJ1 = PinDaughterToMassShell(out.fP, out.fC, kf_1.fMassHypo, kf_1.fSumDaughterMass);
    const std::optional<Eigen::Matrix<double, 4, 4>> mJ2 = PinDaughterToMassShell(m, mV, kf_2.fMassHypo, kf_2.fSumDaughterMass);

    // -- rotate into the pinned frame, mDf -> J2 * mDf * J1'. Both Jacobians are the identity outside their
    //    (px,py,pz,E) block, so the position columns only ever pick up the left factor; and a pin that didn't fire
    //    *is* the identity, so its side is skipped rather than multiplied through.
    if (mJ2) {
        mDf.leftCols<3>() = (*mJ2 * mDf.leftCols<3>()).eval();
        mDf.rightCols<4>() = (*mJ2 * mDf.rightCols<4>()).eval();
    }
    if (mJ1) mDf.rightCols<4>() = (mDf.rightCols<4>() * mJ1->transpose()).eval();

    // keep the two shell-pinned 4-momenta -- these are exactly what gets summed just below //

    res.dau_1 = out.fP.segment<4>(3);
    res.dau_2 = m.segment<4>(3);

    // add the daughter 4-momentum to the mother's //

    out.fP.segment<4>(3).noalias() += m.segment<4>(3);
    out.fC.block<4, 4>(3, 3).noalias() += mV.block<4, 4>(3, 3);

    // fold in cross-covariance //

    out.fC.block<4, 3>(3, 0) += mDf.leftCols<3>();
    out.fC.block<4, 4>(3, 3) += mDf.rightCols<4>() + mDf.rightCols<4>().transpose();

    // correlation correction //

    Algebra::ApplyCrossCorrection(out.fC.topLeftCorner<3, 3>(), pair.cross, mK.topRows<3>());

    // restore symmetry //

    out.fC.block<3, 4>(0, 3) = out.fC.block<4, 3>(3, 0).transpose();  // NOTE: it's ok to just mirror this block

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
    const Eigen::Matrix<double, 4, 4> id_4 = Eigen::Matrix<double, 4, 4>::Identity();  // -- stands in for a pin that didn't fire
    Logger::Debug(__FUNCTION__, "mJ1    = {}", mJ1.value_or(id_4));
    Logger::Debug(__FUNCTION__, "mJ2    = {}", mJ2.value_or(id_4));
    Logger::Debug(__FUNCTION__, "mDf    = {}", mDf);
    Logger::Debug(__FUNCTION__, "out.fP = {}", out.fP);
    Logger::Debug(__FUNCTION__, "out.fC = {}", out.fC);
#endif

    return res;
}

}  // namespace T2DS::KF::Internal
