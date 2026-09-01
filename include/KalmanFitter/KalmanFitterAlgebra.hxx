#pragma once

#include <Eigen/Core>
#include <Eigen/LU>

#if T2DS_DEBUG
#include "Utils/Logger.hxx"
#include "Utils/Utilities.hxx"
#endif

namespace T2DS::KF::Internal::Algebra {

// Fold a position cross-covariance back into the updated position block:
//     C[0:3,0:3] += M + M', with M = K * cross * (1 - K')
// `cross` is whatever correlates the two things that were just merged -- the two daughters in the vertex updates,
// the candidate and the vertex in `AtProductionVertex`.
inline void ApplyCrossCorrection(Eigen::Ref<Eigen::Matrix<double, 3, 3>> c_pos, const Eigen::Matrix<double, 3, 3>& cross,
                                 const Eigen::Matrix<double, 3, 3>& k) {

    const Eigen::Matrix<double, 3, 3> mA = cross * (Eigen::Matrix<double, 3, 3>::Identity() - k.transpose());
    const Eigen::Matrix<double, 3, 3> mM = k * mA;

    c_pos += mM + mM.transpose();

#if T2DS_DEBUG
    Logger::Debug(__FUNCTION__, "cross = {}", cross);
    Logger::Debug(__FUNCTION__, "K     = {}", k);
    Logger::Debug(__FUNCTION__, "A     = {}", mA);
    Logger::Debug(__FUNCTION__, "M     = {}", mM);
#endif
}

// Invert a symmetric 3x3, forcing the result symmetric.
// Every fit below solves against the same matrix for a whole batch of right-hand sides at once (8, 15 and 8
// columns respectively), so one closed-form inverse spends a single division where an LDLT solve spends three per
// column. The input is always a sum of position covariances -- SPD, and better conditioned than either term -- so
// the cofactor form keeps far more precision than the float Cholesky inverse of the original.
// NOTE: a singular input yields inf/nan, unguarded. So does an LDLT solve, so this is no worse; but neither is it
//       a check, and nothing downstream performs one.
[[nodiscard]] inline Eigen::Matrix<double, 3, 3> InvertSym3(const Eigen::Matrix<double, 3, 3>& m) {

    Eigen::Matrix<double, 3, 3> out = m.inverse();

    // -- force symmetry
    out(0, 1) = out(1, 0);
    out(0, 2) = out(2, 0);
    out(1, 2) = out(2, 1);

    return out;
}

}  // namespace T2DS::KF::Internal::Algebra
