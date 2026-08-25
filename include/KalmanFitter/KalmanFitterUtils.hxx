#pragma once

#include <cstddef>

#include <Eigen/Core>

// ## Utilities ## //

namespace T2DS::KF {

// == Constants == //

namespace Constants {

inline constexpr double Initial_Css = 1.;
inline constexpr int Initial_NDF = -1;

}  // namespace Constants

// == Covariance Packing == //

// Index of the (i,j) entry of a symmetric matrix stored as a lower-triangular packed array, which is how every
// `POD::` type carries its covariance.
constexpr std::size_t IJ(std::size_t i, std::size_t j) { return (j <= i) ? i * (i + 1) / 2 + j : j * (j + 1) / 2 + i; }

// Unpack a `POD::` covariance into the top-left NxN corner of `out`. Fill both halves = store full symmetric.
template <unsigned int N, typename Derived, typename Packed>
void UnpackSym(Eigen::MatrixBase<Derived> &out, const Packed &packed) {
    for (unsigned int i = 0; i < N; ++i) {
        for (unsigned int j = 0; j <= i; ++j) {
            out(j, i) = out(i, j) = static_cast<double>(packed[IJ(i, j)]);
        }
    }
}

// == Mass Shell Rescaling == //

// The two factors a mass-shell pin applies: p -> p * p_scale, E -> E * e_scale. Both come out of the same
// lambda, so applying them to a set of 4-momenta rescales their sum by exactly the same amounts.
struct MassScale {
    double p_scale{1.};  // 1/(1 - lambda)
    double e_scale{1.};  // 1/(1 + lambda)
};

}  // namespace T2DS::KF
