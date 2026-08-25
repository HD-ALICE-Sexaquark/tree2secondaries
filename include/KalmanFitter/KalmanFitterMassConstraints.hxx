#pragma once

#include <optional>

#include <Eigen/Core>

#include "KalmanFitter/KalmanFitterParticle.hxx"
#include "KalmanFitter/KalmanFitterUtils.hxx"

// ## Mass Constraints ## //

namespace T2DS::KF::Internal {

// Pick the mass shell a daughter should be pinned to, and pin it.
// - an exact hypothesis always wins;
// - otherwise the constraint acts as a clamp, and only fires once the mass has fallen below the physical minimum.
// Return the pin's d(px,py,pz,E)'/d(px,py,pz,E). The position rows of the full 7x7 Jacobian are the identity and its
// off-diagonal blocks are zero, so only that corner is handed back.
// Return `std::nullopt` when no constraint applies and when one was attempted but couldn't be solved, leaving `p` and `c`
// untouched.
[[nodiscard]] std::optional<Eigen::Matrix<double, 4, 4>> PinDaughterToMassShell(Eigen::Vector<double, 8>& p, Eigen::Matrix<double, 8, 8>& c,
                                                                                const std::optional<double>& mass_hypo, double sum_daughter_mass);

// Pin the fitted mother `part` onto the mass shell defined by `mass`, updating chi2, NDF (+= 1) and the mass bookkeeping.
// It always fires, when enabled.
// Return the rescaling that was applied, so the caller can pass it on to the daughters; `std::nullopt` means nothing was touched.
[[nodiscard]] std::optional<MassScale> PinMotherMass(Particle& part, double mass);

}  // namespace T2DS::KF::Internal
