#pragma once

#include "common/DB_Particles.hpp"
#include "common/POD_PreFoundLambda.hpp"
#include "common/POD_Track.hpp"
#include "common/POD_V0.hpp"

#include "KalmanFitter/KalmanFitterFitTypes.hxx"
#include "KalmanFitter/KalmanFitterParticle.hxx"
#include "Seeder/SeederTypes.hxx"

namespace T2DS::KF {

// == Component Factories == //

// Build one side of a fit out of a `POD::` type, its particle hypothesis and its seed.
// `composite_on_shell` asserts that a *composite* daughter arrives carrying an exact mass hypothesis. There is no `on_shell` for tracks.

[[nodiscard]] inline Component MakeComponent(const POD::Track& track, const DB::Particles::Definition& pid, const Seeder::Result& seed) {
    return {Particle::FromTrack(track, pid), seed};
}

[[nodiscard]] inline Component MakeComponent(const POD::V0& v0, const DB::Particles::Definition& pid, const Seeder::Result& seed,
                                             bool composite_on_shell) {
    return {Particle::FromV0(v0, pid, composite_on_shell), seed};
}

[[nodiscard]] inline Component MakeComponent(const POD::Extended::PreFoundLambda& lambda, const DB::Particles::Definition& pid,
                                             const Seeder::Result& seed, bool composite_on_shell) {
    return {Particle::FromPreFoundLambda(lambda, pid, composite_on_shell), seed};
}

// == Fit == //

// Fit a two-prong decay vertex out of two seeded particles, subject to a certain fit policy.
[[nodiscard]] FitResult FitVertex(const Component& p1, const Component& p2, double bz = 0., const FitPolicy& policy = {});

}  // namespace T2DS::KF
