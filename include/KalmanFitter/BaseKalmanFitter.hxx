#pragma once

#include <cstddef>

#include <Math/Point3Dfwd.h>
#include <Math/Vector3Dfwd.h>

#include <Eigen/Eigen>

#include "common/POD_PreFoundLambda.hpp"
#include "common/POD_Track.hpp"
#include "common/POD_V0.hpp"

#include "KalmanFitter/KalmanFitterParticle.hxx"
#include "Seeder/BaseSeeder.hxx"

namespace T2DS::KF {

static constexpr double Initial_Css = 1.;

constexpr std::size_t IJ(std::size_t i, std::size_t j) { return (j <= i) ? i * (i + 1) / 2 + j : j * (j + 1) / 2 + i; }

// Daughter Struct //

struct Daughter : Particle {
    explicit Daughter(const Particle& kf) : Particle{kf}, bt_cov{kf.fC} {}

    void PrepareJacobAndCorr(const Seeder::Result& s, double bz = 0.);
    void Transport(const Seeder::PCA& pca, const Eigen::Ref<const Eigen::Matrix<double, 8, 8>>& other_bt_cov);

    Eigen::Matrix<double, 8, 8> bt_cov;  // covariance before transport
    Eigen::Matrix<double, 8, 8> jacob;   // NOTE: non-initialized on purpose
    Eigen::Matrix<double, 6, 6> corr;    // NOTE: non-initialized on purpose
};

// Main Fitting Method //

Particle GetUpdated(const Daughter& kf_1, const Daughter& kf_2);

Particle FitVertex(const KF::Particle& part_1, const KF::Particle& part_2, const Seeder::Result& s_1, const Seeder::Result& s_2, double bz = 0.);

// Inline Methods //

inline Particle FitVertex(const POD::Track& track_1, const POD::Track& track_2, double mass_1, double mass_2, const Seeder::Result& s_1,
                          const Seeder::Result& s_2, double bz) {
    return FitVertex(Particle::FromTrack(track_1, mass_1), Particle::FromTrack(track_2, mass_2), s_1, s_2, bz);
}

inline Particle FitVertex(const POD::Track& track, const POD::V0& v0, double mass_track, const Seeder::Result& s_track, const Seeder::Result& s_v0,
                          double bz) {
    return FitVertex(Particle::FromTrack(track, mass_track), Particle::FromV0(v0), s_track, s_v0, bz);
}

inline Particle FitVertex(const POD::V0& v0_1, const POD::V0& v0_2, const Seeder::Result& s_1, const Seeder::Result& s_2) {
    return FitVertex(Particle::FromV0(v0_1), Particle::FromV0(v0_2), s_1, s_2);
}

inline Particle FitVertex(const POD::Extended::PreFoundLambda& l1, const POD::Extended::PreFoundLambda& l2, const Seeder::Result& s1,
                          const Seeder::Result& s2) {
    return FitVertex(Particle::FromPreFoundLambda(l1), Particle::FromPreFoundLambda(l2), s1, s2);
}

}  // namespace T2DS::KF
