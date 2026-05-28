#pragma once

#include <cstddef>

#include <Math/Point3D.h>
#include <Math/Vector3D.h>

#include <Eigen/Eigen>

#include "common/VC_TrackView.hpp"
#include "common/VC_V0View.hpp"

#include "KalmanFitter/KalmanFitterParticle.hxx"
#include "Seeder/BaseSeeder.hxx"

namespace R2DS::KalmanFitter {

static constexpr double Initial_Css = 1.;

static std::size_t IJ(std::size_t i, std::size_t j) { return (j <= i) ? i * (i + 1) / 2 + j : j * (j + 1) / 2 + i; }

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

Particle FitVertex(const KalmanFitter::Particle& part_1, const KalmanFitter::Particle& part_2, const Seeder::Result& s_1, const Seeder::Result& s_2,
                   double bz = 0.);

// Inline Methods //

inline Particle FitVertex(const Vector::TrackView& track_1, const Vector::TrackView& track_2, double mass_1, double mass_2, const Seeder::Result& s_1,
                          const Seeder::Result& s_2, double bz) {
    return FitVertex(Particle::FromTrack(track_1, mass_1), Particle::FromTrack(track_2, mass_2), s_1, s_2, bz);
}

inline Particle FitVertex(const Vector::TrackView& track, const Vector::V0View& v0, double mass_track, const Seeder::Result& s_track,
                          const Seeder::Result& s_v0, double bz) {
    return FitVertex(Particle::FromTrack(track, mass_track), Particle::FromV0(v0), s_track, s_v0, bz);
}

inline Particle FitVertex(const Vector::V0View& v0_1, const Vector::V0View& v0_2, const Seeder::Result& s_1, const Seeder::Result& s_2) {
    return FitVertex(Particle::FromV0(v0_1), Particle::FromV0(v0_2), s_1, s_2);
}

}  // namespace R2DS::KalmanFitter
