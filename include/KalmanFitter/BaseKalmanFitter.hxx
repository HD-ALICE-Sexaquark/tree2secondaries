#pragma once

#include <Math/Point3D.h>
#include <Math/Vector3D.h>

#include <Eigen/Eigen>

#include "KalmanFitter/KalmanFitterParticle.hxx"
#include "Seeder/BaseSeeder.hxx"
#include "View/Reconstructed/ViewRecTrack.hxx"
#include "View/Reconstructed/ViewRecV0.hxx"

namespace Tree2Secondaries::KalmanFitter {

// Daughter Struct //

struct Daughter : Particle {
    explicit Daughter(const Particle& kf) : Particle{kf}, bt_cov{kf.fC} {}

    void PrepareJacobAndCorr(const Seeder::Result& s, double bz = 0.);
    void Transport(const Seeder::PCA& pca, const Eigen::Ref<const Eigen::Matrix<double, 8, 8>>& other_bt_cov);

    Eigen::Matrix<double, 8, 8> bt_cov;  // covariance before transport
    Eigen::Matrix<double, 8, 8> jacob;   // NOTE: uninitialized on purpose
    Eigen::Matrix<double, 6, 6> corr;    // NOTE: uninitialized on purpose
};

// Main Fitting Method //

Particle GetUpdated(const Daughter& kf_1, const Daughter& kf_2);

Particle FitVertex(const KalmanFitter::Particle& part_1, const KalmanFitter::Particle& part_2, const Seeder::Result& s_1, const Seeder::Result& s_2,
                   double bz = 0.);

// Inline Methods //

inline Particle FitVertex(const View::Rec::Track& track_1, const View::Rec::Track& track_2, double mass_1, double mass_2, const Seeder::Result& s_1,
                          const Seeder::Result& s_2, double bz) {
    return FitVertex(Particle::FromTrack(track_1, mass_1), Particle::FromTrack(track_2, mass_2), s_1, s_2, bz);
}

inline Particle FitVertex(const View::Rec::Track& track, const View::Rec::V0& v0, double mass_track, const Seeder::Result& s_track,
                          const Seeder::Result& s_v0, double bz) {
    return FitVertex(Particle::FromTrack(track, mass_track), Particle::FromV0(v0), s_track, s_v0, bz);
}

inline Particle FitVertex(const View::Rec::V0& v0_1, const View::Rec::V0& v0_2, const Seeder::Result& s_1, const Seeder::Result& s_2) {
    return FitVertex(Particle::FromV0(v0_1), Particle::FromV0(v0_2), s_1, s_2);
}

}  // namespace Tree2Secondaries::KalmanFitter
