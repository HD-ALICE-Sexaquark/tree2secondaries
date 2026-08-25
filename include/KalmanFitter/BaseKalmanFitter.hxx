#pragma once

#include <cstddef>
#include <optional>

#include <Math/Point3Dfwd.h>
#include <Math/Vector3Dfwd.h>

#include <Eigen/Eigen>

#include "common/DB_Particles.hpp"
#include "common/POD_PreFoundLambda.hpp"
#include "common/POD_Track.hpp"
#include "common/POD_V0.hpp"

#include "KalmanFitter/KalmanFitterParticle.hxx"
#include "Seeder/SeederTypes.hxx"

namespace T2DS::KF {

constexpr std::size_t IJ(std::size_t i, std::size_t j) { return (j <= i) ? i * (i + 1) / 2 + j : j * (j + 1) / 2 + i; }

// Mass Shell Rescaling //

// The two factors a mass-shell pin applies: p -> p * p_scale, E -> E * e_scale. Both come out of the same
// lambda, so applying them to a set of 4-momenta rescales their sum by exactly the same amounts.
struct MassScale {
    double p_scale{1.};  // 1/(1 - lambda)
    double e_scale{1.};  // 1/(1 + lambda)
};

// Fit Result //

struct FitResult {
    Particle mother;                                                   // at decay/secondary vertex; closed tree values
    std::optional<Particle> at_pv;                                     // available when pinned to production/primary vertex
    Eigen::Vector<double, 4> dau_1{Eigen::Vector<double, 4>::Zero()};  // (px, py, pz, E)
    Eigen::Vector<double, 4> dau_2{Eigen::Vector<double, 4>::Zero()};

    // Pass a mother pin's rescaling down to both daughters, keeping the 4-momentum sum exact.
    void RescaleDaughters(const MassScale& s) {
        dau_1.head<3>() *= s.p_scale;
        dau_1(3) *= s.e_scale;
        dau_2.head<3>() *= s.p_scale;
        dau_2(3) *= s.e_scale;
    }

    [[nodiscard]] double Dau1_Px() const { return dau_1(0); }
    [[nodiscard]] double Dau1_Py() const { return dau_1(1); }
    [[nodiscard]] double Dau1_Pz() const { return dau_1(2); }
    [[nodiscard]] double Dau1_E() const { return dau_1(3); }

    [[nodiscard]] double Dau2_Px() const { return dau_2(0); }
    [[nodiscard]] double Dau2_Py() const { return dau_2(1); }
    [[nodiscard]] double Dau2_Pz() const { return dau_2(2); }
    [[nodiscard]] double Dau2_E() const { return dau_2(3); }
};

// Daughter Struct //

struct Daughter : Particle {
    explicit Daughter(const Particle& kf) : Particle{kf}, bt_cov{kf.fC} {}

    void PrepareJacobAndCorr(const Seeder::Result& s, double bz = 0.);
    void Transport(const Seeder::PCA& pca, const Eigen::Ref<const Eigen::Matrix<double, 8, 8>>& other_bt_cov);

    Eigen::Matrix<double, 8, 8> bt_cov;  // covariance before transport
    Eigen::Matrix<double, 8, 8> jacob;   // NOTE: non-initialized on purpose
    Eigen::Matrix<double, 6, 6> corr;    // NOTE: non-initialized on purpose
};

// Fit Policy //

// Toggle fit constraints.
// - `pin_daughters` -- acts *during* the 4-momentum sum (this is KFParticle's `fConstructMethod == 2`): each daughter
//                      is put back on its mass shell first, which guarantees mass(mother) >= sum(mass(daughters)). Off, momentum and
//                      energy are treated as independent and the daughter masses may drift off shell -- KFParticle's method 0.
// - `daughters_already_pinned` -- composite daughters arrive carrying an exact mass hypothesis, which does
//                                 not survive the POD round-trip
//                                 NOTE: inert unless `pin_daughters` = true
// - `mother_mass` -- acts on the fitted mother *afterwards*, it always fires when enabled.
//                    IMPORTANT: this constraint makes `Mass()` a delta, however it pays for itself in both chi2 and NDF.
// - `prod_vertex` -- pins fit to a certain production vertex
struct FitPolicy {
    bool pin_daughters{false};
    bool daughters_already_pinned{false};
    std::optional<double> mother_mass;
    std::optional<Vertex> prod_vertex;
};

// Mass Constraint //

[[nodiscard]] std::optional<MassScale> SetMassConstraint(Eigen::Vector<double, 8>& p, Eigen::Matrix<double, 8, 8>& c, Eigen::Matrix<double, 7, 7>& j,
                                                         double mass);
[[nodiscard]] std::optional<MassScale> SetNonlinearMassConstraint(Particle& part, double mass);

// Pick the mass shell a daughter should be pinned to, and pin it.
// - an exact hypothesis always wins;
// - otherwise the constraint acts as a clamp, and only fires once the mass has fallen below the physical minimum.
// `j` is left as identity, when no constraint applies and when one was attempted but couldn't be solved,
// leaving `p` and `c` untouched
inline void ConstrainToMassShell(Eigen::Vector<double, 8>& p, Eigen::Matrix<double, 8, 8>& c, Eigen::Matrix<double, 7, 7>& j,
                                 const std::optional<double>& mass_hypo, double sum_daughter_mass) {
    j.setIdentity();

    if (mass_hypo) {
        static_cast<void>(SetMassConstraint(p, c, j, *mass_hypo));
        return;
    }

    const double m2 = p(6) * p(6) - p.segment<3>(3).squaredNorm();
    if (m2 < sum_daughter_mass * sum_daughter_mass) {
        static_cast<void>(SetMassConstraint(p, c, j, sum_daughter_mass));
    }
}

// Production Vertex //

[[nodiscard]] Particle SetProductionVertex(const Particle& part, const Vertex& vtx, double bz = 0.);

// Main Fitting Methods //

FitResult GetUpdated(const Daughter& kf_1, const Daughter& kf_2);
FitResult GetUpdatedMC(const Daughter& kf_1, const Daughter& kf_2);

FitResult FitVertex(const KF::Particle& part_1, const KF::Particle& part_2, const Seeder::Result& s_1, const Seeder::Result& s_2, double bz = 0.,
                    const FitPolicy& policy = {});

// Inline Methods //

inline FitResult FitVertex(const POD::Track& track_1, const POD::Track& track_2, const DB::Particles::Definition& pid_1,
                           const DB::Particles::Definition& pid_2, const Seeder::Result& s_1, const Seeder::Result& s_2, double bz,
                           const FitPolicy& policy = {}) {
    return FitVertex(Particle::FromTrack(track_1, pid_1), Particle::FromTrack(track_2, pid_2), s_1, s_2, bz, policy);
}

inline FitResult FitVertex(const POD::Track& track, const POD::V0& v0, const DB::Particles::Definition& pid_track,
                           const DB::Particles::Definition& pid_v0, const Seeder::Result& s_track, const Seeder::Result& s_v0, double bz,
                           const FitPolicy& policy = {}) {
    return FitVertex(Particle::FromTrack(track, pid_track), Particle::FromV0(v0, pid_v0, policy.daughters_already_pinned), s_track, s_v0, bz, policy);
}

inline FitResult FitVertex(const POD::V0& v0_1, const POD::V0& v0_2, const DB::Particles::Definition& pid_1, const DB::Particles::Definition& pid_2,
                           const Seeder::Result& s_1, const Seeder::Result& s_2, const FitPolicy& policy = {}) {
    return FitVertex(Particle::FromV0(v0_1, pid_1, policy.daughters_already_pinned), Particle::FromV0(v0_2, pid_2, policy.daughters_already_pinned),
                     s_1, s_2, 0., policy);
}

inline FitResult FitVertex(const POD::Extended::PreFoundLambda& l1, const POD::Extended::PreFoundLambda& l2, const DB::Particles::Definition& pid_1,
                           const DB::Particles::Definition& pid_2, const Seeder::Result& s1, const Seeder::Result& s2, const FitPolicy& policy = {}) {
    // `bz`=0 is safe, because `l1` and `l2` are neutral particles; equivalent to q=0 route for both particles
    return FitVertex(Particle::FromPreFoundLambda(l1, pid_1, policy.daughters_already_pinned),
                     Particle::FromPreFoundLambda(l2, pid_2, policy.daughters_already_pinned), s1, s2, 0., policy);
}

}  // namespace T2DS::KF
