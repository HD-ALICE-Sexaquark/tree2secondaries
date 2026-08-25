#pragma once

#include <optional>

#include <Eigen/Core>

#include "KalmanFitter/KalmanFitterParticle.hxx"
#include "KalmanFitter/KalmanFitterUtils.hxx"
#include "KalmanFitter/KalmanFitterVertex.hxx"
#include "Seeder/SeederTypes.hxx"

namespace T2DS::KF {

// == Fit Component == //

struct Component {
    Particle part;
    Seeder::Result seed;
};

// == Fit Result == //

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

// == Fit Policy == //

// Toggle fit constraints.
// - `pin_daughters` -- acts *during* the 4-momentum sum (this is KFParticle's `fConstructMethod == 2`): each daughter
//                      is put back on its mass shell first, which guarantees mass(mother) >= sum(mass(daughters)). Off, momentum and
//                      energy are treated as independent and the daughter masses may drift off shell -- KFParticle's method 0.
// - `mother_mass` -- acts on the fitted mother *afterwards*, it always fires when enabled.
//                    IMPORTANT: this constraint makes `Mass()` a delta, however it pays for itself in both chi2 and NDF.
// - `prod_vertex` -- pins fit to a certain production vertex
struct FitPolicy {
    bool pin_daughters{false};
    std::optional<double> mother_mass;
    std::optional<Vertex> prod_vertex;
};

}  // namespace T2DS::KF
