#pragma once

#include <cmath>

#include <Math/Point3Dfwd.h>
#include <Math/Vector3Dfwd.h>
#include <optional>

#include "common/Math.hpp"
#include "common/POD_Track.hpp"

#include "KalmanFitter/KalmanFitterParticle.hxx"
#include "Seeder/BaseSeeder.hxx"

namespace R2DS::KF {

struct V0 : Particle {

    V0() = delete;
    V0(const Particle& fit, const Seeder::PCA& pca_neg, const Seeder::PCA& pca_pos, const POD::Track& neg, const POD::Track& pos)
        : Particle{fit},  //
          Neg_PCAwrtV0{pca_neg},
          Pos_PCAwrtV0{pca_pos},
          Neg{&neg},
          Pos{&pos} {}

    // Member functions //

    [[nodiscard]] double SquaredDCA_Daughters() const { return Common::Math::SquaredDistance(Neg_PCAwrtV0.xyz, Pos_PCAwrtV0.xyz); }
    [[nodiscard]] double SquaredDCA_Neg_V0() const { return Common::Math::SquaredDistance(Neg_PCAwrtV0.xyz, GetXYZ()); }
    [[nodiscard]] double SquaredDCA_Pos_V0() const { return Common::Math::SquaredDistance(Pos_PCAwrtV0.xyz, GetXYZ()); }

    [[nodiscard]] double DCA_Daughters() const { return std::sqrt(SquaredDCA_Daughters()); }
    [[nodiscard]] double DCA_Neg_V0() const { return std::sqrt(SquaredDCA_Neg_V0()); }
    [[nodiscard]] double DCA_Pos_V0() const { return std::sqrt(SquaredDCA_Pos_V0()); }

    [[nodiscard]] double DCA_Vertex(double x, double y, double z) const {
        return Common::Math::FastDCA_LineVertex(GetPxPyPz_AsROOT(), GetXYZ_AsROOT(), {x, y, z});
    }
    [[nodiscard]] double CPA_Vertex(double x, double y, double z) const {
        return Common::Math::CosinePointingAngle(GetPxPyPz_AsROOT(), GetXYZ_AsROOT(), {x, y, z});
    }

    [[nodiscard]] double ArmenterosQt() const {
        return Common::Math::ArmenterosQt(Px(), Py(), Pz(), Neg_PCAwrtV0.Px(), Neg_PCAwrtV0.Py(), Neg_PCAwrtV0.Pz());
    }
    [[nodiscard]] std::optional<double> ArmenterosAlpha() const {
        return Common::Math::ArmenterosAlpha(Px(), Py(), Pz(),                                         //
                                             Neg_PCAwrtV0.Px(), Neg_PCAwrtV0.Py(), Neg_PCAwrtV0.Pz(),  //
                                             Pos_PCAwrtV0.Px(), Pos_PCAwrtV0.Py(), Pos_PCAwrtV0.Pz());
    }

    // Member variables //

    Seeder::PCA Neg_PCAwrtV0;
    Seeder::PCA Pos_PCAwrtV0;
    const POD::Track* Neg;
    const POD::Track* Pos;
};

}  // namespace R2DS::KF
