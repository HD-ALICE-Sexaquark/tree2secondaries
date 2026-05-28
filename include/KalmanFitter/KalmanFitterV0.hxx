#pragma once

#include <cmath>

#include <Math/Point3D.h>
#include <Math/Vector3D.h>

#include "common/Constants.hpp"
#include "common/Math.hpp"
#include "common/VC_TrackView.hpp"

#include "KalmanFitter/KalmanFitterParticle.hxx"
#include "Seeder/BaseSeeder.hxx"

namespace R2DS::KalmanFitter {

struct V0 : Particle {

    V0() = delete;
    V0(const Particle& fit, const Seeder::PCA& pca_neg, const Seeder::PCA& pca_pos, const Vector::TrackView& neg, const Vector::TrackView& pos)
        : Particle{fit},  //
          Neg_PCAwrtV0{pca_neg},
          Pos_PCAwrtV0{pca_pos},
          Neg{neg},
          Pos{pos} {}

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
    [[nodiscard]] double ArmenterosAlpha() const {
        return Common::Math::ArmenterosAlpha(Px(), Py(), Pz(),                                         //
                                             Neg_PCAwrtV0.Px(), Neg_PCAwrtV0.Py(), Neg_PCAwrtV0.Pz(),  //
                                             Pos_PCAwrtV0.Px(), Pos_PCAwrtV0.Py(), Pos_PCAwrtV0.Pz());
    }
    [[nodiscard]] double AbsArmQtOverAlpha() const {
        double abs_alpha = std::abs(ArmenterosAlpha());
        if (abs_alpha < Common::AbsAlmostZero) return Common::BigNumber;  // protection
        return ArmenterosQt() / abs_alpha;
    }

    // Member variables //

    Seeder::PCA Neg_PCAwrtV0;
    Seeder::PCA Pos_PCAwrtV0;
    Vector::TrackView Neg;
    Vector::TrackView Pos;
};

}  // namespace R2DS::KalmanFitter
