#pragma once

#include <cmath>

#include <Math/Point3D.h>
#include <Math/Vector3D.h>

#include "KalmanFitter/KalmanFitterParticle.hxx"
#include "Math/BaseMath.hxx"
#include "Math/Constants.hxx"
#include "Seeder/BaseSeeder.hxx"
#include "View/ViewVectorTracks.hxx"

namespace Tree2Secondaries::KalmanFitter {

struct alignas(T2S_SIMD_ALIGN) V0 : Particle {

    V0() = delete;
    V0(const Particle& fit, const Seeder::PCA& pca_neg, const Seeder::PCA& pca_pos, const View::VecTracks& track_neg,
       const View::VecTracks& track_pos)
        : Particle{fit},  //
          Neg_PCAwrtV0{pca_neg},
          Pos_PCAwrtV0{pca_pos},
          Neg{track_neg},
          Pos{track_pos} {}

    // Member functions //

    [[nodiscard]] double SquaredDCA_Daughters() const { return Math::SquaredDistance(Neg_PCAwrtV0.xyz, Pos_PCAwrtV0.xyz); }
    [[nodiscard]] double SquaredDCA_Neg_V0() const { return Math::SquaredDistance(Neg_PCAwrtV0.xyz, GetXYZ()); }
    [[nodiscard]] double SquaredDCA_Pos_V0() const { return Math::SquaredDistance(Pos_PCAwrtV0.xyz, GetXYZ()); }

    [[nodiscard]] double DCA_Daughters() const { return std::sqrt(SquaredDCA_Daughters()); }
    [[nodiscard]] double DCA_Neg_V0() const { return std::sqrt(SquaredDCA_Neg_V0()); }
    [[nodiscard]] double DCA_Pos_V0() const { return std::sqrt(SquaredDCA_Pos_V0()); }

    [[nodiscard]] double DCA_Vertex(double x, double y, double z) const {
        return Math::FastDCALineVertex(GetPxPyPz_AsROOT(), GetXYZ_AsROOT(), {x, y, z});
    }
    [[nodiscard]] double CPA_Vertex(double x, double y, double z) const {
        return Math::CosinePointingAngle(GetPxPyPz_AsROOT(), GetXYZ_AsROOT(), {x, y, z});
    }

    [[nodiscard]] double ArmenterosQt() const {
        return Math::ArmenterosQt(Px(), Py(), Pz(),  //
                                  Neg_PCAwrtV0.Px(), Neg_PCAwrtV0.Py(), Neg_PCAwrtV0.Pz());
    }
    [[nodiscard]] double ArmenterosAlpha() const {
        return Math::ArmenterosAlpha(Px(), Py(), Pz(),                                         //
                                     Neg_PCAwrtV0.Px(), Neg_PCAwrtV0.Py(), Neg_PCAwrtV0.Pz(),  //
                                     Pos_PCAwrtV0.Px(), Pos_PCAwrtV0.Py(), Pos_PCAwrtV0.Pz());
    }
    [[nodiscard]] double AbsArmQtOverAlpha() const {
        double abs_alpha = std::abs(ArmenterosAlpha());
        if (abs_alpha < Const::AbsAlmostZero) return Const::BigNumber;  // protection
        return ArmenterosQt() / abs_alpha;
    }

    // Member variables //

    Seeder::PCA Neg_PCAwrtV0;
    Seeder::PCA Pos_PCAwrtV0;
    View::VecTracks Neg;
    View::VecTracks Pos;
};

}  // namespace Tree2Secondaries::KalmanFitter
