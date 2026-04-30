#pragma once

#include <Math/Point3D.h>
#include <Math/Vector3D.h>

#include "App/DB_Particles.hxx"
#include "KalmanFitter/KalmanFitterParticle.hxx"
#include "Math/BaseMath.hxx"
#include "Math/Constants.hxx"
#include "Seeder/BaseSeeder.hxx"
#include "View/ViewVectorV0s.hxx"

namespace Tree2Secondaries::KalmanFitter {

struct alignas(T2S_SIMD_ALIGN) ChannelA : Particle {

    ChannelA() = delete;
    ChannelA(const Particle& fit, const Seeder::PCA& pca_v0a, const Seeder::PCA& pca_v0b, const View::VecV0s& v0a, const View::VecV0s& v0b)
        : Particle{fit},  //
          V0A_PCAwrtSV{pca_v0a},
          V0B_PCAwrtSV{pca_v0b},
          V0A{v0a},
          V0B{v0b} {}

    // Member functions //

    [[nodiscard]] double E_MinusNucleon() const { return E() - Particles::Particle("Neutron").mass; }

    [[nodiscard]] double SquaredDCA_V0s() const { return Math::SquaredDistance(V0A_PCAwrtSV.xyz, V0B_PCAwrtSV.xyz); }
    [[nodiscard]] double SquaredDCA_V0A_SV() const { return Math::SquaredDistance(V0A_PCAwrtSV.xyz, GetXYZ()); }
    [[nodiscard]] double SquaredDCA_V0B_SV() const { return Math::SquaredDistance(V0B_PCAwrtSV.xyz, GetXYZ()); }

    [[nodiscard]] double DCA_V0s() const { return std::sqrt(SquaredDCA_V0s()); }
    [[nodiscard]] double DCA_V0A_SV() const { return std::sqrt(SquaredDCA_V0A_SV()); }
    [[nodiscard]] double DCA_V0B_SV() const { return std::sqrt(SquaredDCA_V0B_SV()); }

    [[nodiscard]] double CPA_V0A_SV() const {
        return Math::CosinePointingAngle(V0A_PCAwrtSV.GetPxPyPz_AsROOT(), V0A_PCAwrtSV.GetXYZ_AsROOT(), GetXYZ_AsROOT());
    }
    [[nodiscard]] double CPA_V0B_SV() const {
        return Math::CosinePointingAngle(V0B_PCAwrtSV.GetPxPyPz_AsROOT(), V0B_PCAwrtSV.GetXYZ_AsROOT(), GetXYZ_AsROOT());
    }

    [[nodiscard]] double DCA_Vertex(double x, double y, double z) const {
        return Math::FastDCALineVertex(GetPxPyPz_AsROOT(), GetXYZ_AsROOT(), {x, y, z});
    }
    [[nodiscard]] double CPA_Vertex(double x, double y, double z) const {
        return Math::CosinePointingAngle(GetPxPyPz_AsROOT(), GetXYZ_AsROOT(), {x, y, z});
    }

    // Member variables //

    Seeder::PCA V0A_PCAwrtSV;
    Seeder::PCA V0B_PCAwrtSV;
    View::VecV0s V0A;
    View::VecV0s V0B;
};

}  // namespace Tree2Secondaries::KalmanFitter
