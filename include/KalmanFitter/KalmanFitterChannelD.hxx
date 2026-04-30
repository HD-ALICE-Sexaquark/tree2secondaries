#pragma once

#include <Math/Point3D.h>
#include <Math/Vector3D.h>

#include "App/DB_Particles.hxx"
#include "KalmanFitter/KalmanFitterParticle.hxx"
#include "Math/BaseMath.hxx"
#include "Math/Constants.hxx"
#include "Seeder/BaseSeeder.hxx"
#include "View/ViewVectorTracks.hxx"
#include "View/ViewVectorV0s.hxx"

namespace Tree2Secondaries::KalmanFitter {

struct alignas(T2S_SIMD_ALIGN) ChannelD : Particle {

    ChannelD() = delete;
    ChannelD(const Particle& fit, const Seeder::PCA& pca_ka, const Seeder::PCA& pca_v0, const View::VecTracks& ka, const View::VecV0s& v0)
        : Particle{fit},  //
          Kaon_PCAwrtSV{pca_ka},
          V0_PCAwrtSV{pca_v0},
          V0{v0},
          Kaon{ka} {}

    // Member functions //

    [[nodiscard]] double E_MinusNucleon() const { return E() - Particles::Particle("Proton").mass; }

    [[nodiscard]] double SquaredDCA_Kaon_V0() const { return Math::SquaredDistance(Kaon_PCAwrtSV.xyz, V0_PCAwrtSV.xyz); }
    [[nodiscard]] double SquaredDCA_V0_SV() const { return Math::SquaredDistance(V0_PCAwrtSV.xyz, GetXYZ()); }
    [[nodiscard]] double SquaredDCA_Kaon_SV() const { return Math::SquaredDistance(Kaon_PCAwrtSV.xyz, GetXYZ()); }

    [[nodiscard]] double DCA_Kaon_V0() const { return std::sqrt(SquaredDCA_Kaon_V0()); }
    [[nodiscard]] double DCA_V0_SV() const { return std::sqrt(SquaredDCA_V0_SV()); }
    [[nodiscard]] double DCA_Kaon_SV() const { return std::sqrt(SquaredDCA_Kaon_SV()); }

    [[nodiscard]] double CPA_V0_SV() const {
        return Math::CosinePointingAngle(V0_PCAwrtSV.GetPxPyPz_AsROOT(), V0_PCAwrtSV.GetXYZ_AsROOT(), GetXYZ_AsROOT());
    }

    [[nodiscard]] double DCA_Vertex(double x, double y, double z) const {
        return Math::FastDCALineVertex(GetPxPyPz_AsROOT(), GetXYZ_AsROOT(), {x, y, z});
    }
    [[nodiscard]] double CPA_Vertex(double x, double y, double z) const {
        return Math::CosinePointingAngle(GetPxPyPz_AsROOT(), GetXYZ_AsROOT(), {x, y, z});
    }

    // Member variables //

    Seeder::PCA Kaon_PCAwrtSV;
    Seeder::PCA V0_PCAwrtSV;
    View::VecV0s V0;
    View::VecTracks Kaon;
};

}  // namespace Tree2Secondaries::KalmanFitter
