#pragma once

#include <Math/Point3D.h>
#include <Math/Vector3D.h>

#include "KalmanFitter/KalmanFitterParticle.hxx"
#include "Math/BaseMath.hxx"
#include "Math/Constants.hxx"
#include "Seeder/BaseSeeder.hxx"
#include "View/Reconstructed/ViewRecTrack.hxx"
#include "View/Reconstructed/ViewRecV0.hxx"

namespace Tree2Secondaries::KalmanFitter {

struct alignas(T2S_SIMD_ALIGN) ChannelD : KalmanFitter::Particle {

    ChannelD() = delete;
    ChannelD(const Particle& fit, const Seeder::PCA& pca_ka, const Seeder::PCA& pca_v0, const View::Rec::Track& ka, const View::Rec::V0& v0)
        : Particle{fit}, Kaon_at_PCA{pca_ka}, V0_at_PCA{pca_v0}, V0{v0}, Kaon{ka} {}

    // Member functions //

    [[nodiscard]] double E_MinusNucleon() const { return E() - Const::Particle_Mass[PID_StableParticle::Proton]; }

    [[nodiscard]] double SquaredDCA_Kaon_V0() const { return Math::SquaredDistance(Kaon_at_PCA.xyz, V0_at_PCA.xyz); }
    [[nodiscard]] double SquaredDCA_V0_SV() const { return Math::SquaredDistance(V0_at_PCA.xyz, GetXYZ()); }
    [[nodiscard]] double SquaredDCA_Kaon_SV() const { return Math::SquaredDistance(Kaon_at_PCA.xyz, GetXYZ()); }

    [[nodiscard]] double DCA_Kaon_V0() const { return std::sqrt(SquaredDCA_Kaon_V0()); }
    [[nodiscard]] double DCA_V0_SV() const { return std::sqrt(SquaredDCA_V0_SV()); }
    [[nodiscard]] double DCA_Kaon_SV() const { return std::sqrt(SquaredDCA_Kaon_SV()); }

    [[nodiscard]] double CPA_V0_SV() const {
        return Math::CosinePointingAngle(V0_at_PCA.GetPxPyPz_AsROOT(), V0_at_PCA.GetXYZ_AsROOT(), GetXYZ_AsROOT());
    };

    [[nodiscard]] double DCA_Vertex(double x, double y, double z) const {
        return Math::FastDCALineVertex(GetPxPyPz_AsROOT(), GetXYZ_AsROOT(), {x, y, z});
    }
    [[nodiscard]] double CPA_Vertex(double x, double y, double z) const {
        return Math::CosinePointingAngle(GetPxPyPz_AsROOT(), GetXYZ_AsROOT(), {x, y, z});
    }

    // Member variables //

    Seeder::PCA Kaon_at_PCA;
    Seeder::PCA V0_at_PCA;
    View::Rec::V0 V0;
    View::Rec::Track Kaon;
};

}  // namespace Tree2Secondaries::KalmanFitter
