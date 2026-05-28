#pragma once

#include <Math/Point3D.h>
#include <Math/Vector3D.h>

#include "common/Math.hpp"
#include "common/VC_TrackView.hpp"
#include "common/VC_V0View.hpp"

#include "KalmanFitter/KalmanFitterParticle.hxx"
#include "Seeder/BaseSeeder.hxx"

namespace R2DS::KalmanFitter {

struct ChannelD : Particle {

    ChannelD() = delete;
    ChannelD(const Particle& fit, const Seeder::PCA& pca_v0, const Seeder::PCA& pca_ka, const Vector::V0View& v0, const Vector::TrackView& ka)
        : Particle{fit},  //
          V0_PCAwrtSV{pca_v0},
          Kaon_PCAwrtSV{pca_ka},
          V0{v0},
          Kaon{ka} {}

    // Member functions //

    [[nodiscard]] double SquaredDCA_V0_SV() const { return Common::Math::SquaredDistance(V0_PCAwrtSV.xyz, GetXYZ()); }
    [[nodiscard]] double SquaredDCA_Kaon_SV() const { return Common::Math::SquaredDistance(Kaon_PCAwrtSV.xyz, GetXYZ()); }
    [[nodiscard]] double SquaredDCA_Kaon_V0() const { return Common::Math::SquaredDistance(Kaon_PCAwrtSV.xyz, V0_PCAwrtSV.xyz); }

    [[nodiscard]] double DCA_V0_SV() const { return std::sqrt(SquaredDCA_V0_SV()); }
    [[nodiscard]] double DCA_Kaon_SV() const { return std::sqrt(SquaredDCA_Kaon_SV()); }
    [[nodiscard]] double DCA_Kaon_V0() const { return std::sqrt(SquaredDCA_Kaon_V0()); }

    [[nodiscard]] double CPA_V0_SV() const {
        return Common::Math::CosinePointingAngle(V0_PCAwrtSV.GetPxPyPz_AsROOT(), V0_PCAwrtSV.GetXYZ_AsROOT(), GetXYZ_AsROOT());
    }

    [[nodiscard]] double DCA_Vertex(double x, double y, double z) const {
        return Common::Math::FastDCA_LineVertex(GetPxPyPz_AsROOT(), GetXYZ_AsROOT(), {x, y, z});
    }
    [[nodiscard]] double CPA_Vertex(double x, double y, double z) const {
        return Common::Math::CosinePointingAngle(GetPxPyPz_AsROOT(), GetXYZ_AsROOT(), {x, y, z});
    }

    // Member variables //

    Seeder::PCA V0_PCAwrtSV;
    Seeder::PCA Kaon_PCAwrtSV;
    Vector::V0View V0;
    Vector::TrackView Kaon;
};

}  // namespace R2DS::KalmanFitter
