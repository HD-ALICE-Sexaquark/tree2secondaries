#pragma once

#include <Math/Point3D.h>
#include <Math/Vector3D.h>

#include "common/Math.hpp"
#include "common/VC_V0View.hpp"

#include "KalmanFitter/KalmanFitterParticle.hxx"
#include "Seeder/BaseSeeder.hxx"

namespace R2DS::KalmanFitter {

struct ChannelA : Particle {

    ChannelA() = delete;
    ChannelA(const Particle& fit, const Seeder::PCA& pca_v0a, const Seeder::PCA& pca_v0b, const Vector::V0View& v0a, const Vector::V0View& v0b)
        : Particle{fit},  //
          V0A_PCAwrtSV{pca_v0a},
          V0B_PCAwrtSV{pca_v0b},
          V0A{v0a},
          V0B{v0b} {}

    // Member functions //

    [[nodiscard]] double SquaredDCA_V0s() const { return Common::Math::SquaredDistance(V0A_PCAwrtSV.xyz, V0B_PCAwrtSV.xyz); }
    [[nodiscard]] double SquaredDCA_V0A_SV() const { return Common::Math::SquaredDistance(V0A_PCAwrtSV.xyz, GetXYZ()); }
    [[nodiscard]] double SquaredDCA_V0B_SV() const { return Common::Math::SquaredDistance(V0B_PCAwrtSV.xyz, GetXYZ()); }

    [[nodiscard]] double DCA_V0s() const { return std::sqrt(SquaredDCA_V0s()); }
    [[nodiscard]] double DCA_V0A_SV() const { return std::sqrt(SquaredDCA_V0A_SV()); }
    [[nodiscard]] double DCA_V0B_SV() const { return std::sqrt(SquaredDCA_V0B_SV()); }

    [[nodiscard]] double CPA_V0A_SV() const {
        return Common::Math::CosinePointingAngle(V0A_PCAwrtSV.GetPxPyPz_AsROOT(), V0A_PCAwrtSV.GetXYZ_AsROOT(), GetXYZ_AsROOT());
    }
    [[nodiscard]] double CPA_V0B_SV() const {
        return Common::Math::CosinePointingAngle(V0B_PCAwrtSV.GetPxPyPz_AsROOT(), V0B_PCAwrtSV.GetXYZ_AsROOT(), GetXYZ_AsROOT());
    }

    [[nodiscard]] double DCA_Vertex(double x, double y, double z) const {
        return Common::Math::FastDCA_LineVertex(GetPxPyPz_AsROOT(), GetXYZ_AsROOT(), {x, y, z});
    }
    [[nodiscard]] double CPA_Vertex(double x, double y, double z) const {
        return Common::Math::CosinePointingAngle(GetPxPyPz_AsROOT(), GetXYZ_AsROOT(), {x, y, z});
    }

    // Member variables //

    Seeder::PCA V0A_PCAwrtSV;
    Seeder::PCA V0B_PCAwrtSV;
    Vector::V0View V0A;
    Vector::V0View V0B;
};

}  // namespace R2DS::KalmanFitter
