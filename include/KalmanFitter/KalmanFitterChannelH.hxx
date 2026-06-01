#pragma once

#include <Math/Point3Dfwd.h>
#include <Math/Vector3Dfwd.h>

#include "common/Math.hpp"
#include "common/POD_Track.hpp"

#include "KalmanFitter/KalmanFitterParticle.hxx"
#include "Seeder/BaseSeeder.hxx"

namespace R2DS::KF {

struct ChannelH : Particle {

    ChannelH() = delete;
    ChannelH(const Particle& fit, const Seeder::PCA& pca_ka1, const Seeder::PCA& pca_ka2, const POD::Track& ka1, const POD::Track& ka2)
        : Particle{fit},  //
          Kaon1_PCAwrtSV{pca_ka1},
          Kaon2_PCAwrtSV{pca_ka2},
          Kaon1{&ka1},
          Kaon2{&ka2} {}

    // Member functions //

    [[nodiscard]] double SquaredDCA_Kaons() const { return Common::Math::SquaredDistance(Kaon1_PCAwrtSV.xyz, Kaon2_PCAwrtSV.xyz); }
    [[nodiscard]] double SquaredDCA_Kaon1_SV() const { return Common::Math::SquaredDistance(Kaon1_PCAwrtSV.xyz, GetXYZ()); }
    [[nodiscard]] double SquaredDCA_Kaon2_SV() const { return Common::Math::SquaredDistance(Kaon2_PCAwrtSV.xyz, GetXYZ()); }

    [[nodiscard]] double DCA_Kaons() const { return std::sqrt(SquaredDCA_Kaons()); }
    [[nodiscard]] double DCA_Kaon1_SV() const { return std::sqrt(SquaredDCA_Kaon1_SV()); }
    [[nodiscard]] double DCA_Kaon2_SV() const { return std::sqrt(SquaredDCA_Kaon2_SV()); }

    [[nodiscard]] double DCA_Vertex(double x, double y, double z) const {
        return Common::Math::FastDCA_LineVertex(GetPxPyPz_AsROOT(), GetXYZ_AsROOT(), {x, y, z});
    }
    [[nodiscard]] double CPA_Vertex(double x, double y, double z) const {
        return Common::Math::CosinePointingAngle(GetPxPyPz_AsROOT(), GetXYZ_AsROOT(), {x, y, z});
    }

    // Member variables //

    Seeder::PCA Kaon1_PCAwrtSV;
    Seeder::PCA Kaon2_PCAwrtSV;
    const POD::Track* Kaon1;
    const POD::Track* Kaon2;
};

}  // namespace R2DS::KF
