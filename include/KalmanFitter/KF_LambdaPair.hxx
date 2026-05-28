#pragma once

#include <Math/Point3D.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>

#include "common/Math.hpp"
#include "common/VC_OnTheFlyLambdaView.hpp"

#include "Seeder/BaseSeeder.hxx"

namespace R2DS::KalmanFitter {

// NOTE: not really a `KF::Particle`, but this could change in the future
// For now, it's just a simple wrapper
struct LambdaPair : ROOT::Math::PxPyPzEVector {

    // Constructors //

    LambdaPair() = delete;
    LambdaPair(const Vector::OnTheFlyLambdaView& lambda1, const Vector::OnTheFlyLambdaView& lambda2, const Seeder::PCA& pca_lambda1,
               const Seeder::PCA& pca_lambda2)
        : ROOT::Math::PxPyPzEVector{lambda1.lv + lambda2.lv},
          Lambda1{lambda1},  //
          Lambda2{lambda2},
          Lambda1_PCAwrtDV{pca_lambda1},
          Lambda2_PCAwrtDV{pca_lambda2} {}

    // Member functions //

    [[nodiscard]] double SquaredDCA_Lambdas() const { return Common::Math::SquaredDistance(Lambda1_PCAwrtDV.xyz, Lambda2_PCAwrtDV.xyz); }
    [[nodiscard]] double SquaredDCA_Lambda1_DV() const { return Common::Math::SquaredDistance(Lambda1_PCAwrtDV.xyz, {DV.X(), DV.Y(), DV.Z()}); }
    [[nodiscard]] double SquaredDCA_Lambda2_DV() const { return Common::Math::SquaredDistance(Lambda2_PCAwrtDV.xyz, {DV.X(), DV.Y(), DV.Z()}); }

    [[nodiscard]] double DCA_Lambdas() const { return std::sqrt(SquaredDCA_Lambdas()); }
    [[nodiscard]] double DCA_Lambda1_DV() const { return std::sqrt(SquaredDCA_Lambda1_DV()); }
    [[nodiscard]] double DCA_Lambda2_DV() const { return std::sqrt(SquaredDCA_Lambda2_DV()); }

    [[nodiscard]] double CPA_Lambda1_DV() const {
        return Common::Math::CosinePointingAngle(Lambda1_PCAwrtDV.GetPxPyPz_AsROOT(), Lambda1_PCAwrtDV.GetXYZ_AsROOT(), DV);
    }
    [[nodiscard]] double CPA_Lambda2_DV() const {
        return Common::Math::CosinePointingAngle(Lambda2_PCAwrtDV.GetPxPyPz_AsROOT(), Lambda2_PCAwrtDV.GetXYZ_AsROOT(), DV);
    }

    [[nodiscard]] double DCA_Vertex(double x, double y, double z) const {  //
        return Common::Math::FastDCA_LineVertex(Vect(), DV, {x, y, z});
    }
    [[nodiscard]] double CPA_Vertex(double x, double y, double z) const {  //
        return Common::Math::CosinePointingAngle(Vect(), DV, {x, y, z});
    }

    // Member variables //

    Vector::OnTheFlyLambdaView Lambda1;
    Vector::OnTheFlyLambdaView Lambda2;
    Seeder::PCA Lambda1_PCAwrtDV;
    Seeder::PCA Lambda2_PCAwrtDV;
    ROOT::Math::XYZPoint DV;  // decay vertex
};

}  // namespace R2DS::KalmanFitter
