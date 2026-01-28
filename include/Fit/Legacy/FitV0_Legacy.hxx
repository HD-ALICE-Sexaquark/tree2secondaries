#pragma once

#include <array>
#include <cmath>
#include <utility>

#include <Math/Point3D.h>
#include <Math/Vector3D.h>

#include <KFPVertex.h>
#include <KFParticle.h>

#include "App/Logger.hxx"
#include "Fit/Legacy/FitTrack_Legacy.hxx"
#include "Math/Constants.hxx"
#include "Math/Legacy/BaseMath_Legacy.hxx"
#include "View/Reconstructed/Legacy/ViewV0_Legacy.hxx"

namespace Tree2Secondaries::Fit {

struct alignas(T2S_SIMD_ALIGN) V0 : KFParticle {
    // constructors //
    V0() = delete;
    V0(int entry, Fit::Track neg, Fit::Track pos)  //
        : Neg{std::move(neg)}, Pos{std::move(pos)}, View{nullptr, entry} {}
    V0(const View::Rec::V0& view, double mass_neg, double mass_pos)  //
        : Neg{view.Neg, mass_neg},                                   //
          Pos{view.Pos, mass_pos},
          View{view} {

        std::array<float, 7> state{view.State()};
        std::array<float, 28> cov_matrix{view.CovMatrix()};

        double mass_v0{std::sqrt(view.Energy() * view.Energy() - view.Px() * view.Px() - view.Py() * view.Py() - view.Pz() * view.Px())};
        Create(state.data(), cov_matrix.data(), 0, static_cast<float>(mass_v0));

        Parameter(6) = state[6];

        Covariance(21) = cov_matrix[21];
        Covariance(22) = cov_matrix[22];
        Covariance(23) = cov_matrix[23];
        Covariance(24) = cov_matrix[24];
        Covariance(25) = cov_matrix[25];
        Covariance(26) = cov_matrix[26];
        Covariance(27) = cov_matrix[27];
    }

    // main //
    void DoFit(float bz, const double* mass = nullptr) {
        KFParticle::SetField(bz);
        AddDaughter(Neg);
        AddDaughter(Pos);
        if (mass != nullptr) SetMassConstraint(static_cast<float>(*mass));

        // Logger::Info(__FUNCTION__, "before transport, (px,py,pz)(neg)=({},{},{})", Neg.Px(), Neg.Py(), Neg.Pz());
        Neg.TransportToParticle(Pos);
        // Logger::Info(__FUNCTION__, "after transport, (px,py,pz)(neg)=({},{},{})", Neg.Px(), Neg.Py(), Neg.Pz());

        // Logger::Info(__FUNCTION__, "before transport, (px,py,pz)(pos)=({},{},{})", Pos.Px(), Pos.Py(), Pos.Pz());
        Pos.TransportToParticle(Neg);
        // Logger::Info(__FUNCTION__, "after transport, (px,py,pz)(pos)=({},{},{})", Pos.Px(), Pos.Py(), Pos.Pz());
    }

    // utilities //
    [[nodiscard]] std::array<float, 3> Neg_PCA_XYZ() const { return {Neg.X(), Neg.Y(), Neg.Z()}; }
    [[nodiscard]] std::array<float, 3> Pos_PCA_XYZ() const { return {Pos.X(), Pos.Y(), Pos.Z()}; }
    [[nodiscard]] std::array<float, 3> Neg_PCA_PxPyPz() const { return {Neg.Px(), Neg.Py(), Neg.Pz()}; }
    [[nodiscard]] std::array<float, 3> Pos_PCA_PxPyPz() const { return {Neg.Px(), Neg.Py(), Neg.Pz()}; }
    [[nodiscard]] ROOT::Math::XYZPoint GetXYZ_AsROOT() const { return {X(), Y(), Z()}; }
    [[nodiscard]] ROOT::Math::XYZVector GetPxPyPz_AsROOT() const { return {Px(), Py(), Pz()}; }

    [[nodiscard]] double ArmenterosQt() const {
        return Tree2Secondaries::Math::ArmenterosQt(Px(), Py(), Pz(),  //
                                                    Neg_PCA_PxPyPz()[0], Neg_PCA_PxPyPz()[1], Neg_PCA_PxPyPz()[2]);
    }
    [[nodiscard]] double ArmenterosAlpha() const {
        return Tree2Secondaries::Math::ArmenterosAlpha(Px(), Py(), Pz(),                                               //
                                                       Neg_PCA_PxPyPz()[0], Neg_PCA_PxPyPz()[1], Neg_PCA_PxPyPz()[2],  //
                                                       Pos_PCA_PxPyPz()[0], Pos_PCA_PxPyPz()[1], Pos_PCA_PxPyPz()[2]);
    }

    // cuts //
    [[nodiscard]] double AbsZ() const { return std::abs(Z()); }
    [[nodiscard]] double AbsEta() const { return std::abs(GetEta()); }
    [[nodiscard]] double DCA_Daughters() const { return Neg.GetDistanceFromParticle(Pos); }
    [[nodiscard]] double DCA_Neg_V0() const { return GetDistanceFromParticle(Neg); }
    [[nodiscard]] double DCA_Pos_V0() const { return GetDistanceFromParticle(Pos); }
    [[nodiscard]] double AbsArmQtOverAlpha() const { return ArmenterosQt() / std::abs(ArmenterosAlpha()); };
    [[nodiscard]] double DCA_Point(float x, float y, float z) const {
        return Tree2Secondaries::Math::FastDCALineVertex({Px(), Py(), Pz()}, {X(), Y(), Z()}, {x, y, z});
    }
    [[nodiscard]] double CPA_Point(float x, float y, float z) const {
        return Tree2Secondaries::Math::CosinePointingAngle({Px(), Py(), Pz()}, {X(), Y(), Z()}, {x, y, z});
    }

    // member vars //
    Fit::Track Neg;
    Fit::Track Pos;
    View::Rec::V0 View;
};

}  // namespace Tree2Secondaries::Fit
