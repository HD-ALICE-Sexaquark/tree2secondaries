#pragma once

#include <Math/Point3D.h>
#include <Math/Vector3D.h>

#include <KFParticle.hxx>
#include <KFParticle_Math.hxx>

#include "Fit/FitTrack.hxx"
#include "Math/BaseMath.hxx"
#include "Math/Constants.hxx"
#include "View/Reconstructed/ViewV0.hxx"

namespace Tree2Secondaries::Fit {

struct alignas(T2S_SIMD_ALIGN) V0 : KF::Particle {
    // constructors //
    V0() = delete;
    V0(int entry, const Fit::Track& neg, const Fit::Track& pos)  //
        : Neg{neg}, Pos{pos}, View{nullptr, entry} {}
    V0(const View::Rec::V0& view, double mass_neg, double mass_pos)  //
        : KF::Particle{view.State(), view.CovMatrix(), 0},           //
          Neg{view.Neg, mass_neg},
          Pos{view.Pos, mass_pos},
          View{view} {}

    // main //
    void DoFit(float bz, const double* mass = nullptr) {
        AddDaughter(Neg, bz);
        AddDaughter(Pos, bz);
        if (mass != nullptr) AddMassConstraint(*mass);
    }

    // utilities //
    [[nodiscard]] KF::Vector<3> Neg_PCA_XYZ() const { return GetPCA(0).xyz; }
    [[nodiscard]] KF::Vector<3> Pos_PCA_XYZ() const { return GetPCA(1).xyz; }
    [[nodiscard]] KF::Vector<3> Neg_PCA_PxPyPz() const { return GetPCA(0).dir; }
    [[nodiscard]] KF::Vector<3> Pos_PCA_PxPyPz() const { return GetPCA(1).dir; }
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
    [[nodiscard]] double AbsEta() const { return std::abs(Eta()); }
    [[nodiscard]] double DCA_Daughters() const { return GetDCA(0, 1); }
    [[nodiscard]] double DCA_Neg_V0() const { return GetDCA(0); }
    [[nodiscard]] double DCA_Pos_V0() const { return GetDCA(1); }
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
