#pragma once

#include <KFParticle.hxx>
#include <KFParticle_Math.hxx>

#include "Fit/Track.hxx"
#include "Math/Constants.hxx"
#include "Math/Math.hxx"

namespace Tree2Secondaries::Fit {

struct alignas(T2S_SIMD_ALIGN) V0 : KF::Particle {
    // constructors //
    V0() = delete;
    V0(int entry, const Fit::Track& neg, const Fit::Track& pos) : Neg(neg), Pos(pos), Entry{entry} {}
    V0(int entry, const KF::Vector<7>& p, const KF::SymMatrix<7>& cov, const Fit::Track& neg, const Fit::Track& pos)
        : KF::Particle{p, cov, 0}, Neg(neg), Pos(pos), Entry{entry} {}

    // main //
    void DoFit(float bz, const double* mass = nullptr) {
        AddDaughter(Neg, bz);
        AddDaughter(Pos, bz);
        if (mass != nullptr) AddMassConstraint(*mass);
    }

    // utilities //
    [[nodiscard]] KF::Vector<3> Neg_PCA_XYZ() const { return GetPCA(0).xyz; };
    [[nodiscard]] KF::Vector<3> Pos_PCA_XYZ() const { return GetPCA(1).xyz; };
    [[nodiscard]] KF::Vector<3> Neg_PCA_PxPyPz() const { return GetPCA(0).dir; };
    [[nodiscard]] KF::Vector<3> Pos_PCA_PxPyPz() const { return GetPCA(1).dir; };

    [[nodiscard]] double ArmenterosQt() const {
        return Math::ArmenterosQt(Px(), Py(), Pz(),  //
                                  Neg_PCA_PxPyPz()[0], Neg_PCA_PxPyPz()[1], Neg_PCA_PxPyPz()[2]);
    }
    [[nodiscard]] double ArmenterosAlpha() const {
        return Math::ArmenterosAlpha(Px(), Py(), Pz(),                                               //
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
        return Math::FastDCALineVertex({Px(), Py(), Pz()}, {X(), Y(), Z()}, {x, y, z});
    }
    [[nodiscard]] double CPA_Point(float x, float y, float z) const {
        return Math::CosinePointingAngle({Px(), Py(), Pz()}, {X(), Y(), Z()}, {x, y, z});
    }

    // member vars //
    Fit::Track Neg;
    Fit::Track Pos;
    int Entry{};
};

}  // namespace Tree2Secondaries::Fit
