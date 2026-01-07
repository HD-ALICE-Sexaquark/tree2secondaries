#pragma once

#include <KFParticle.hxx>
#include <KFParticle_Math.hxx>

#include "Fit/FitSexaquark.hxx"
#include "Fit/FitV0.hxx"
#include "Math/BaseMath.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::Fit {

struct alignas(T2S_SIMD_ALIGN) ChannelA : Fit::Sexaquark {
    // constructors //
    ChannelA() = delete;
    ChannelA(const Fit::V0& v0a, const Fit::V0& v0b)  //
        : Sexaquark{Const::Particle_Mass[EParticle::Neutron]}, V0A{v0a}, V0B{v0b} {}

    // main //
    void DoFit(float bz) {
        AddDaughter(V0A, bz);
        AddDaughter(V0B, bz);
    }

    // utilities //
    [[nodiscard]] KF::Vector<3> V0A_PCA_XYZ() const { return GetPCA(0).xyz; }
    [[nodiscard]] KF::Vector<3> V0A_PCA_PxPyPz() const { return GetPCA(0).dir; }
    [[nodiscard]] KF::Vector<3> V0B_PCA_XYZ() const { return GetPCA(1).xyz; }
    [[nodiscard]] KF::Vector<3> V0B_PCA_PxPyPz() const { return GetPCA(1).dir; }

    // cuts //
    [[nodiscard]] double DecayLength_V0A() const {
        return KF::Math::Norm(KF::Vector<3>{V0A.X() - V0A_PCA_XYZ()[0], V0A.Y() - V0A_PCA_XYZ()[1], V0A.Z() - V0A_PCA_XYZ()[2]});
    }
    [[nodiscard]] double DecayLength_V0B() const {
        return KF::Math::Norm(KF::Vector<3>{V0B.X() - V0B_PCA_XYZ()[0], V0B.Y() - V0B_PCA_XYZ()[1], V0B.Z() - V0B_PCA_XYZ()[2]});
    }
    [[nodiscard]] double DCA_btw_V0s() const { return GetDCA(0, 1); };
    [[nodiscard]] double DCA_V0A_wrt_SV() const { return GetDCA(0); };
    [[nodiscard]] double DCA_V0B_wrt_SV() const { return GetDCA(1); };
    [[nodiscard]] double CPA_V0A_wrt_SV() const {
        return Tree2Secondaries::Math::CosinePointingAngle(V0A.GetPxPyPz_AsROOT(), V0A.GetXYZ_AsROOT(), GetXYZ_AsROOT());
    }
    [[nodiscard]] double CPA_V0B_wrt_SV() const {
        return Tree2Secondaries::Math::CosinePointingAngle(V0A.GetPxPyPz_AsROOT(), V0A.GetXYZ_AsROOT(), GetXYZ_AsROOT());
    }

    // daughters //
    Fit::V0 V0A;
    Fit::V0 V0B;
};

}  // namespace Tree2Secondaries::Fit
