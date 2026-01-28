#pragma once

#include <array>
#include <utility>

#include <KFParticle.h>
#include <KFParticleBase.h>

#include "Fit/Legacy/FitSexaquark_Legacy.hxx"
#include "Fit/Legacy/FitV0_Legacy.hxx"
#include "Math/Constants.hxx"
#include "Math/Legacy/BaseMath_Legacy.hxx"

namespace Tree2Secondaries::Fit {

struct alignas(T2S_SIMD_ALIGN) ChannelA : Fit::Sexaquark {
    // constructors //
    ChannelA() = delete;
    ChannelA(Fit::V0 v0a, Fit::V0 v0b)  //
        : Sexaquark{Const::Particle_Mass[EParticle::Neutron]}, V0A{std::move(v0a)}, V0B{std::move(v0b)} {}

    // main //
    void DoFit(float bz) {
        KFParticle::SetField(bz);
        AddDaughter(V0A);
        AddDaughter(V0B);
    }

    // utilities //
    [[nodiscard]] std::array<float, 3> V0A_PCA_XYZ() const { return {0., 0., 0.}; }     // PENDING
    [[nodiscard]] std::array<float, 3> V0A_PCA_PxPyPz() const { return {0., 0., 0.}; }  // PENDING
    [[nodiscard]] std::array<float, 3> V0B_PCA_XYZ() const { return {0., 0., 0.}; }     // PENDING
    [[nodiscard]] std::array<float, 3> V0B_PCA_PxPyPz() const { return {0., 0., 0.}; }  // PENDING

    // cuts //
    [[nodiscard]] double DecayLength_V0A() const { return 0.; }  // PENDING
    [[nodiscard]] double DecayLength_V0B() const { return 0.; }  // PENDING
    [[nodiscard]] double DCA_btw_V0s() const { return V0A.GetDistanceFromParticle(V0B); }
    [[nodiscard]] double DCA_V0A_wrt_SV() const { return V0A.GetDistanceFromVertex(*this); }
    [[nodiscard]] double DCA_V0B_wrt_SV() const { return V0B.GetDistanceFromVertex(*this); }
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
