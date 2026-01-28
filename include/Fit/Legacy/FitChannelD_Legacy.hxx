#pragma once

#include <array>
#include <utility>

#include <KFParticle.h>

#include "Fit/Legacy/FitSexaquark_Legacy.hxx"
#include "Fit/Legacy/FitTrack_Legacy.hxx"
#include "Fit/Legacy/FitV0_Legacy.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::Fit {

struct alignas(T2S_SIMD_ALIGN) ChannelD : Fit::Sexaquark {
    // constructors //
    ChannelD() = delete;
    ChannelD(Fit::V0 v0, Fit::Track kaon)                      //
        : Sexaquark{Const::Particle_Mass[EParticle::Proton]},  //
          V0{std::move(v0)},
          Kaon{std::move(kaon)} {}

    // main //
    void DoFit(float bz) {
        KFParticle::SetField(bz);
        AddDaughter(V0);
        AddDaughter(Kaon);
    }

    // utilities //
    [[nodiscard]] std::array<float, 3> V0_PCA_XYZ() const { return {0., 0., 0.}; }       // PENDING
    [[nodiscard]] std::array<float, 3> V0_PCA_PxPyPz() const { return {0., 0., 0.}; }    // PENDING
    [[nodiscard]] std::array<float, 3> Kaon_PCA_XYZ() const { return {0., 0., 0.}; }     // PENDING
    [[nodiscard]] std::array<float, 3> Kaon_PCA_PxPyPz() const { return {0., 0., 0.}; }  // PENDING
    [[nodiscard]] double Kaon_PCA_E() const {
        std::array<float, 3> pca_p{Kaon_PCA_PxPyPz()};
        return std::sqrt(Kaon.GetMass() * Kaon.GetMass() + pca_p[0] * pca_p[0] + pca_p[1] * pca_p[1] + pca_p[2] * pca_p[2]);
    }

    // cuts //
    [[nodiscard]] double DCA_btw_V0_Kaon() const { return V0.GetDistanceFromParticle(Kaon); }
    [[nodiscard]] double DCA_V0_wrt_SV() const { return V0.GetDistanceFromVertex(*this); }
    [[nodiscard]] double DCA_Kaon_wrt_SV() const { return Kaon.GetDistanceFromVertex(*this); }
    [[nodiscard]] double DCA_V0Neg_wrt_SV(float bz) const { return Tree2Secondaries::Math::FastDCAHelixVertex(V0.Neg, X(), Y(), Z(), bz); }
    [[nodiscard]] double DCA_V0Pos_wrt_SV(float bz) const { return Tree2Secondaries::Math::FastDCAHelixVertex(V0.Pos, X(), Y(), Z(), bz); }
    [[nodiscard]] double DCA_btw_V0Neg_Kaon(float bz) const { return Tree2Secondaries::Math::FastDCAHelixHelix(V0.Neg, Kaon, bz); }
    [[nodiscard]] double DCA_btw_V0Pos_Kaon(float bz) const { return Tree2Secondaries::Math::FastDCAHelixHelix(V0.Pos, Kaon, bz); }

    // daughters //
    Fit::V0 V0;
    Fit::Track Kaon;
};

}  // namespace Tree2Secondaries::Fit
