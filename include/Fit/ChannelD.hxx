#pragma once

#include <KFParticle.hxx>
#include <KFParticle_Math.hxx>

#include "Fit/Sexaquark.hxx"
#include "Fit/Track.hxx"
#include "Fit/V0.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::Fit {

struct alignas(T2S_SIMD_ALIGN) ChannelD : Fit::Sexaquark {
    // constructors //
    ChannelD() = delete;
    ChannelD(const Fit::V0& v0, const Fit::Track& kaon)        //
        : Sexaquark{Const::Particle_Mass[EParticle::Proton]},  //
          V0{v0},
          Kaon{kaon} {};

    // main //
    void DoFit(float bz) {
        AddDaughter(V0, bz);
        AddDaughter(Kaon, bz);
    }

    // utilities //
    [[nodiscard]] KF::Vector<3> V0_PCA_XYZ() const { return GetPCA(0).xyz; };
    [[nodiscard]] KF::Vector<3> V0_PCA_PxPyPz() const { return GetPCA(0).dir; };
    [[nodiscard]] KF::Vector<3> Kaon_PCA_XYZ() const { return GetPCA(1).xyz; };
    [[nodiscard]] KF::Vector<3> Kaon_PCA_PxPyPz() const { return GetPCA(1).dir; };
    [[nodiscard]] double Kaon_PCA_E() const {
        KF::Vector<3> pca_p{Kaon_PCA_PxPyPz()};
        return std::sqrt(Kaon.Mass() * Kaon.Mass() + pca_p[0] * pca_p[0] + pca_p[1] * pca_p[1] + pca_p[2] * pca_p[2]);
    };

    // cuts //
    [[nodiscard]] double DCA_btw_V0_Kaon() const { return GetDCA(0, 1); };
    [[nodiscard]] double DCA_V0_wrt_SV() const { return GetDCA(0); };
    [[nodiscard]] double DCA_Kaon_wrt_SV() const { return GetDCA(1); };
    [[nodiscard]] double DCA_V0Neg_wrt_SV(float bz) const { return Tree2Secondaries::Math::FastDCAHelixVertex(V0.Neg, X(), Y(), Z(), bz); };
    [[nodiscard]] double DCA_V0Pos_wrt_SV(float bz) const { return Tree2Secondaries::Math::FastDCAHelixVertex(V0.Pos, X(), Y(), Z(), bz); };
    [[nodiscard]] double DCA_btw_V0Neg_Kaon(float bz) const { return Tree2Secondaries::Math::FastDCAHelixHelix(V0.Neg, Kaon, bz); };
    [[nodiscard]] double DCA_btw_V0Pos_Kaon(float bz) const { return Tree2Secondaries::Math::FastDCAHelixHelix(V0.Pos, Kaon, bz); };

    // daughters //
    Fit::V0 V0;
    Fit::Track Kaon;
};

}  // namespace Tree2Secondaries::Fit
