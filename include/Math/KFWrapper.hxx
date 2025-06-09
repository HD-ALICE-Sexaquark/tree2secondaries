#ifndef T2S_KF_WRAPPER_HXX
#define T2S_KF_WRAPPER_HXX

#include <array>
#include <cmath>
#include <cstddef>

#include "KFParticle.hxx"
#include "KFParticle_Math.hxx"

#include "Math/Constants.hxx"
#include "Math/Math.hxx"
#include "Structures/Events.hxx"
#include "Structures/PackedEvents.hxx"

namespace KF {

// Extensions of `KF::Particle` //

struct alignas(32) Track : Particle {
    // constructors //
    Track() = default;
    Track(const Vector<7>& p, const SymMatrix<7>& cov, int charge, size_t idx_track) : Particle{p, cov, charge}, idx{idx_track} {}

    // member vars //
    size_t idx{};
};

struct alignas(32) V0 : Particle {
    // constructors //
    V0() = default;
    V0(const Vector<7>& p, const SymMatrix<7>& cov, int charge, size_t idx_v0, size_t idx_v0_neg, size_t idx_v0_pos, int pdg_code_v0_hyp)
        : Particle{p, cov, charge}, idx{idx_v0}, idx_neg{idx_v0_neg}, idx_pos{idx_v0_pos}, pdg_code_hyp{pdg_code_v0_hyp} {}

    // utilities //
    void SetIndices(const std::array<size_t, 3>& indices) {
        idx = indices[0];
        idx_neg = indices[1];
        idx_pos = indices[2];
    }
    KF::Vector<3> PCA_Neg() const { return GetPCA(0).xyz; };
    KF::Vector<3> PCA_Pos() const { return GetPCA(1).xyz; };
    KF::Vector<3> Mom_Neg() const { return GetPCA(0).dir; };
    KF::Vector<3> Mom_Pos() const { return GetPCA(1).dir; };
    double ArmenterosQt() const {
        return Tree2Secondaries::Math::ArmenterosQt(Px(), Py(), Pz(),  //
                                                    Mom_Neg()[0], Mom_Neg()[1], Mom_Neg()[2]);
    }
    double ArmenterosAlpha() const {
        return Tree2Secondaries::Math::ArmenterosAlpha(Px(), Py(), Pz(),                          //
                                                       Mom_Neg()[0], Mom_Neg()[1], Mom_Neg()[2],  //
                                                       Mom_Pos()[0], Mom_Pos()[1], Mom_Pos()[2]);
    }

    // cuts //
    double AbsZ() const { return std::abs(Z()); }
    double AbsEta() const { return std::abs(Eta()); }
    double DCA_Daughters() const { return GetDCA(0, 1); }
    double DCA_Neg_V0() const { return GetDCA(0); }
    double DCA_Pos_V0() const { return GetDCA(1); }
    double AbsArmQtOverAlpha() const { return ArmenterosQt() / std::abs(ArmenterosAlpha()); };
    double DCA_Point(float x, float y, float z) const {
        return Tree2Secondaries::Math::FastDCALineVertex({Px(), Py(), Pz()}, {X(), Y(), Z()}, {x, y, z});
    }
    double CPA_Point(float x, float y, float z) const {
        return Tree2Secondaries::Math::CosinePointingAngle({Px(), Py(), Pz()}, {X(), Y(), Z()}, {x, y, z});
    }

    // member vars //
    size_t idx{};
    size_t idx_neg{};
    size_t idx_pos{};
    double Neg_Energy{};
    double Pos_Energy{};
    int pdg_code_hyp{};  // hypothesis
};

struct alignas(32) Sexaquark : Particle {
    // utilities //
    double E_MinusNucleon() const { return E() - Nucleon_Mass; };
    double Mass_MinusNucleon() const {
        double mass2{E_MinusNucleon() * E_MinusNucleon() - P2()};
        if (mass2 > 0.) return std::sqrt(mass2);
        return Tree2Secondaries::Const::DummyDouble;
    };

    // cuts //
    double AbsRapidity_MinusNucleon() const { return std::abs(std::log((E_MinusNucleon() + Pz()) / (E_MinusNucleon() - Pz())) / 2.); };
    double CPA_Point(float x, float y, float z) const {
        return Tree2Secondaries::Math::CosinePointingAngle({Px(), Py(), Pz()}, {X(), Y(), Z()}, {x, y, z});
    }

    // member vars //
    size_t idx{};
    double Nucleon_Mass{};
};

struct alignas(32) ChannelA : Sexaquark {
    // utilities //
    void SetIndices(const std::array<size_t, 7>& indices) {
        idx = indices[0];
        idx_lambda = indices[1];
        idx_lambda_neg = indices[2];
        idx_lambda_pos = indices[3];
        idx_k0s = indices[4];
        idx_k0s_neg = indices[5];
        idx_k0s_pos = indices[6];
    }
    void SetAdditionalInfo_V0A(const KF::Vector<3>& decay_vtx, double energy) {
        V0A_DecayVtx = decay_vtx;
        V0A_Energy = energy;
    }
    void SetAdditionalInfo_V0B(const KF::Vector<3>& decay_vtx, double energy) {
        V0B_DecayVtx = decay_vtx;
        V0B_Energy = energy;
    }
    KF::Vector<3> PCA_V0A() const { return GetPCA(0).xyz; };
    KF::Vector<3> PCA_V0B() const { return GetPCA(1).xyz; };
    KF::Vector<3> Mom_V0A() const { return GetPCA(0).dir; };
    KF::Vector<3> Mom_V0B() const { return GetPCA(1).dir; };

    // cuts //
    double DecayLength_V0A() const {
        KF::Vector<3> diff{};
        for (int i{0}; i < 3; ++i) diff[i] += V0A_DecayVtx[i] - PCA_V0A()[i];
        return KF::Math::SquaredNorm(diff);
    };
    double DecayLength_V0B() const {
        KF::Vector<3> diff{};
        for (int i{0}; i < 3; ++i) diff[i] += V0B_DecayVtx[i] - PCA_V0B()[i];
        return KF::Math::SquaredNorm(diff);
    };
    double DCA_V0s() const { return GetDCA(0, 1); };
    double DCA_V0A() const { return GetDCA(0); };
    double DCA_V0B() const { return GetDCA(1); };
    double DCA_V0A_Neg() const;  // PENDING: not trivial
    double DCA_V0A_Pos() const;  // PENDING: not trivial
    double DCA_V0B_Neg() const;  // PENDING: not trivial
    double DCA_V0B_Pos() const;  // PENDING: not trivial

    // member vars //
    KF::Vector<3> V0A_DecayVtx{};
    KF::Vector<3> V0B_DecayVtx{};
    double V0A_Energy{};
    double V0B_Energy{};
    size_t idx_lambda{};
    size_t idx_lambda_neg{};
    size_t idx_lambda_pos{};
    size_t idx_k0s{};
    size_t idx_k0s_neg{};
    size_t idx_k0s_pos{};
};

struct alignas(32) ChannelD : Sexaquark {
    // utilities //
    void SetIndices(const std::array<size_t, 5>& indices) {
        idx = indices[0];
        idx_lambda = indices[1];
        idx_lambda_neg = indices[2];
        idx_lambda_pos = indices[3];
        idx_kaon = indices[4];
    }
    void SetAdditionalInfo_V0(const KF::Vector<3>& decay_vtx, double energy) {
        V0_DecayVtx = decay_vtx;
        V0_Energy = energy;
    }
    KF::Vector<3> PCA_V0() const { return GetPCA(0).xyz; };
    KF::Vector<3> PCA_Kaon() const { return GetPCA(1).xyz; };
    KF::Vector<3> Mom_V0() const { return GetPCA(0).dir; };
    KF::Vector<3> Mom_Kaon() const { return GetPCA(1).dir; };

    // cuts //
    double DCA_V0_Kaon() const { return GetDCA(0, 1); };
    double DCA_V0() const { return GetDCA(0); };
    double DCA_Kaon() const { return GetDCA(1); };
    double DCA_V0_Neg() const;       // PENDING: not trivial
    double DCA_V0_Pos() const;       // PENDING: not trivial
    double DCA_V0_Neg_Kaon() const;  // PENDING: not trivial
    double DCA_V0_Pos_Kaon() const;  // PENDING: not trivial

    // member vars //
    KF::Vector<3> V0_DecayVtx{};
    double V0_Energy{};
    double Kaon_Energy{};
    size_t idx_lambda{};
    size_t idx_lambda_neg{};
    size_t idx_lambda_pos{};
    size_t idx_kaon{};
};

struct alignas(32) ChannelE : Sexaquark {
    size_t idx_lambda{};
    size_t idx_lambda_neg{};
    size_t idx_lambda_pos{};
    size_t idx_kaon{};
    size_t idx_pim{};
    size_t idx_pip{};
};

struct alignas(32) ChannelH : Sexaquark {
    size_t idx_kaon1{};
    size_t idx_kaon2{};
};

// Utitilies //

inline KF::Vector<6> PackParams(const Tree2Secondaries::Events::Tracks& soa, size_t esd_idx) {
    return {soa.X->at(esd_idx),  soa.Y->at(esd_idx),  soa.Z->at(esd_idx),  //
            soa.Px->at(esd_idx), soa.Py->at(esd_idx), soa.Pz->at(esd_idx)};
}

inline KF::Vector<7> UnpackParams(const Tree2Secondaries::PackedEvents::Particle& sov, size_t entry) {
    return {sov.X->at(entry),  sov.Y->at(entry),  sov.Z->at(entry),  //
            sov.Px->at(entry), sov.Py->at(entry), sov.Pz->at(entry), sov.Pz->at(entry)};
}

inline KF::SymMatrix<7> UnpackCovMatrix(const Tree2Secondaries::PackedEvents::Particle& sov, size_t entry) {
    return {sov.Sigma.X2->at(entry),                                                                                                            //
            sov.Sigma.XY->at(entry),  sov.Sigma.Y2->at(entry),                                                                                  //
            sov.Sigma.XZ->at(entry),  sov.Sigma.YZ->at(entry),  sov.Sigma.Z2->at(entry),                                                        //
            sov.Sigma.XPx->at(entry), sov.Sigma.YPx->at(entry), sov.Sigma.ZPx->at(entry), sov.Sigma.Px2->at(entry),                             //
            sov.Sigma.XPy->at(entry), sov.Sigma.YPy->at(entry), sov.Sigma.ZPy->at(entry), sov.Sigma.PxPy->at(entry), sov.Sigma.Py2->at(entry),  //
            sov.Sigma.XPz->at(entry), sov.Sigma.YPz->at(entry), sov.Sigma.ZPz->at(entry), sov.Sigma.PxPz->at(entry), sov.Sigma.PyPz->at(entry),
            sov.Sigma.Pz2->at(entry),  //
            sov.Sigma.XE->at(entry),  sov.Sigma.YE->at(entry),  sov.Sigma.ZE->at(entry),  sov.Sigma.PxE->at(entry),  sov.Sigma.PyE->at(entry),
            sov.Sigma.PzE->at(entry), sov.Sigma.E2->at(entry)};
}

inline Particle CreateParticle(const KF::Vector<6>& kf_params, const std::array<float, 5>& alice_params, const std::array<float, 15>& alice_cov,
                               float alpha, int charge, double mass) {

    double pt{1. / std::abs(alice_params[4]) * std::abs(charge)};
    double cs{std::cos(alpha)};
    double sn{std::sin(alpha)};
    double r{std::sqrt((1. - alice_params[2]) * (1. + alice_params[2]))};

    double m00{-sn};
    double m10{cs};
    double m23{-pt * (sn + alice_params[2] * cs / r)};
    double m43{-pt * pt * (r * cs - alice_params[2] * sn)};
    double m24{pt * (cs - alice_params[2] * sn / r)};
    double m44{-pt * pt * (r * sn + alice_params[2] * cs)};
    double m35{pt};
    double m45{-pt * pt * alice_params[3]};

    m43 *= charge;
    m44 *= charge;
    m45 *= charge;

    SymMatrix<6> kf_cov{alice_cov[0] * m00 * m00,
                        alice_cov[0] * m00 * m10,
                        alice_cov[0] * m10 * m10,
                        alice_cov[1] * m00,
                        alice_cov[1] * m10,
                        alice_cov[2],
                        m00 * (alice_cov[3] * m23 + alice_cov[10] * m43),
                        m10 * (alice_cov[3] * m23 + alice_cov[10] * m43),
                        alice_cov[4] * m23 + alice_cov[11] * m43,
                        m23 * (alice_cov[5] * m23 + alice_cov[12] * m43) + m43 * (alice_cov[12] * m23 + alice_cov[14] * m43),
                        m00 * (alice_cov[3] * m24 + alice_cov[10] * m44),
                        m10 * (alice_cov[3] * m24 + alice_cov[10] * m44),
                        alice_cov[4] * m24 + alice_cov[11] * m44,
                        m23 * (alice_cov[5] * m24 + alice_cov[12] * m44) + m43 * (alice_cov[12] * m24 + alice_cov[14] * m44),
                        m24 * (alice_cov[5] * m24 + alice_cov[12] * m44) + m44 * (alice_cov[12] * m24 + alice_cov[14] * m44),
                        m00 * (alice_cov[6] * m35 + alice_cov[10] * m45),
                        m10 * (alice_cov[6] * m35 + alice_cov[10] * m45),
                        alice_cov[7] * m35 + alice_cov[11] * m45,
                        m23 * (alice_cov[8] * m35 + alice_cov[12] * m45) + m43 * (alice_cov[13] * m35 + alice_cov[14] * m45),
                        m24 * (alice_cov[8] * m35 + alice_cov[12] * m45) + m44 * (alice_cov[13] * m35 + alice_cov[14] * m45),
                        m35 * (alice_cov[9] * m35 + alice_cov[13] * m45) + m45 * (alice_cov[13] * m35 + alice_cov[14] * m45)};

    return {kf_params, kf_cov, charge, mass};
}

}  // namespace KF

namespace ALICE {

inline std::array<float, 5> PackParams(const Tree2Secondaries::Events::Tracks& soa, size_t esd_idx) {
    return {soa.Y->at(esd_idx), soa.Z->at(esd_idx), soa.Snp->at(esd_idx),  //
            soa.Tgl->at(esd_idx), soa.Signed1Pt->at(esd_idx)};
};

inline std::array<float, 15> PackCovMatrix(const Tree2Secondaries::Events::Tracks& soa, size_t esd_idx) {
    return {soa.SigmaY2->at(esd_idx),                                                                                          //
            soa.SigmaZY->at(esd_idx),   soa.SigmaZ2->at(esd_idx),                                                              //
            soa.SigmaSnpY->at(esd_idx), soa.SigmaSnpZ->at(esd_idx), soa.SigmaSnp2->at(esd_idx),                                //
            soa.SigmaTglY->at(esd_idx), soa.SigmaTglZ->at(esd_idx), soa.SigmaTglSnp->at(esd_idx), soa.SigmaTgl2->at(esd_idx),  //
            soa.Sigma1PtY->at(esd_idx), soa.Sigma1PtZ->at(esd_idx), soa.Sigma1PtSnp->at(esd_idx), soa.Sigma1PtTgl->at(esd_idx),
            soa.Sigma1Pt2->at(esd_idx)};
};

}  // namespace ALICE

#endif  // T2S_KF_WRAPPER_HXX
