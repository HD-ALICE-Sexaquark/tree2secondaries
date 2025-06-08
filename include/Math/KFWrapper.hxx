#ifndef T2S_KF_WRAPPER_HXX
#define T2S_KF_WRAPPER_HXX

#include <array>
#include <cmath>
#include <cstddef>

#include "Math/Math.hxx"
#include "Structures/Events.hxx"

#include "KFParticle.hxx"
#include "KFParticle_Math.hxx"

namespace KF {

struct alignas(32) Track : Particle {
    size_t idx{};
};

struct alignas(32) V0 : Particle {
    size_t idx{};
    size_t idx_neg{};
    size_t idx_pos{};
    int pdg_code_hyp{};  // hypothesis

    double DCA_Daughters() const { return GetDCA(0, 1); }
    double DCA_Neg_V0() const { return GetDCA(0); }
    double DCA_Pos_V0() const { return GetDCA(1); }
    KF::Vector<3> PCA_Neg() const { return GetPCA(0).xyz; };
    KF::Vector<3> PCA_Pos() const { return GetPCA(1).xyz; };
    KF::Vector<3> Mom_Neg() const { return GetPCA(0).dir; };
    KF::Vector<3> Mom_Pos() const { return GetPCA(1).dir; };
    double AbsZ() const { return std::abs(Z()); }
    double AbsEta() const { return std::abs(Eta()); }
    double ArmenterosQt() const {
        return Tree2Secondaries::Math::ArmenterosQt(Px(), Py(), Pz(),  //
                                                    Mom_Neg()[0], Mom_Neg()[1], Mom_Neg()[2]);
    }
    double ArmenterosAlpha() const {
        return Tree2Secondaries::Math::ArmenterosAlpha(Px(), Py(), Pz(),                          //
                                                       Mom_Neg()[0], Mom_Neg()[1], Mom_Neg()[2],  //
                                                       Mom_Pos()[0], Mom_Pos()[1], Mom_Pos()[2]);
    }
    double AbsArmQtOverAlpha() const { return ArmenterosQt() / std::abs(ArmenterosAlpha()); };
    double DCA_Point(float x, float y, float z) const {
        return Tree2Secondaries::Math::FastDCALineVertex({Px(), Py(), Pz()}, {X(), Y(), Z()}, {x, y, z});
    }
    double CPA_Point(float x, float y, float z) const {
        return Tree2Secondaries::Math::CosinePointingAngle({Px(), Py(), Pz()}, {X(), Y(), Z()}, {x, y, z});
    }
};

struct alignas(32) ChannelA : Particle {
    // PENDING
};

struct alignas(32) ChannelD : Particle {
    // PENDING
};

struct alignas(32) ChannelE : Particle {
    // PENDING
};

struct alignas(32) ChannelH : Particle {
    // PENDING
};

struct Extended : Particle {
    long mc_idx{};
    unsigned long neg_idx{};
    unsigned long pos_idx{};
    unsigned int reaction_id{};
};

inline std::array<double, 6> PackParams(const Tree2Secondaries::Events::Tracks& soa, size_t esd_idx) {
    return {soa.X->at(esd_idx),  soa.Y->at(esd_idx),  soa.Z->at(esd_idx),  //
            soa.Px->at(esd_idx), soa.Py->at(esd_idx), soa.Pz->at(esd_idx)};
}

inline Particle CreateParticle(const std::array<double, 6>& kf_params, const std::array<float, 5>& alice_params,
                               const std::array<float, 15>& alice_cov, float alpha, int charge, double mass) {

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
