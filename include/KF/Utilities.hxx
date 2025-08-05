#ifndef T2S_KF_UTILS_HXX
#define T2S_KF_UTILS_HXX

#include <array>
#include <cmath>

#include <KFParticle.hxx>
#include <KFParticle_Math.hxx>

#include "DataFormats/Events.hxx"
#include "DataFormats/PackedEvents.hxx"
#include "KF/Extensions.hxx"
#include "Math/Constants.hxx"

// Utitilies //

namespace KF {

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

inline KF::Vector<6> PackParams(const Tree2Secondaries::Events::Tracks& soa, int esd_idx) {
    return {soa.X->at(esd_idx),  soa.Y->at(esd_idx),  soa.Z->at(esd_idx),  //
            soa.Px->at(esd_idx), soa.Py->at(esd_idx), soa.Pz->at(esd_idx)};
}

inline KF::Vector<7> UnpackParams(const Tree2Secondaries::PackedEvents::State& sov, int idx) {
    return {sov.X->at(idx),  sov.Y->at(idx),  sov.Z->at(idx),  //
            sov.Px->at(idx), sov.Py->at(idx), sov.Pz->at(idx), sov.E->at(idx)};
}

inline KF::SymMatrix<7> UnpackCovMatrix(const Tree2Secondaries::PackedEvents::Particle& sov, int idx) {
    return {sov.Sigma.X2->at(idx),                                                                           //
            sov.Sigma.XY->at(idx),  sov.Sigma.Y2->at(idx),                                                   //
            sov.Sigma.XZ->at(idx),  sov.Sigma.YZ->at(idx),  sov.Sigma.Z2->at(idx),                           //
            sov.Sigma.XPx->at(idx), sov.Sigma.YPx->at(idx), sov.Sigma.ZPx->at(idx), sov.Sigma.Px2->at(idx),  //
            sov.Sigma.XPy->at(idx), sov.Sigma.YPy->at(idx), sov.Sigma.ZPy->at(idx), sov.Sigma.PxPy->at(idx),
            sov.Sigma.Py2->at(idx),  //
            sov.Sigma.XPz->at(idx), sov.Sigma.YPz->at(idx), sov.Sigma.ZPz->at(idx), sov.Sigma.PxPz->at(idx), sov.Sigma.PyPz->at(idx),
            sov.Sigma.Pz2->at(idx),  //
            sov.Sigma.XE->at(idx),  sov.Sigma.YE->at(idx),  sov.Sigma.ZE->at(idx),  sov.Sigma.PxE->at(idx),  sov.Sigma.PyE->at(idx),
            sov.Sigma.PzE->at(idx), sov.Sigma.E2->at(idx)};
}

inline KF::V0 UnpackV0(const Tree2Secondaries::PackedEvents::V0s& sov, int idx, Tree2Secondaries::PdgCode pdg_code_hyp) {

    auto param_v0 = UnpackParams(sov, idx);
    auto cov_v0 = UnpackCovMatrix(sov, idx);

    auto param_neg = UnpackParams(sov.Neg, idx);
    KF::Track neg{param_neg, {}, -1, sov.Neg.Entry->at(idx)};

    auto param_pos = UnpackParams(sov.Pos, idx);
    KF::Track pos{param_pos, {}, +1, sov.Pos.Entry->at(idx)};

    return {idx, pdg_code_hyp, param_v0, cov_v0, neg, pos};
}

inline KF::Track UnpackTrack(const Tree2Secondaries::PackedEvents::Tracks& sov, int idx, int charge) {

    auto param = UnpackParams(sov, idx);
    auto cov = UnpackCovMatrix(sov, idx);

    return {param, cov, charge, idx};
}

}  // namespace KF

#endif  // T2S_KF_UTILS_HXX
