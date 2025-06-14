#ifndef T2S_KF_UTILS_HXX
#define T2S_KF_UTILS_HXX

#include <array>
#include <cmath>
#include <cstddef>

#include <KFParticle.hxx>
#include <KFParticle_Math.hxx>

#include "DataFormats/Events.hxx"
#include "DataFormats/PackedEvents.hxx"

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

inline KF::Vector<6> PackParams(const Tree2Secondaries::Events::Tracks& soa, size_t esd_idx) {
    return {soa.X->at(esd_idx),  soa.Y->at(esd_idx),  soa.Z->at(esd_idx),  //
            soa.Px->at(esd_idx), soa.Py->at(esd_idx), soa.Pz->at(esd_idx)};
}

inline std::array<float, 5> PackParams_ALICE(const Tree2Secondaries::Events::Tracks& soa, size_t esd_idx) {
    return {soa.Y->at(esd_idx), soa.Z->at(esd_idx), soa.Snp->at(esd_idx),  //
            soa.Tgl->at(esd_idx), soa.Signed1Pt->at(esd_idx)};
};

inline std::array<float, 15> PackCovMatrix_ALICE(const Tree2Secondaries::Events::Tracks& soa, size_t esd_idx) {
    return {soa.SigmaY2->at(esd_idx),                                                                                          //
            soa.SigmaZY->at(esd_idx),   soa.SigmaZ2->at(esd_idx),                                                              //
            soa.SigmaSnpY->at(esd_idx), soa.SigmaSnpZ->at(esd_idx), soa.SigmaSnp2->at(esd_idx),                                //
            soa.SigmaTglY->at(esd_idx), soa.SigmaTglZ->at(esd_idx), soa.SigmaTglSnp->at(esd_idx), soa.SigmaTgl2->at(esd_idx),  //
            soa.Sigma1PtY->at(esd_idx), soa.Sigma1PtZ->at(esd_idx), soa.Sigma1PtSnp->at(esd_idx), soa.Sigma1PtTgl->at(esd_idx),
            soa.Sigma1Pt2->at(esd_idx)};
};

inline KF::Vector<7> UnpackParams(const Tree2Secondaries::PackedEvents::Particle& sov, size_t esd_idx) {
    return {sov.X->at(esd_idx),  sov.Y->at(esd_idx),  sov.Z->at(esd_idx),  //
            sov.Px->at(esd_idx), sov.Py->at(esd_idx), sov.Pz->at(esd_idx), sov.Pz->at(esd_idx)};
}

inline KF::SymMatrix<7> UnpackCovMatrix(const Tree2Secondaries::PackedEvents::Particle& sov, size_t esd_idx) {
    return {sov.Sigma.X2->at(esd_idx),                                                                                        //
            sov.Sigma.XY->at(esd_idx),   sov.Sigma.Y2->at(esd_idx),                                                           //
            sov.Sigma.XZ->at(esd_idx),   sov.Sigma.YZ->at(esd_idx),  sov.Sigma.Z2->at(esd_idx),                               //
            sov.Sigma.XPx->at(esd_idx),  sov.Sigma.YPx->at(esd_idx), sov.Sigma.ZPx->at(esd_idx), sov.Sigma.Px2->at(esd_idx),  //
            sov.Sigma.XPy->at(esd_idx),  sov.Sigma.YPy->at(esd_idx), sov.Sigma.ZPy->at(esd_idx), sov.Sigma.PxPy->at(esd_idx),
            sov.Sigma.Py2->at(esd_idx),  //
            sov.Sigma.XPz->at(esd_idx),  sov.Sigma.YPz->at(esd_idx), sov.Sigma.ZPz->at(esd_idx), sov.Sigma.PxPz->at(esd_idx),
            sov.Sigma.PyPz->at(esd_idx),
            sov.Sigma.Pz2->at(esd_idx),  //
            sov.Sigma.XE->at(esd_idx),   sov.Sigma.YE->at(esd_idx),  sov.Sigma.ZE->at(esd_idx),  sov.Sigma.PxE->at(esd_idx),
            sov.Sigma.PyE->at(esd_idx),  sov.Sigma.PzE->at(esd_idx), sov.Sigma.E2->at(esd_idx)};
}

}  // namespace KF

#endif  // T2S_KF_UTILS_HXX
