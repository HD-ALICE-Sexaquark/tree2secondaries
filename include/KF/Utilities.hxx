#ifndef T2S_KF_UTILS_HXX
#define T2S_KF_UTILS_HXX

#include <cmath>

#include <KFParticle.hxx>
#include <KFParticle_Math.hxx>

#include "DataFormats/Events.hxx"
#include "DataFormats/PackedEvents.hxx"
#include "KF/Extensions.hxx"
#include "Math/Constants.hxx"

// Utitilies //

namespace KF {

// # From `Events` //

inline KF::Vector<6> Pack_XYZPxPyPz(const Tree2Secondaries::Events::Tracks& soa, int iESD) {
    return {soa.X->at(iESD),  soa.Y->at(iESD),  soa.Z->at(iESD),  //
            soa.Px->at(iESD), soa.Py->at(iESD), soa.Pz->at(iESD)};
}

inline KF::SymMatrix<6> Pack_CovXYZPxPyPz(const Tree2Secondaries::Events::Tracks& soa, int iESD) {
    return {soa.SigmaX2->at(iESD),                                                                            //
            soa.SigmaXY->at(iESD),   soa.SigmaY2->at(iESD),                                                   //
            soa.SigmaXZ->at(iESD),   soa.SigmaYZ->at(iESD),  soa.SigmaZ2->at(iESD),                           //
            soa.SigmaXPx->at(iESD),  soa.SigmaYPx->at(iESD), soa.SigmaZPx->at(iESD), soa.SigmaPx2->at(iESD),  //
            soa.SigmaXPy->at(iESD),  soa.SigmaYPy->at(iESD), soa.SigmaZPy->at(iESD), soa.SigmaPxPy->at(iESD),
            soa.SigmaPy2->at(iESD),  //
            soa.SigmaXPz->at(iESD),  soa.SigmaYPz->at(iESD), soa.SigmaZPz->at(iESD), soa.SigmaPxPz->at(iESD),
            soa.SigmaPyPz->at(iESD), soa.SigmaPz2->at(iESD)};
}

// # From `PackedEvents` //

inline KF::Vector<7> Unpack_XYZPxPyPzE(const Tree2Secondaries::PackedEvents::State& sov, int idx) {
    return {sov.X->at(idx),  sov.Y->at(idx),  sov.Z->at(idx),  //
            sov.Px->at(idx), sov.Py->at(idx), sov.Pz->at(idx), sov.E->at(idx)};
}

inline KF::SymMatrix<7> Unpack_CovXYZPxPyPzE(const Tree2Secondaries::PackedEvents::Particle& sov, int idx) {
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

inline KF::V0 UnpackV0(const Tree2Secondaries::PackedEvents::V0s& sov, int idx, Tree2Secondaries::EParticle v0_hypothesis) {

    auto param_v0 = Unpack_XYZPxPyPzE(sov, idx);
    auto cov_v0 = Unpack_CovXYZPxPyPzE(sov, idx);

    auto param_neg = Unpack_XYZPxPyPzE(sov.Neg, idx);
    KF::Track neg{param_neg, {}, -1, sov.Neg.Entry->at(idx)};

    auto param_pos = Unpack_XYZPxPyPzE(sov.Pos, idx);
    KF::Track pos{param_pos, {}, +1, sov.Pos.Entry->at(idx)};

    return {idx, v0_hypothesis, param_v0, cov_v0, neg, pos};
}

inline KF::Track UnpackTrack(const Tree2Secondaries::PackedEvents::Tracks& sov, int idx, int charge) {

    auto param = Unpack_XYZPxPyPzE(sov, idx);
    auto cov = Unpack_CovXYZPxPyPzE(sov, idx);

    return {param, cov, charge, idx};
}

}  // namespace KF

#endif  // T2S_KF_UTILS_HXX
