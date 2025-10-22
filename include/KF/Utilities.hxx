#ifndef T2S_KF_UTILS_HXX
#define T2S_KF_UTILS_HXX

#include <cmath>

#include <KFParticle.hxx>
#include <KFParticle_Math.hxx>

#include "DataFormats/DataFormats.hxx"
#include "DataFormats/PackedEvents.hxx"
#include "KF/Extensions.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::KF {

inline ::KF::Vector<6> IntoKF_States_NoE(const DF::SOV::States_NoE& df, int iESD) {
    return {df.X->at(iESD),  df.Y->at(iESD),  df.Z->at(iESD),  //
            df.Px->at(iESD), df.Py->at(iESD), df.Pz->at(iESD)};
}

inline ::KF::SymMatrix<6> IntoKF_CovMatrices_NoE(const DF::SOV::CovMatrices_NoE& df, int iESD) {
    return {df.SigmaX2->at(iESD),                                                                         //
            df.SigmaXY->at(iESD),   df.SigmaY2->at(iESD),                                                 //
            df.SigmaXZ->at(iESD),   df.SigmaYZ->at(iESD),  df.SigmaZ2->at(iESD),                          //
            df.SigmaXPx->at(iESD),  df.SigmaYPx->at(iESD), df.SigmaZPx->at(iESD), df.SigmaPx2->at(iESD),  //
            df.SigmaXPy->at(iESD),  df.SigmaYPy->at(iESD), df.SigmaZPy->at(iESD), df.SigmaPxPy->at(iESD),
            df.SigmaPy2->at(iESD),  //
            df.SigmaXPz->at(iESD),  df.SigmaYPz->at(iESD), df.SigmaZPz->at(iESD), df.SigmaPxPz->at(iESD),
            df.SigmaPyPz->at(iESD), df.SigmaPz2->at(iESD)};
}

inline ::KF::Vector<7> IntoKF_States(const DF::SOV::States& df, int idx) {
    return {df.X->at(idx),  df.Y->at(idx),  df.Z->at(idx),  //
            df.Px->at(idx), df.Py->at(idx), df.Pz->at(idx), df.Energy->at(idx)};
}

inline ::KF::SymMatrix<7> IntoKF_CovMatrices(const DF::SOV::CovMatrices& df, int idx) {
    return {df.SigmaX2->at(idx),                                                                     //
            df.SigmaXY->at(idx),  df.SigmaY2->at(idx),                                               //
            df.SigmaXZ->at(idx),  df.SigmaYZ->at(idx),  df.SigmaZ2->at(idx),                         //
            df.SigmaXPx->at(idx), df.SigmaYPx->at(idx), df.SigmaZPx->at(idx), df.SigmaPx2->at(idx),  //
            df.SigmaXPy->at(idx), df.SigmaYPy->at(idx), df.SigmaZPy->at(idx), df.SigmaPxPy->at(idx),
            df.SigmaPy2->at(idx),  //
            df.SigmaXPz->at(idx), df.SigmaYPz->at(idx), df.SigmaZPz->at(idx), df.SigmaPxPz->at(idx), df.SigmaPyPz->at(idx),
            df.SigmaPz2->at(idx),  //
            df.SigmaXE->at(idx),  df.SigmaYE->at(idx),  df.SigmaZE->at(idx),  df.SigmaPxE->at(idx),  df.SigmaPyE->at(idx),
            df.SigmaPzE->at(idx), df.SigmaE2->at(idx)};
}

inline KF::Track UnpackTrack(const DF::Packed::Tracks& df, int charge, double mass, int idx) {

    auto p = IntoKF_States_NoE(df, idx);
    auto cov = IntoKF_CovMatrices_NoE(df, idx);

    return {p, cov, charge, mass, idx};
}

inline KF::V0 UnpackV0(const DF::Packed::V0s& df, int idx, EParticle v0_pid, EParticle neg_pid, EParticle pos_pid) {

    auto param_v0 = IntoKF_States(df, idx);
    auto cov_v0 = IntoKF_CovMatrices(df, idx);

    KF::Track neg{UnpackTrack(df.Neg, -1, Particle::Mass[neg_pid], df.Neg.Entry->at(idx))};
    KF::Track pos{UnpackTrack(df.Pos, +1, Particle::Mass[pos_pid], df.Pos.Entry->at(idx))};

    return {idx, v0_pid, param_v0, cov_v0, neg, pos};
}

}  // namespace Tree2Secondaries::KF

#endif  // T2S_KF_UTILS_HXX
