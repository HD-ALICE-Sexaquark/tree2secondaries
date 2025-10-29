#pragma once

#include <cmath>

#include <KFParticle.hxx>
#include <KFParticle_Math.hxx>

#include "DataFormats/Events.hxx"
#include "DataFormats/PackedEvents.hxx"
#include "Fit/Track.hxx"
#include "Fit/V0.hxx"

namespace Tree2Secondaries::Fit {

inline KF::Vector<6> IntoKF_States_NoE(const DF::SOV::States_NoE& df, int idx) {
    return {df.X->at(idx),  df.Y->at(idx),  df.Z->at(idx),  //
            df.Px->at(idx), df.Py->at(idx), df.Pz->at(idx)};
}

inline KF::SymMatrix<6> IntoKF_CovMatrices_NoE(const DF::SOV::CovMatrices_NoE& df, int idx) {
    return {df.SigmaX2->at(idx),                                                                     //
            df.SigmaXY->at(idx),  df.SigmaY2->at(idx),                                               //
            df.SigmaXZ->at(idx),  df.SigmaYZ->at(idx),  df.SigmaZ2->at(idx),                         //
            df.SigmaXPx->at(idx), df.SigmaYPx->at(idx), df.SigmaZPx->at(idx), df.SigmaPx2->at(idx),  //
            df.SigmaXPy->at(idx), df.SigmaYPy->at(idx), df.SigmaZPy->at(idx), df.SigmaPxPy->at(idx),
            df.SigmaPy2->at(idx),  //
            df.SigmaXPz->at(idx), df.SigmaYPz->at(idx), df.SigmaZPz->at(idx), df.SigmaPxPz->at(idx), df.SigmaPyPz->at(idx), df.SigmaPz2->at(idx)};
}

inline Fit::Track CreateTrack(const DF::Events::Tracks& df, int idx, int charge, double mass) {

    KF::Vector<6> params = IntoKF_States_NoE(df, idx);
    KF::SymMatrix<6> cov = IntoKF_CovMatrices_NoE(df, idx);

    return {params, cov, charge, mass, idx};
}

inline KF::Vector<7> IntoKF_States(const DF::SOV::States& df, int idx) {
    return {df.X->at(idx),  df.Y->at(idx),  df.Z->at(idx),  //
            df.Px->at(idx), df.Py->at(idx), df.Pz->at(idx), df.Energy->at(idx)};
}

inline KF::SymMatrix<7> IntoKF_CovMatrices(const DF::SOV::CovMatrices& df, int idx) {
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

inline Fit::Track UnpackTrack(const DF::Packed::Tracks& df, int idx, int charge, double mass) {

    auto p = IntoKF_States_NoE(df, idx);
    auto cov = IntoKF_CovMatrices_NoE(df, idx);

    return {p, cov, charge, mass, idx};
}

inline Fit::V0 UnpackV0(const DF::Packed::V0s& df, int idx, const Fit::Track& neg, const Fit::Track& pos) {

    auto param_v0 = IntoKF_States(df, idx);
    auto cov_v0 = IntoKF_CovMatrices(df, idx);

    return {idx, param_v0, cov_v0, neg, pos};
}

}  // namespace Tree2Secondaries::Fit
