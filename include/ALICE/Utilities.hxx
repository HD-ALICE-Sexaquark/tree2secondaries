#ifndef ALICE_UTILS_HXX
#define ALICE_UTILS_HXX

#include <array>

#include "DataFormats/Events.hxx"

namespace ALICE {

template <typename T>
inline std::array<T, 5> PackParams(const Tree2Secondaries::Events::Tracks& soa, int esd_idx) {
    return {static_cast<T>(soa.Y->at(esd_idx)), static_cast<T>(soa.Z->at(esd_idx)), static_cast<T>(soa.Snp->at(esd_idx)),  //
            static_cast<T>(soa.Tgl->at(esd_idx)), static_cast<T>(soa.Signed1Pt->at(esd_idx))};
};

template <typename T>
inline std::array<T, 15> PackCovMatrix(const Tree2Secondaries::Events::Tracks& soa, int esd_idx) {
    return {static_cast<T>(soa.SigmaY2->at(esd_idx)),                                                                                              //
            static_cast<T>(soa.SigmaZY->at(esd_idx)),     static_cast<T>(soa.SigmaZ2->at(esd_idx)),                                                //
            static_cast<T>(soa.SigmaSnpY->at(esd_idx)),   static_cast<T>(soa.SigmaSnpZ->at(esd_idx)), static_cast<T>(soa.SigmaSnp2->at(esd_idx)),  //
            static_cast<T>(soa.SigmaTglY->at(esd_idx)),   static_cast<T>(soa.SigmaTglZ->at(esd_idx)), static_cast<T>(soa.SigmaTglSnp->at(esd_idx)),
            static_cast<T>(soa.SigmaTgl2->at(esd_idx)),  //
            static_cast<T>(soa.Sigma1PtY->at(esd_idx)),   static_cast<T>(soa.Sigma1PtZ->at(esd_idx)), static_cast<T>(soa.Sigma1PtSnp->at(esd_idx)),
            static_cast<T>(soa.Sigma1PtTgl->at(esd_idx)), static_cast<T>(soa.Sigma1Pt2->at(esd_idx))};
};

}  // namespace ALICE

#endif  // ALICE_UTILS_HXX
