#ifndef T2S_DF_INJECTED_HXX
#define T2S_DF_INJECTED_HXX

#include <vector>

#include "DataFormats/DataFormats.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::DF {

namespace SOV {
struct alignas(T2S_SIMD_ALIGN) Injected {
    std::vector<int> *ReactionID{nullptr};
    std::vector<float> *X{nullptr};  // NOTE: not used by packager
    std::vector<float> *Y{nullptr};  // NOTE: not used by packager
    std::vector<float> *Z{nullptr};  // NOTE: not used by packager
    std::vector<float> *Px{nullptr};
    std::vector<float> *Py{nullptr};
    std::vector<float> *Pz{nullptr};
    std::vector<float> *Nucleon_Px{nullptr};
    std::vector<float> *Nucleon_Py{nullptr};
    std::vector<float> *Nucleon_Pz{nullptr};
    void Clear_Injected() {
        ReactionID->clear();
        X->clear();
        Y->clear();
        Z->clear();
        Px->clear();
        Py->clear();
        Pz->clear();
        Nucleon_Px->clear();
        Nucleon_Py->clear();
        Nucleon_Pz->clear();
    }
};
}  // namespace SOV

namespace Flat {
struct alignas(T2S_SIMD_ALIGN) Injected : Flat::State {
    // event properties
    unsigned int RunNumber{};
    unsigned int DirNumber{};
    unsigned int EventNumber{};
    // reaction id
    int ReactionID{};
};
}  // namespace Flat

}  // namespace Tree2Secondaries::DF

#endif  // T2S_DF_INJECTED_HXX
