#ifndef T2S_DF_SEXAQUARK_HXX
#define T2S_DF_SEXAQUARK_HXX

#include "DataFormats/DataFormats.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::DF {

struct alignas(T2S_SIMD_ALIGN) Sexaquark : Flat::State {
    // event properties
    unsigned int RunNumber{};
    unsigned int DirNumber{};
    unsigned int DirNumberB{};  // note: not used when analyzing real data
    unsigned int EventNumber{};
    float MagneticField{};
    // -- primary vertex
    float PV_Xv{};
    float PV_Yv{};
    float PV_Zv{};
    // fit info
    float Chi2NDF{};
    // extra kinematics
    float E_MinusNucleon{};
};

struct alignas(T2S_SIMD_ALIGN) MC_Sexaquark {
    Flat::LorentzVector Before{};
    Flat::LorentzVector After{};
    Flat::LorentzVector Nucleon{};
    // secondary vertex
    float X{};
    float Y{};
    float Z{};
    // event properties
    float PV_Xv{};
    float PV_Yv{};
    float PV_Zv{};
    // reaction id + flags
    int ReactionID{};
    bool IsSignal{};
    bool IsHybrid{};
};

}  // namespace Tree2Secondaries::DF

#endif  // T2S_DF_SEXAQUARK_HXX
