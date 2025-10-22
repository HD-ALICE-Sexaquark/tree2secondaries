#ifndef T2S_DF_CHANNEL_D_HXX
#define T2S_DF_CHANNEL_D_HXX

#include "DataFormats/Sexaquark.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::DF {

struct alignas(T2S_SIMD_ALIGN) ChannelD : Sexaquark {
    // V0 : anti-lambda, lambda
    Flat::State V0_atPCA{};
    Flat::State V0_atDecay{};
    // Kaon : positive kaon, negative kaon
    Flat::State Kaon_atPCA{};
    // V0_Neg : anti-proton from anti-lambda, pi-minus from lambda
    Flat::State V0_Neg_atPCA{};
    Flat::State V0_Neg_atV0{};
    // V0_Pos : pi-plus from anti-lambda, proton from lambda
    Flat::State V0_Pos_atPCA{};
    Flat::State V0_Pos_atV0{};
    // indices
    int V0_Entry{};
    int Kaon_Entry{};
    int V0_Neg_Entry{};
    int V0_Pos_Entry{};
};

struct alignas(T2S_SIMD_ALIGN) MC_ChannelD : MC_Sexaquark {
    // reaction products
    Flat::LorentzVector V0_atSV{};
    Flat::LorentzVector Kaon_atSV{};
    // V0 info
    int V0_Entry{};
    int V0_PdgCode{};
    int V0_Mother_Entry{};
    int V0_Mother_PdgCode{};
    int V0_ReactionID{};
    // V0 daughters
    int V0_Neg_Entry{};
    int V0_Neg_PdgCode{};
    int V0_Pos_Entry{};
    int V0_Pos_PdgCode{};
    // Kaon info
    int Kaon_Entry{};
    int Kaon_PdgCode{};
    int Kaon_Mother_Entry{};
    int Kaon_Mother_PdgCode{};
    int Kaon_GrandMother_Entry{};
    int Kaon_GrandMother_PdgCode{};
    int Kaon_ReactionID{};
    // V0 flags
    bool V0_IsTrue{};
    bool V0_IsSignal{};
    bool V0_IsSecondary{};
    bool V0_IsHybrid{};
    // Kaon flags
    bool Kaon_IsTrue{};
    bool Kaon_IsSignal{};
    bool Kaon_IsSecondary{};
};

}  // namespace Tree2Secondaries::DF

#endif  // T2S_DF_CHANNEL_D_HXX
