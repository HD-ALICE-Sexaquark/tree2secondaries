#ifndef T2S_DF_CHANNEL_A_HXX
#define T2S_DF_CHANNEL_A_HXX

#include "DataFormats/DataFormats.hxx"
#include "DataFormats/Sexaquark.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::DF {

struct alignas(T2S_SIMD_ALIGN) ChannelA : Sexaquark {
    // V0A : anti-lambda, lambda
    Flat::State V0A_atPCA{};
    Flat::State V0A_atDecay{};
    // V0B : kaon-zero-short
    Flat::State V0B_atPCA{};
    Flat::State V0B_atDecay{};
    // V0A Neg : anti-proton from anti-lambda, pi-minus from lambda
    Flat::State V0A_Neg_atPCA{};
    Flat::State V0A_Neg_atV0{};
    // V0A Pos : pi-plus from anti-lambda, proton from lambda
    Flat::State V0A_Pos_atPCA{};
    Flat::State V0A_Pos_atV0{};
    // V0B Neg : pi-minus from kaon-zero-short
    Flat::State V0B_Neg_atPCA{};
    Flat::State V0B_Neg_atV0{};
    // V0B Pos : pi-plus from kaon-zero-short
    Flat::State V0B_Pos_atPCA{};
    Flat::State V0B_Pos_atV0{};
    // indices
    int V0A_Entry{};
    int V0B_Entry{};
    int V0A_Neg_Entry{};
    int V0A_Pos_Entry{};
    int V0B_Neg_Entry{};
    int V0B_Pos_Entry{};
};

struct alignas(T2S_SIMD_ALIGN) MC_ChannelA : MC_Sexaquark {
    Flat::LorentzVector V0A_atSV{};
    Flat::LorentzVector V0B_atSV{};
    // V0A
    int V0A_Entry{};
    int V0A_PdgCode{};
    int V0A_Mother_Entry{};
    int V0A_Mother_PdgCode{};
    int V0A_ReactionID{};
    // V0A daughters
    int V0A_Neg_Entry{};
    int V0A_Neg_PdgCode{};
    int V0A_Pos_Entry{};
    int V0A_Pos_PdgCode{};
    // V0B
    int V0B_Entry{};
    int V0B_PdgCode{};
    int V0B_Mother_Entry{};
    int V0B_Mother_PdgCode{};
    // V0B daughters
    int V0B_Neg_Entry{};
    int V0B_Neg_PdgCode{};
    int V0B_Pos_Entry{};
    int V0B_Pos_PdgCode{};
    int V0B_ReactionID{};
    // V0A flags
    bool V0A_IsTrue{};
    bool V0A_IsSignal{};
    bool V0A_IsSecondary{};
    bool V0A_IsHybrid{};
    // V0B flags
    bool V0B_IsTrue{};
    bool V0B_IsSignal{};
    bool V0B_IsSecondary{};
    bool V0B_IsHybrid{};
};

}  // namespace Tree2Secondaries::DF

#endif  // T2S_DF_CHANNEL_A_HXX
