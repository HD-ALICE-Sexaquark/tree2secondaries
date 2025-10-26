#ifndef T2S_DF_CHANNEL_A_HXX
#define T2S_DF_CHANNEL_A_HXX

#include <TTree.h>

#include "App/Utilities.hxx"
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

    void CreateBranches_ChannelA(TTree* tree, bool is_mc) {
        CreateBranches_Sexaquark(tree, is_mc);
        // V0A
        V0A_atPCA.CreateBranches_State(tree, "V0A_atPCA");
        V0A_atDecay.CreateBranches_State(tree, "V0A_atDecay");
        // V0B
        V0B_atPCA.CreateBranches_State(tree, "V0B_atPCA");
        V0B_atDecay.CreateBranches_State(tree, "V0B_atDecay");
        // V0A Neg
        V0A_Neg_atPCA.CreateBranches_State(tree, "V0A_Neg_atPCA");
        V0A_Neg_atV0.CreateBranches_State(tree, "V0A_Neg_atV0");
        // V0A Pos
        V0A_Pos_atPCA.CreateBranches_State(tree, "V0A_Pos_atPCA");
        V0A_Pos_atV0.CreateBranches_State(tree, "V0A_Pos_atV0");
        // V0B Neg
        V0B_Neg_atPCA.CreateBranches_State(tree, "V0B_Neg_atPCA");
        V0B_Neg_atV0.CreateBranches_State(tree, "V0B_Neg_atV0");
        // V0B Pos
        V0B_Pos_atPCA.CreateBranches_State(tree, "V0B_Pos_atPCA");
        V0B_Pos_atV0.CreateBranches_State(tree, "V0B_Pos_atV0");
        // indices
        Utils::CreateBranch(tree, "V0A_Entry", &V0A_Entry);
        Utils::CreateBranch(tree, "V0B_Entry", &V0B_Entry);
        Utils::CreateBranch(tree, "V0A_Neg_Entry", &V0A_Neg_Entry);
        Utils::CreateBranch(tree, "V0A_Pos_Entry", &V0A_Pos_Entry);
        Utils::CreateBranch(tree, "V0B_Neg_Entry", &V0B_Neg_Entry);
        Utils::CreateBranch(tree, "V0B_Pos_Entry", &V0B_Pos_Entry);
    }
};

struct alignas(T2S_SIMD_ALIGN) MC_ChannelA : MC_Sexaquark {
    Flat::TrueInfo V0A{};
    Flat::TrueInfo V0B{};
    Flat::TrueInfo V0A_Neg{};
    Flat::TrueInfo V0A_Pos{};
    Flat::TrueInfo V0B_Neg{};
    Flat::TrueInfo V0B_Pos{};
    Flat::Coordinates V0A_atDecay{};
    Flat::Coordinates V0B_atDecay{};
    bool V0A_IsHybrid{};
    bool V0B_IsHybrid{};

    void CreateBranches_MC_ChannelA(TTree* tree) {
        CreateBranches_MC_Sexaquark(tree);
        V0A.CreateBranches_TrueInfo(tree, "V0A");
        V0B.CreateBranches_TrueInfo(tree, "V0B");
        V0A_Neg.CreateBranches_TrueInfo(tree, "V0A_Neg");
        V0A_Pos.CreateBranches_TrueInfo(tree, "V0A_Pos");
        V0B_Neg.CreateBranches_TrueInfo(tree, "V0B_Neg");
        V0B_Pos.CreateBranches_TrueInfo(tree, "V0B_Pos");
        V0A_atDecay.CreateBranches_Coordinates(tree, "V0A_atDecay");
        V0B_atDecay.CreateBranches_Coordinates(tree, "V0B_atDecay");
        Utils::CreateBranch(tree, "V0A_IsHybrid", &V0A_IsHybrid);
        Utils::CreateBranch(tree, "V0B_IsHybrid", &V0B_IsHybrid);
    }
};

}  // namespace Tree2Secondaries::DF

#endif  // T2S_DF_CHANNEL_A_HXX
