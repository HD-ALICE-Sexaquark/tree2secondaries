#ifndef T2S_DF_CHANNEL_D_HXX
#define T2S_DF_CHANNEL_D_HXX

#include <TTree.h>

#include "App/Utilities.hxx"
#include "DataFormats/DataFormats.hxx"
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

    void CreateBranches_ChannelD(TTree* tree, bool is_mc) {
        CreateBranches_Sexaquark(tree, is_mc);
        // V0
        V0_atPCA.CreateBranches_State(tree, "V0_atPCA");
        V0_atDecay.CreateBranches_State(tree, "V0_atDecay");
        // Kaon
        Kaon_atPCA.CreateBranches_State(tree, "Kaon_atPCA");
        // V0 Neg
        V0_Neg_atPCA.CreateBranches_State(tree, "V0_Neg_atPCA");
        V0_Neg_atV0.CreateBranches_State(tree, "V0_Neg_atV0");
        // V0 Pos
        V0_Pos_atPCA.CreateBranches_State(tree, "V0_Pos_atPCA");
        V0_Pos_atV0.CreateBranches_State(tree, "V0_Pos_atV0");
        // indices
        Utils::CreateBranch(tree, "V0_Entry", &V0_Entry);
        Utils::CreateBranch(tree, "Kaon_Entry", &Kaon_Entry);
        Utils::CreateBranch(tree, "V0_Neg_Entry", &V0_Neg_Entry);
        Utils::CreateBranch(tree, "V0_Pos_Entry", &V0_Pos_Entry);
    }
};

struct alignas(T2S_SIMD_ALIGN) MC_ChannelD : MC_Sexaquark {
    Flat::TrueInfo V0{};
    Flat::TrueInfo Kaon{};
    Flat::TrueInfo V0_Neg{};
    Flat::TrueInfo V0_Pos{};
    Flat::Coordinates V0_atDecay{};
    int Kaon_GrandMother_Entry{};
    int Kaon_GrandMother_PdgCode{};
    bool V0_IsHybrid{};

    void CreateBranches_MC_ChannelD(TTree* tree) {
        CreateBranches_MC_Sexaquark(tree);
        V0.CreateBranches_TrueInfo(tree, "V0A");
        Kaon.CreateBranches_TrueInfo(tree, "V0B");
        V0_Neg.CreateBranches_TrueInfo(tree, "V0_Neg");
        V0_Pos.CreateBranches_TrueInfo(tree, "V0_Pos");
        V0_atDecay.CreateBranches_Coordinates(tree, "V0_atDecay");
        Utils::CreateBranch(tree, "Kaon_GrandMother_Entry", &Kaon_GrandMother_Entry);
        Utils::CreateBranch(tree, "Kaon_GrandMother_PdgCode", &Kaon_GrandMother_PdgCode);
        Utils::CreateBranch(tree, "V0_IsHybrid", &V0_IsHybrid);
    }
};

}  // namespace Tree2Secondaries::DF

#endif  // T2S_DF_CHANNEL_D_HXX
