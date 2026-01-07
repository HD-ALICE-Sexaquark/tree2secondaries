#pragma once

#include <TTree.h>

#include "Math/Constants.hxx"
#include "Storage/Flat/BaseFlat.hxx"
#include "Storage/Flat/FlatSexaquark.hxx"
#include "Storage/Flat/FlatV0.hxx"

namespace Tree2Secondaries::Storage::Flat {

// `Flat::Sexaquark` +
// `V0A` (`Flat::V0`) + `V0B` (`Flat::V0`) +
// `V0A_atPCA` (`Flat::State_NoE`) + `V0B_atPCA` (`Flat::State_NoE`).
struct alignas(T2S_SIMD_ALIGN) ChannelA : Flat::Sexaquark {
    Flat::V0 V0A{};
    Flat::V0 V0B{};
    Flat::State_NoE V0A_atPCA{};
    Flat::State_NoE V0B_atPCA{};

    void CreateBranches_FlatChannelA(TTree* tree, bool is_mc) {
        CreateBranches_FlatSexaquark(tree, is_mc);
        V0A.CreateBranches_FlatV0(tree, "V0A");
        V0B.CreateBranches_FlatV0(tree, "V0B");
        V0A_atPCA.CreateBranches_FlatState_NoE(tree, "V0A_atPCA");
        V0B_atPCA.CreateBranches_FlatState_NoE(tree, "V0B_atPCA");
    }
};

// `Flat::MC_Sexaquark` +
// `V0A` (`Flat::MC_V0`) + `V0B` (`Flat::MC_V0`).
struct alignas(T2S_SIMD_ALIGN) MC_ChannelA : Flat::MC_Sexaquark {
    Flat::MC_V0 V0A{};
    Flat::MC_V0 V0B{};

    void CreateBranches_FlatMC_ChannelA(TTree* tree) {
        CreateBranches_FlatMC_Sexaquark(tree);
        V0A.CreateBranches_FlatMC_V0(tree, "V0A");
        V0B.CreateBranches_FlatMC_V0(tree, "V0B");
    }
};

}  // namespace Tree2Secondaries::Storage::Flat
