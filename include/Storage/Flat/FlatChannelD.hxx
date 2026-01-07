#pragma once

#include <TTree.h>

#include "Math/Constants.hxx"
#include "Storage/Flat/BaseFlat.hxx"
#include "Storage/Flat/FlatSexaquark.hxx"
#include "Storage/Flat/FlatTrack.hxx"
#include "Storage/Flat/FlatV0.hxx"

namespace Tree2Secondaries::Storage::Flat {

// `Flat::Sexaquark` +
// `V0` (`Flat::V0`) + `Kaon` (`Flat::Track`) +
// `V0_atPCA` (`Flat::State_NoE`) + `Kaon_atPCA` (`Flat::State_NoE`).
struct alignas(T2S_SIMD_ALIGN) ChannelD : Flat::Sexaquark {
    Flat::V0 V0{};
    Flat::Track Kaon{};
    Flat::State_NoE V0_atPCA{};
    Flat::State_NoE Kaon_atPCA{};

    void CreateBranches_FlatChannelD(TTree* tree, bool is_mc) {
        CreateBranches_FlatSexaquark(tree, is_mc);
        V0.CreateBranches_FlatV0(tree, "V0");
        Kaon.CreateBranches_FlatTrack(tree, "Kaon");
        V0_atPCA.CreateBranches_FlatState_NoE(tree, "V0_atPCA");
        Kaon_atPCA.CreateBranches_FlatState_NoE(tree, "Kaon_atPCA");
    }
};

// `Flat::MC_Sexaquark` +
// `V0` (`Flat::MC_V0`) + `Kaon` (`Flat::MC_Track`).
struct alignas(T2S_SIMD_ALIGN) MC_ChannelD : Flat::MC_Sexaquark {
    Flat::MC_V0 V0{};
    Flat::MC_Track Kaon{};

    void CreateBranches_FlatMC_ChannelD(TTree* tree) {
        CreateBranches_FlatMC_Sexaquark(tree);
        V0.CreateBranches_FlatMC_V0(tree, "V0");
        Kaon.CreateBranches_FlatMC_Track(tree, "Kaon");
    }
};

}  // namespace Tree2Secondaries::Storage::Flat
