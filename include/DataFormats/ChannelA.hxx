#pragma once

#include <TTree.h>

#include "DataFormats/Flat.hxx"
#include "DataFormats/Found.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::DF::Found {

// `Found::Sexaquark` + `V0A` (`Flat::V0`) + `V0B` (`Flat::V0`) + `V0A_atPCA` (`Flat::Coordinates`) + `V0B_atPCA` (`Flat::Coordinates`)
struct alignas(T2S_SIMD_ALIGN) ChannelA : Found::Sexaquark {
    Flat::V0 V0A;
    Flat::V0 V0B;
    Flat::Coordinates V0A_atPCA{};  // MAYBE: could add momentum later
    Flat::Coordinates V0B_atPCA{};  // MAYBE: could add momentum later

    void CreateBranches_ChannelA(TTree* tree, bool is_mc) {
        CreateBranches_Sexaquark(tree, is_mc);
        V0A.CreateBranches_V0(tree, "V0A");
        V0B.CreateBranches_V0(tree, "V0B");
        V0A_atPCA.CreateBranches_Coordinates(tree, "V0A_atPCA", "");
        V0B_atPCA.CreateBranches_Coordinates(tree, "V0B_atPCA", "");
    }
};

// `Found::MC_Sexaquark` + `V0A` (`Flat::MCInfo_V0`) + `V0B` (`Flat::MCInfo_V0`).
struct alignas(T2S_SIMD_ALIGN) MC_ChannelA : Found::MC_Sexaquark {
    Flat::MCInfo_V0 V0A{};
    Flat::MCInfo_V0 V0B{};

    void CreateBranches_MC_ChannelA(TTree* tree) {
        CreateBranches_MC_Sexaquark(tree);
        V0A.CreateBranches_MCInfo_V0(tree, "V0A");
        V0B.CreateBranches_MCInfo_V0(tree, "V0B");
    }
};

}  // namespace Tree2Secondaries::DF::Found
