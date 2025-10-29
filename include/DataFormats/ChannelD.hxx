#pragma once

#include <TTree.h>

#include "DataFormats/Flat.hxx"
#include "DataFormats/Found.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::DF::Found {

// `Found::Sexaquark` + `V0` (`Flat::V0`) + `Kaon` (`Flat::Track`) + `V0_atPCA` (`Flat::Coordinates `) + `Kaon_atPCA` (`Flat::Coordinates `)
struct alignas(T2S_SIMD_ALIGN) ChannelD : Found::Sexaquark {
    Flat::V0 V0{};
    Flat::Track Kaon{};
    Flat::Coordinates V0_atPCA{};    // MAYBE: could add momentum later
    Flat::Coordinates Kaon_atPCA{};  // MAYBE: could add momentum later

    void CreateBranches_ChannelD(TTree* tree, bool is_mc) {
        CreateBranches_Sexaquark(tree, is_mc);
        V0.CreateBranches_V0(tree, "V0");
        Kaon.CreateBranches_Track(tree, "Kaon");
        V0_atPCA.CreateBranches_Coordinates(tree, "V0_atPCA", "");
        Kaon_atPCA.CreateBranches_Coordinates(tree, "Kaon_atPCA", "");
    }
};

// `Found::MC_Sexaquark` + `V0` (`Flat::MCInfo_V0`) + `Kaon` (`Flat::MCInfo_Bachelor`).
struct alignas(T2S_SIMD_ALIGN) MC_ChannelD : Found::MC_Sexaquark {
    Flat::MCInfo_V0 V0{};
    Flat::MCInfo_Bachelor Kaon{};

    void CreateBranches_MC_ChannelD(TTree* tree) {
        CreateBranches_MC_Sexaquark(tree);
        V0.CreateBranches_MCInfo_V0(tree, "V0");
        Kaon.CreateBranches_MCInfo_Bachelor(tree, "Kaon");
    }
};

}  // namespace Tree2Secondaries::DF::Found
