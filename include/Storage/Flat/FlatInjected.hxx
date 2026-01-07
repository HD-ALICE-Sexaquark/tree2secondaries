#pragma once

#include <TTree.h>

#include "Math/Constants.hxx"
#include "Storage/Flat/BaseFlat.hxx"

namespace Tree2Secondaries::Storage::Flat {

// Written by `Finder`.
// Flattened version of `Vector::Injected`, +5 new properties.
// `Flat::LV`
// `SV` (`Flat::Coordinates`) + `Nucleon` (`Flat::LV`) +
// `RunNumber` + `DirNumber` + `EventNumber` + `ReactionID`.
struct alignas(T2S_SIMD_ALIGN) Injected : Flat::LorentzVector {
    Flat::LorentzVector Nucleon{};
    Flat::Coordinates SV{};
    unsigned int RunNumber{};
    unsigned int DirNumber{};
    unsigned int EventNumber{};
    int ReactionID{};

    void CreateBranches_FlatInjected(TTree *tree) {
        CreateBranches_FlatLV(tree, "Sexa");
        Nucleon.CreateBranches_FlatLV(tree, "Nucleon");
        SV.CreateBranches_FlatCoordinates(tree, "SV", "");
        tree->Branch("RunNumber", &RunNumber);
        tree->Branch("DirNumber", &DirNumber);
        tree->Branch("EventNumber", &EventNumber);
        tree->Branch("ReactionID", &ReactionID);
    }
};

}  // namespace Tree2Secondaries::Storage::Flat
