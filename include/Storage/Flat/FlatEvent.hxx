#pragma once

#include <TTree.h>

#include "App/Utilities.hxx"
#include "Math/Constants.hxx"
#include "Storage/Flat/BaseFlat.hxx"

namespace Tree2Secondaries::Storage::Flat {

// `PV` (`Flat::Coordinates`) +
// `RunNumber` + `DirNumber` + `DirNumberB` + `EventNumber` +
// `Centrality` + `MagneticField`.
struct alignas(T2S_SIMD_ALIGN) Event {
    Flat::Coordinates PV{};
    unsigned int RunNumber{0};
    unsigned int DirNumber{0};
    unsigned int DirNumberB{0};  // NOTE: only read when analyzing RD
    unsigned int EventNumber{0};
    float Centrality{0.};
    float MagneticField{0.};

    void CreateBranches_FlatEvent(TTree *tree, bool is_mc) {
        PV.CreateBranches_FlatCoordinates(tree, "PV", "v");
        tree->Branch("RunNumber", &RunNumber);
        tree->Branch("DirNumber", &DirNumber);
        if (!is_mc) tree->Branch("DirNumberB", &DirNumberB);
        tree->Branch("EventNumber", &EventNumber);
        tree->Branch("Centrality", &Centrality);
        tree->Branch("MagneticField", &MagneticField);
    }
    void ReadBranches_FlatEvent(TTree *tree, bool is_mc) {
        PV.ReadBranches_FlatCoordinates(tree, "PV", "v");
        Utils::ReadBranch(tree, "RunNumber", &RunNumber);
        Utils::ReadBranch(tree, "DirNumber", &DirNumber);
        if (!is_mc) Utils::ReadBranch(tree, "DirNumberB", &DirNumberB);
        Utils::ReadBranch(tree, "EventNumber", &EventNumber);
        Utils::ReadBranch(tree, "Centrality", &Centrality);
        Utils::ReadBranch(tree, "MagneticField", &MagneticField);
    }
};

}  // namespace Tree2Secondaries::Storage::Flat
