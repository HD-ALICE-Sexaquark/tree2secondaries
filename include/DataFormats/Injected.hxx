#pragma once

#include <vector>

#include <TTree.h>

#include "App/Utilities.hxx"
#include "DataFormats/Flat.hxx"
#include "DataFormats/StructsOfVectors.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::DF {

// Read and written by `Packager` //

namespace SOV {

// `SOV::Coordinates` + `SOV::PxPyPz` + `ReactionID` + `Nucleon` (`SOV::PxPyPz`).
struct alignas(T2S_SIMD_ALIGN) Injected : SOV::Coordinates, SOV::PxPyPz {
    SOV::PxPyPz Nucleon;
    std::vector<int> *ReactionID{nullptr};

    void Clear_SOV_Injected(bool include_coord) {
        if (include_coord) Clear_Coordinates();
        Clear_PxPyPz();
        Nucleon.Clear_PxPyPz();
        ReactionID->clear();
    }
    void ReadBranches_SOV_Injected(TTree *tree, bool include_coord) {
        if (include_coord) ReadBranches_Coordinates(tree, "SV", "");
        ReadBranches_PxPyPz(tree, "Sexaquark");
        Nucleon.ReadBranches_PxPyPz(tree, "Nucleon");
        Utils::ReadBranch(tree, "ReactionID", &ReactionID);
    }
    void CreateBranches_SOV_Injected(TTree *tree, bool include_coord) {
        if (include_coord) CreateBranches_Coordinates(tree, "SV", "");
        CreateBranches_PxPyPz(tree, "Sexaquark");
        Nucleon.CreateBranches_PxPyPz(tree, "Nucleon");
        tree->Branch("ReactionID", &ReactionID);
    }
};
}  // namespace SOV

// Written by `Finder` //

namespace Flat {

// `Flat::Coordinates` + `Flat::LorentzVector` + `Nucleon` (`Flat::LorentzVector`) + `RunNumber` + `DirNumber` + `EventNumber` + `ReactionID`.
struct alignas(T2S_SIMD_ALIGN) Injected : Flat::Coordinates, Flat::LorentzVector {
    Flat::LorentzVector Nucleon;
    unsigned int RunNumber{};
    unsigned int DirNumber{};
    unsigned int EventNumber{};
    int ReactionID{};

    void CreateBranches_Flat_Injected(TTree *tree) {
        CreateBranches_Coordinates(tree, "SV", "");
        CreateBranches_LorentzVector(tree, "Sexaquark");
        Nucleon.CreateBranches_LorentzVector(tree, "Nucleon");
        tree->Branch("RunNumber", &RunNumber);
        tree->Branch("DirNumber", &DirNumber);
        tree->Branch("EventNumber", &EventNumber);
        tree->Branch("ReactionID", &ReactionID);
    }
};
}  // namespace Flat

}  // namespace Tree2Secondaries::DF
