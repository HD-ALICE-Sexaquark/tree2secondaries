#pragma once

#include <vector>

#include <TTree.h>

#include "App/Utilities.hxx"
#include "Storage/Vector/BaseVector.hxx"

namespace Tree2Secondaries::Storage::Vector {

// Read and written by `Packager`.
// `Vector::Coordinates` + `Vector::PxPyPz` + `ReactionID` + `Nucleon` (`Vector::PxPyPz`).
struct Injected : Vector::PxPyPz {
    Vector::Coordinates SV{};
    Vector::PxPyPz Nucleon{};
    std::vector<int> *ReactionID{nullptr};

    void Clear_VectorInjected(bool include_coord) {
        Clear_VectorPxPyPz();
        if (include_coord) SV.Clear_VectorCoordinates();
        Nucleon.Clear_VectorPxPyPz();
        ReactionID->clear();
    }
    void ReadBranches_VectorInjected(TTree *tree, bool include_coord) {
        ReadBranches_VectorPxPyPz(tree, "Sexaquark");
        if (include_coord) SV.ReadBranches_VectorCoordinates(tree, "SV", "");
        Nucleon.ReadBranches_VectorPxPyPz(tree, "Nucleon");
        Utils::ReadBranch(tree, "ReactionID", &ReactionID);
    }
    void CreateBranches_VectorInjected(TTree *tree, bool include_coord) {
        CreateBranches_VectorPxPyPz(tree, "Sexaquark");
        if (include_coord) SV.CreateBranches_VectorCoordinates(tree, "SV", "");
        Nucleon.CreateBranches_VectorPxPyPz(tree, "Nucleon");
        tree->Branch("ReactionID", &ReactionID);
    }
};

}  // namespace Tree2Secondaries::Storage::Vector
