#pragma once

#include <cstdlib>
#include <format>
#include <string_view>
#include <vector>

#include <TTree.h>

#include "App/Utilities.hxx"
#include "Storage/Vector/BaseVector.hxx"

namespace Tree2Secondaries::Storage::Vector {

// `Vector::States_NoE` + `Vector::CovMatrices_NoE` +
// `Charge` + `DCAxy` + `DCAz` + `TPCSignal` + `NSigmaPion` + `NSigmaKaon` + `NSigmaProton` +
// `McEntry`.
struct Tracks : Vector::States_NoE, Vector::CovMatrices_NoE {
    std::vector<int> *Charge{nullptr};
    std::vector<float> *DCAxy{nullptr};
    std::vector<float> *DCAz{nullptr};
    // pid info
    std::vector<float> *TPCSignal{nullptr};
    std::vector<float> *NSigmaPion{nullptr};
    std::vector<float> *NSigmaKaon{nullptr};
    std::vector<float> *NSigmaProton{nullptr};
    // mc info
    std::vector<int> *McEntry{nullptr};  // NOTE: only read when analyzing MC
    // esd index
    std::vector<size_t> *Index{nullptr};  // NOTE: written by `Packager`, read by `Finder`

    void Clear_VectorTracks(bool include_mc, bool include_esd) {
        Clear_VectorStates_NoE();
        Clear_VectorCovMatrices_NoE();
        Charge->clear();
        DCAxy->clear();
        DCAz->clear();
        TPCSignal->clear();
        NSigmaPion->clear();
        NSigmaKaon->clear();
        NSigmaProton->clear();
        if (include_mc) McEntry->clear();
        if (include_esd) Index->clear();
    }
    void CreateBranches_VectorTracks(TTree *tree, bool include_mc, bool include_esd, std::string_view prefix) {
        CreateBranches_VectorStates_NoE(tree, prefix);
        CreateBranches_VectorCovMatrices_NoE(tree, prefix);
        tree->Branch(std::format("{}_Charge", prefix).c_str(), &Charge);
        tree->Branch(std::format("{}_DCAxy", prefix).c_str(), &DCAxy);
        tree->Branch(std::format("{}_DCAz", prefix).c_str(), &DCAz);
        tree->Branch(std::format("{}_TPCSignal", prefix).c_str(), &TPCSignal);
        tree->Branch(std::format("{}_NSigmaPion", prefix).c_str(), &NSigmaPion);
        tree->Branch(std::format("{}_NSigmaKaon", prefix).c_str(), &NSigmaKaon);
        tree->Branch(std::format("{}_NSigmaProton", prefix).c_str(), &NSigmaProton);
        if (include_mc) tree->Branch(std::format("{}_McEntry", prefix).c_str(), &McEntry);
        if (include_esd) tree->Branch(std::format("{}_Index", prefix).c_str(), &Index);
    }
    void ReadBranches_VectorTracks(TTree *tree, bool include_mc, bool include_esd, std::string_view prefix) {
        ReadBranches_VectorStates_NoE(tree, prefix);
        ReadBranches_VectorCovMatrices_NoE(tree, prefix);
        Utils::ReadBranch(tree, std::format("{}_Charge", prefix), &Charge);
        Utils::ReadBranch(tree, std::format("{}_DCAxy", prefix), &DCAxy);
        Utils::ReadBranch(tree, std::format("{}_DCAz", prefix), &DCAz);
        Utils::ReadBranch(tree, std::format("{}_TPCSignal", prefix), &TPCSignal);
        Utils::ReadBranch(tree, std::format("{}_NSigmaPion", prefix), &NSigmaPion);
        Utils::ReadBranch(tree, std::format("{}_NSigmaKaon", prefix), &NSigmaKaon);
        Utils::ReadBranch(tree, std::format("{}_NSigmaProton", prefix), &NSigmaProton);
        if (include_mc) Utils::ReadBranch(tree, std::format("{}_McEntry", prefix), &McEntry);
        if (include_esd) Utils::ReadBranch(tree, std::format("{}_Index", prefix), &Index);
    }
};

// `Vector::MC` + `Vector::States` +
// `Mother` + `GrandMother`.
struct MC_Tracks : Vector::MC, Vector::States {
    Vector::MC_Id Mother{};
    Vector::MC_Id GrandMother{};

    void Clear_VectorMC_Tracks() {
        Clear_VectorMC();
        Clear_VectorStates();
        Mother.Clear_VectorMC_Id();
        GrandMother.Clear_VectorMC_Id();
    }
    void CreateBranches_VectorMC_Tracks(TTree *tree, std::string_view acronym = "") {
        CreateBranches_VectorMC(tree, acronym);
        CreateBranches_VectorStates(tree, std::format("MC_{}", acronym));
        Mother.CreateBranches_VectorMC_Id(tree, std::format("{}_Mother", acronym));
        GrandMother.CreateBranches_VectorMC_Id(tree, std::format("{}_GrandMother", acronym));
    }
    void ReadBranches_VectorMC_Tracks(TTree *tree, std::string_view acronym = "") {
        ReadBranches_VectorMC(tree, acronym);
        ReadBranches_VectorStates(tree, std::format("MC_{}", acronym));
        Mother.ReadBranches_VectorMC_Id(tree, std::format("{}_Mother", acronym));
        GrandMother.ReadBranches_VectorMC_Id(tree, std::format("{}_GrandMother", acronym));
    }
};

}  // namespace Tree2Secondaries::Storage::Vector
