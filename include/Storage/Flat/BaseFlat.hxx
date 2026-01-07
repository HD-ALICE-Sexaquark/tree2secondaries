#pragma once

#include <format>
#include <string_view>

#include <TTree.h>

#include "App/Utilities.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::Storage::Flat {

struct alignas(T2S_SIMD_ALIGN) Coordinates {
    float X{};
    float Y{};
    float Z{};

    void CreateBranches_FlatCoordinates(TTree* tree, std::string_view prefix, std::string_view suffix) {
        tree->Branch(std::format("{}_X{}", prefix, suffix).c_str(), &X);
        tree->Branch(std::format("{}_Y{}", prefix, suffix).c_str(), &Y);
        tree->Branch(std::format("{}_Z{}", prefix, suffix).c_str(), &Z);
    }
    void ReadBranches_FlatCoordinates(TTree* tree, std::string_view prefix, std::string_view suffix) {
        Utils::ReadBranch(tree, std::format("{}_X{}", prefix, suffix), &X);
        Utils::ReadBranch(tree, std::format("{}_Y{}", prefix, suffix), &Y);
        Utils::ReadBranch(tree, std::format("{}_Z{}", prefix, suffix), &Z);
    }
};

struct alignas(T2S_SIMD_ALIGN) PxPyPz {
    float Px{};
    float Py{};
    float Pz{};

    void CreateBranches_FlatPxPyPz(TTree* tree, std::string_view prefix) {
        tree->Branch(std::format("{}_Px", prefix).c_str(), &Px);
        tree->Branch(std::format("{}_Py", prefix).c_str(), &Py);
        tree->Branch(std::format("{}_Pz", prefix).c_str(), &Pz);
    }
};

struct alignas(T2S_SIMD_ALIGN) LorentzVector : Flat::PxPyPz {
    float Energy{};

    void CreateBranches_FlatLV(TTree* tree, std::string_view prefix) {
        CreateBranches_FlatPxPyPz(tree, prefix);
        tree->Branch(std::format("{}_E", prefix).c_str(), &Energy);
    }
};

struct alignas(T2S_SIMD_ALIGN) State : Flat::Coordinates, Flat::LorentzVector {
    void CreateBranches_FlatState(TTree* tree, std::string_view prefix) {
        CreateBranches_FlatCoordinates(tree, prefix, "");
        CreateBranches_FlatLV(tree, prefix);
    }
};

struct alignas(T2S_SIMD_ALIGN) State_NoE : Flat::Coordinates, Flat::PxPyPz {
    void CreateBranches_FlatState_NoE(TTree* tree, std::string_view prefix) {
        CreateBranches_FlatCoordinates(tree, prefix, "");
        CreateBranches_FlatPxPyPz(tree, prefix);
    }
};

// MC Information //

// `Entry` + `PdgCode`.
struct alignas(T2S_SIMD_ALIGN) MC_Id {
    int McEntry{};
    int PdgCode{};

    // NOTE: will add `MC_` as prefix.
    void CreateBranches_FlatMC_Id(TTree* tree, std::string_view acronym = "") {
        tree->Branch(std::format("MC_{}_McEntry", acronym).c_str(), &McEntry);
        tree->Branch(std::format("MC_{}_PdgCode", acronym).c_str(), &PdgCode);
    }
};

// `ReactionID` + `IsTrue` + `IsSignal` + `IsSecondary`.
struct alignas(T2S_SIMD_ALIGN) MC_Flags {
    int ReactionID{};
    bool IsTrue{};
    bool IsSignal{};
    bool IsSecondary{};

    // NOTE: will add `MC_` as prefix.
    void CreateBranches_FlatMC_Flags(TTree* tree, std::string_view acronym = "") {
        tree->Branch(std::format("MC_{}_ReactionID", acronym).c_str(), &ReactionID);
        tree->Branch(std::format("MC_{}_IsTrue", acronym).c_str(), &IsTrue);
        tree->Branch(std::format("MC_{}_IsSignal", acronym).c_str(), &IsSignal);
        tree->Branch(std::format("MC_{}_IsSecondary", acronym).c_str(), &IsSecondary);
    }
};

// `Flat::MC_Id` + `Flat::MC_Flags`.
struct alignas(T2S_SIMD_ALIGN) MC : Flat::MC_Id, Flat::MC_Flags {
    void CreateBranches_FlatMC(TTree* tree, std::string_view acronym = "") {
        CreateBranches_FlatMC_Id(tree, acronym);
        CreateBranches_FlatMC_Flags(tree, acronym);
    }
};

}  // namespace Tree2Secondaries::Storage::Flat
