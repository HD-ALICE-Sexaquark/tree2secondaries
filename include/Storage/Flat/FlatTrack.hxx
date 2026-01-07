#pragma once

#include <TTree.h>

#include "Math/Constants.hxx"
#include "Storage/Flat/BaseFlat.hxx"

namespace Tree2Secondaries::Storage::Flat {

// `Flat::State_NoE` +
// `Index` + `DCAxy` + `DCAz` + `TPCSignal` + `NSigmaPion` + `NSigmaKaon` + `NSigmaProton`.
struct alignas(T2S_SIMD_ALIGN) Track : Flat::State_NoE {
    int Index{};
    float DCAxy{};
    float DCAz{};
    float TPCSignal{};
    float NSigmaPion{};
    float NSigmaKaon{};
    float NSigmaProton{};

    void CreateBranches_FlatTrack(TTree* tree, std::string_view prefix) {
        CreateBranches_FlatState_NoE(tree, prefix);
        tree->Branch(std::format("{}_Index", prefix).c_str(), &Index);
        tree->Branch(std::format("{}_DCAxy", prefix).c_str(), &DCAxy);
        tree->Branch(std::format("{}_DCAz", prefix).c_str(), &DCAz);
        tree->Branch(std::format("{}_TPCSignal", prefix).c_str(), &TPCSignal);
        tree->Branch(std::format("{}_NSigmaPion", prefix).c_str(), &NSigmaPion);
        tree->Branch(std::format("{}_NSigmaKaon", prefix).c_str(), &NSigmaKaon);
        tree->Branch(std::format("{}_NSigmaProton", prefix).c_str(), &NSigmaProton);
    }
};

// `Flat::MC` + `Flat::LV` +
// `Mother` (`Flat::MC_Id`) + `GrandMother` (`Flat::MC_Id`).
struct alignas(T2S_SIMD_ALIGN) MC_Track : Flat::MC, Flat::LorentzVector {
    Flat::MC_Id Mother{};
    Flat::MC_Id GrandMother{};

    void CreateBranches_FlatMC_Track(TTree* tree, std::string_view acronym = "") {
        CreateBranches_FlatMC(tree, acronym);
        CreateBranches_FlatLV(tree, std::format("MC_{}", acronym));
        Mother.CreateBranches_FlatMC_Id(tree, acronym);
        GrandMother.CreateBranches_FlatMC_Id(tree, acronym);
    }
};

}  // namespace Tree2Secondaries::Storage::Flat
