#pragma once

#include <format>
#include <string_view>
#include <vector>

#include <TTree.h>

#include "App/Utilities.hxx"
#include "DataFormats/StructsOfVectors.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::DF::Packed {

// `SOV::States_NoE` + `SOV::CovMatrices_NoE` + `Index`
struct alignas(T2S_SIMD_ALIGN) Tracks : SOV::States_NoE, SOV::CovMatrices_NoE {
    std::vector<int>* Index{nullptr};

    void Clear_PackedTracks() {
        Clear_States_NoE();
        Clear_CovMatrices_NoE();
        Index->clear();
    }
    void CreateBranches_PackedTracks(TTree* tree, std::string_view prefix) {
        CreateBranches_States_NoE(tree, prefix);
        CreateBranches_CovMatrices_NoE(tree, prefix);
        tree->Branch(std::format("{}_Index", prefix).c_str(), &Index);
    }
    void ReadBranches_PackedTracks(TTree* tree, std::string_view prefix) {
        ReadBranches_States_NoE(tree, prefix);
        ReadBranches_CovMatrices_NoE(tree, prefix);
        Utils::ReadBranch(tree, std::format("{}_Index", prefix), &Index);
    }
};

// `SOV::States` + `SOV::CovMatrices` + `Neg` (`Packed::Tracks`) + `Pos` (`Packed::Tracks`) + `Neg_atPCA` (`States_NoE`) +
// `Pos_atPCA` (`States_NoE`) + `Entry` + `Chi2NDF`
struct alignas(T2S_SIMD_ALIGN) V0s : SOV::States, SOV::CovMatrices {
    Packed::Tracks Neg;
    Packed::Tracks Pos;
    SOV::States_NoE Neg_atPCA;
    SOV::States_NoE Pos_atPCA;
    std::vector<int>* Entry{nullptr};
    std::vector<float>* Chi2NDF{nullptr};

    void Clear_PackedV0s() {
        Clear_States();
        Clear_CovMatrices();
        Neg.Clear_PackedTracks();
        Pos.Clear_PackedTracks();
        Neg_atPCA.Clear_States_NoE();
        Pos_atPCA.Clear_States_NoE();
        Entry->clear();
        Chi2NDF->clear();
    }
    void CreateBranches_PackedV0s(TTree* tree, std::string_view prefix) {
        const std::string& neg_prefix{std::format("{}_Neg", prefix)};
        const std::string& pos_prefix{std::format("{}_Pos", prefix)};
        CreateBranches_States(tree, prefix);
        CreateBranches_CovMatrices(tree, prefix);
        Neg.CreateBranches_PackedTracks(tree, neg_prefix);
        Pos.CreateBranches_PackedTracks(tree, pos_prefix);
        Neg_atPCA.CreateBranches_States_NoE(tree, std::format("{}_atPCA", neg_prefix));
        Pos_atPCA.CreateBranches_States_NoE(tree, std::format("{}_atPCA", pos_prefix));
        tree->Branch(std::format("{}_Entry", prefix).c_str(), &Entry);
        tree->Branch(std::format("{}_Chi2NDF", prefix).c_str(), &Chi2NDF);
    }
    void ReadBranches_PackedV0s(TTree* tree, std::string_view prefix) {
        const std::string& neg_prefix{std::format("{}_Neg", prefix)};
        const std::string& pos_prefix{std::format("{}_Pos", prefix)};
        ReadBranches_States(tree, prefix);
        ReadBranches_CovMatrices(tree, prefix);
        Neg.ReadBranches_PackedTracks(tree, neg_prefix);
        Pos.ReadBranches_PackedTracks(tree, pos_prefix);
        Neg_atPCA.ReadBranches_States_NoE(tree, std::format("{}_atPCA", neg_prefix));
        Pos_atPCA.ReadBranches_States_NoE(tree, std::format("{}_atPCA", pos_prefix));
        Utils::ReadBranch(tree, std::format("{}_Entry", prefix), &Entry);
        Utils::ReadBranch(tree, std::format("{}_Chi2NDF", prefix), &Chi2NDF);
    }
};

// `SOV::MCInfo_States_Mother` + `GrandMother_Entry` + `GrandMother_PdgCode`
struct alignas(T2S_SIMD_ALIGN) LinkedTracks : SOV::MCInfo_States_Mother {
    std::vector<int>* GrandMother_Entry{nullptr};
    std::vector<int>* GrandMother_PdgCode{nullptr};

    void Clear_LinkedTracks() {
        Clear_MCInfo_States_Mother();
        GrandMother_Entry->clear();
        GrandMother_PdgCode->clear();
    }
    void CreateBranches_LinkedTracks(TTree* tree, std::string_view acronym = "") {
        CreateBranches_MCInfo_States_Mother(tree, acronym);
        tree->Branch(std::format("MC_{}_GrandMother_Entry", acronym).c_str(), &GrandMother_Entry);
        tree->Branch(std::format("MC_{}_GrandMother_PdgCode", acronym).c_str(), &GrandMother_PdgCode);
    }
    void ReadBranches_LinkedTracks(TTree* tree, std::string_view acronym = "") {
        ReadBranches_MCInfo_States_Mother(tree, acronym);
        Utils::ReadBranch(tree, std::format("MC_{}_GrandMother_Entry", acronym), &GrandMother_Entry);
        Utils::ReadBranch(tree, std::format("MC_{}_GrandMother_PdgCode", acronym), &GrandMother_PdgCode);
    }
};

// `SOV::MCInfo_States_Mother` + `Neg` (`GrandMother_Entry`) + `Pos` (`GrandMother_PdgCode`) + `DecayX` + `DecayY` + `DecayZ` + `IsHybrid`
struct alignas(T2S_SIMD_ALIGN) LinkedV0s : SOV::MCInfo_States_Mother {
    SOV::MCInfo_PxPyPz Neg;
    SOV::MCInfo_PxPyPz Pos;
    SOV::Coordinates AtDecay;
    std::vector<char>* IsHybrid{nullptr};

    void Clear_LinkedV0s() {
        Clear_MCInfo_States_Mother();
        Neg.Clear_MCInfo_PxPyPz();
        Pos.Clear_MCInfo_PxPyPz();
        AtDecay.Clear_Coordinates();
        IsHybrid->clear();
    }
    void CreateBranches_LinkedV0s(TTree* tree, std::string_view acronym = "") {
        CreateBranches_MCInfo_States_Mother(tree, acronym);
        Neg.CreateBranches_MCInfo_PxPyPz(tree, std::format("{}_Neg", acronym));
        Pos.CreateBranches_MCInfo_PxPyPz(tree, std::format("{}_Pos", acronym));
        AtDecay.CreateBranches_Coordinates(tree, std::format("MC_{}_atDecay", acronym), "");
        tree->Branch(std::format("MC_{}_IsHybrid", acronym).c_str(), &IsHybrid);
    }
    void ReadBranches_LinkedV0s(TTree* tree, std::string_view acronym = "") {
        ReadBranches_MCInfo_States_Mother(tree, acronym);
        Neg.ReadBranches_MCInfo_PxPyPz(tree, std::format("{}_Neg", acronym));
        Pos.ReadBranches_MCInfo_PxPyPz(tree, std::format("{}_Pos", acronym));
        AtDecay.ReadBranches_Coordinates(tree, std::format("MC_{}_atDecay", acronym), "");
        Utils::ReadBranch(tree, std::format("MC_{}_IsHybrid", acronym), &IsHybrid);
    }
};

}  // namespace Tree2Secondaries::DF::Packed
