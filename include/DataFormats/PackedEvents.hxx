#ifndef T2S_DF_PACKED_HXX
#define T2S_DF_PACKED_HXX

#include <format>
#include <string_view>
#include <vector>

#include <TTree.h>

#include "App/Utilities.hxx"
#include "DataFormats/DataFormats.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::DF::Packed {

struct alignas(T2S_SIMD_ALIGN) Tracks : SOV::States_NoE, SOV::CovMatrices_NoE {
    std::vector<int>* Entry{nullptr};

    void Clear_PackedTracks() {
        Clear_States_NoE();
        Clear_CovMatrices_NoE();
        Entry->clear();
    }
    void CreateBranches_PackedTracks(TTree* tree, std::string_view prefix = "") {
        CreateBranches_States_NoE(tree, prefix);
        CreateBranches_CovMatrices_NoE(tree, prefix);
        Utils::CreateBranch(tree, std::format("{}_Entry", prefix), &Entry);
    }
    void ReadBranches_PackedTracks(TTree* tree, std::string_view prefix = "") {
        ReadBranches_States_NoE(tree, prefix);
        ReadBranches_CovMatrices_NoE(tree, prefix);
        Utils::ReadBranch(tree, std::format("{}_Entry", prefix), &Entry);
    }
};

struct alignas(T2S_SIMD_ALIGN) V0s : SOV::States, SOV::CovMatrices {
    Tracks Neg;
    Tracks Pos;
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
    void CreateBranches_PackedV0s(TTree* tree, std::string_view prefix = "") {
        const std::string& neg_prefix{std::format("{}_Neg", prefix)};
        const std::string& pos_prefix{std::format("{}_Pos", prefix)};
        CreateBranches_States(tree, prefix);
        CreateBranches_CovMatrices(tree, prefix);
        Neg.CreateBranches_PackedTracks(tree, neg_prefix);
        Pos.CreateBranches_PackedTracks(tree, pos_prefix);
        Neg_atPCA.CreateBranches_States_NoE(tree, std::format("{}_atPCA", neg_prefix));
        Pos_atPCA.CreateBranches_States_NoE(tree, std::format("{}_atPCA", pos_prefix));
        Utils::CreateBranch(tree, std::format("{}_Entry", prefix), &Entry);
        Utils::CreateBranch(tree, std::format("{}_Chi2NDF", prefix), &Chi2NDF);
    }
    void ReadBranches_PackedV0s(TTree* tree, std::string_view prefix = "") {
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

struct alignas(T2S_SIMD_ALIGN) LinkedTracks : SOV::MCInfo {
    std::vector<int>* GrandMother_Entry{nullptr};
    std::vector<int>* GrandMother_PdgCode{nullptr};

    void Clear_LinkedTracks() {
        Clear_MCInfo();
        GrandMother_Entry->clear();
        GrandMother_PdgCode->clear();
    }
    void CreateBranches_LinkedTracks(TTree* tree, std::string_view acronym = "") {
        CreateBranches_MCInfo(tree, acronym);
        Utils::CreateBranch(tree, std::format("MC_{}_GrandMother_Entry", acronym), &GrandMother_Entry);
        Utils::CreateBranch(tree, std::format("MC_{}_GrandMother_PdgCode", acronym), &GrandMother_PdgCode);
    }
    void ReadBranches_LinkedTracks(TTree* tree, std::string_view acronym = "") {
        ReadBranches_MCInfo(tree, acronym);
        Utils::ReadBranch(tree, std::format("MC_{}_GrandMother_Entry", acronym), &GrandMother_Entry);
        Utils::ReadBranch(tree, std::format("MC_{}_GrandMother_PdgCode", acronym), &GrandMother_PdgCode);
    }
};

struct alignas(T2S_SIMD_ALIGN) LinkedV0s : SOV::MCInfo {
    SOV::MCInfo_Reduced Neg;
    SOV::MCInfo_Reduced Pos;
    std::vector<float>* DecayX{nullptr};
    std::vector<float>* DecayY{nullptr};
    std::vector<float>* DecayZ{nullptr};
    std::vector<char>* IsHybrid{nullptr};

    void Clear_LinkedV0s() {
        Clear_MCInfo();
        Neg.Clear_MCInfo_Reduced();
        Pos.Clear_MCInfo_Reduced();
        DecayX->clear();
        DecayY->clear();
        DecayZ->clear();
        IsHybrid->clear();
    }
    void CreateBranches_LinkedV0s(TTree* tree, std::string_view acronym = "") {
        CreateBranches_MCInfo(tree, acronym);
        Neg.CreateBranches_MCInfo_Reduced(tree, std::format("{}_Neg", acronym));
        Pos.CreateBranches_MCInfo_Reduced(tree, std::format("{}_Pos", acronym));
        Utils::CreateBranch(tree, std::format("MC_{}_DecayX", acronym), &DecayX);
        Utils::CreateBranch(tree, std::format("MC_{}_DecayY", acronym), &DecayY);
        Utils::CreateBranch(tree, std::format("MC_{}_DecayZ", acronym), &DecayZ);
        Utils::CreateBranch(tree, std::format("MC_{}_IsHybrid", acronym), &IsHybrid);
    }
    void ReadBranches_LinkedV0s(TTree* tree, std::string_view acronym = "") {
        ReadBranches_MCInfo(tree, acronym);
        Neg.ReadBranches_MCInfo_Reduced(tree, std::format("{}_Neg", acronym));
        Pos.ReadBranches_MCInfo_Reduced(tree, std::format("{}_Pos", acronym));
        Utils::ReadBranch(tree, std::format("MC_{}_DecayX", acronym), &DecayX);
        Utils::ReadBranch(tree, std::format("MC_{}_DecayY", acronym), &DecayY);
        Utils::ReadBranch(tree, std::format("MC_{}_DecayZ", acronym), &DecayZ);
        Utils::ReadBranch(tree, std::format("MC_{}_IsHybrid", acronym), &IsHybrid);
    }
};

}  // namespace Tree2Secondaries::DF::Packed

#endif  // T2S_DF_PACKED_HXX
