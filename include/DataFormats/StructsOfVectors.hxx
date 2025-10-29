#pragma once

#include <format>
#include <string_view>

#include <TTree.h>

#include "App/Utilities.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::DF::SOV {

// General //

struct alignas(T2S_SIMD_ALIGN) Coordinates {
    std::vector<float>* X{nullptr};
    std::vector<float>* Y{nullptr};
    std::vector<float>* Z{nullptr};

    void Clear_Coordinates() {
        X->clear();
        Y->clear();
        Z->clear();
    }
    void CreateBranches_Coordinates(TTree* tree, std::string_view prefix, std::string_view suffix) {
        tree->Branch(std::format("{}_X{}", prefix, suffix).c_str(), &X);
        tree->Branch(std::format("{}_Y{}", prefix, suffix).c_str(), &Y);
        tree->Branch(std::format("{}_Z{}", prefix, suffix).c_str(), &Z);
    }
    void ReadBranches_Coordinates(TTree* tree, std::string_view prefix, std::string_view suffix) {
        Utils::ReadBranch(tree, std::format("{}_X{}", prefix, suffix), &X);
        Utils::ReadBranch(tree, std::format("{}_Y{}", prefix, suffix), &Y);
        Utils::ReadBranch(tree, std::format("{}_Z{}", prefix, suffix), &Z);
    }
};

struct alignas(T2S_SIMD_ALIGN) PxPyPz {
    std::vector<float>* Px{nullptr};
    std::vector<float>* Py{nullptr};
    std::vector<float>* Pz{nullptr};

    void Clear_PxPyPz() {
        Px->clear();
        Py->clear();
        Pz->clear();
    }
    void CreateBranches_PxPyPz(TTree* tree, std::string_view prefix) {
        tree->Branch(std::format("{}_Px", prefix).c_str(), &Px);
        tree->Branch(std::format("{}_Py", prefix).c_str(), &Py);
        tree->Branch(std::format("{}_Pz", prefix).c_str(), &Pz);
    }
    void ReadBranches_PxPyPz(TTree* tree, std::string_view prefix) {
        Utils::ReadBranch(tree, std::format("{}_Px", prefix), &Px);
        Utils::ReadBranch(tree, std::format("{}_Py", prefix), &Py);
        Utils::ReadBranch(tree, std::format("{}_Pz", prefix), &Pz);
    }
};

struct alignas(T2S_SIMD_ALIGN) States_NoE : SOV::Coordinates, SOV::PxPyPz {
    void Clear_States_NoE() {
        Clear_Coordinates();
        Clear_PxPyPz();
    }
    void CreateBranches_States_NoE(TTree* tree, std::string_view prefix) {
        CreateBranches_Coordinates(tree, prefix, "");
        CreateBranches_PxPyPz(tree, prefix);
    }
    void ReadBranches_States_NoE(TTree* tree, std::string_view prefix) {
        ReadBranches_Coordinates(tree, prefix, "");
        ReadBranches_PxPyPz(tree, prefix);
    }
};

struct alignas(T2S_SIMD_ALIGN) States : SOV::States_NoE {
    std::vector<float>* Energy{nullptr};

    void Clear_States() {
        Clear_States_NoE();
        Energy->clear();
    }
    void CreateBranches_States(TTree* tree, std::string_view prefix) {
        CreateBranches_States_NoE(tree, prefix);
        tree->Branch(std::format("{}_E", prefix).c_str(), &Energy);
    }
    void ReadBranches_States(TTree* tree, std::string_view prefix) {
        ReadBranches_States_NoE(tree, prefix);
        Utils::ReadBranch(tree, std::format("{}_E", prefix), &Energy);
    }
};

struct alignas(T2S_SIMD_ALIGN) CovMatrices_NoE {
    std::vector<float>* SigmaX2{nullptr};
    std::vector<float>* SigmaXY{nullptr};
    std::vector<float>* SigmaY2{nullptr};
    std::vector<float>* SigmaXZ{nullptr};
    std::vector<float>* SigmaYZ{nullptr};
    std::vector<float>* SigmaZ2{nullptr};
    std::vector<float>* SigmaXPx{nullptr};
    std::vector<float>* SigmaYPx{nullptr};
    std::vector<float>* SigmaZPx{nullptr};
    std::vector<float>* SigmaPx2{nullptr};
    std::vector<float>* SigmaXPy{nullptr};
    std::vector<float>* SigmaYPy{nullptr};
    std::vector<float>* SigmaZPy{nullptr};
    std::vector<float>* SigmaPxPy{nullptr};
    std::vector<float>* SigmaPy2{nullptr};
    std::vector<float>* SigmaXPz{nullptr};
    std::vector<float>* SigmaYPz{nullptr};
    std::vector<float>* SigmaZPz{nullptr};
    std::vector<float>* SigmaPxPz{nullptr};
    std::vector<float>* SigmaPyPz{nullptr};
    std::vector<float>* SigmaPz2{nullptr};

    void Clear_CovMatrices_NoE() {
        SigmaX2->clear();
        SigmaXY->clear();
        SigmaY2->clear();
        SigmaXZ->clear();
        SigmaYZ->clear();
        SigmaZ2->clear();
        SigmaXPx->clear();
        SigmaYPx->clear();
        SigmaZPx->clear();
        SigmaPx2->clear();
        SigmaXPy->clear();
        SigmaYPy->clear();
        SigmaZPy->clear();
        SigmaPxPy->clear();
        SigmaPy2->clear();
        SigmaXPz->clear();
        SigmaYPz->clear();
        SigmaZPz->clear();
        SigmaPxPz->clear();
        SigmaPyPz->clear();
        SigmaPz2->clear();
    }
    void CreateBranches_CovMatrices_NoE(TTree* tree, std::string_view prefix) {
        tree->Branch(std::format("{}_SigmaX2", prefix).c_str(), &SigmaX2);
        tree->Branch(std::format("{}_SigmaXY", prefix).c_str(), &SigmaXY);
        tree->Branch(std::format("{}_SigmaY2", prefix).c_str(), &SigmaY2);
        tree->Branch(std::format("{}_SigmaXZ", prefix).c_str(), &SigmaXZ);
        tree->Branch(std::format("{}_SigmaYZ", prefix).c_str(), &SigmaYZ);
        tree->Branch(std::format("{}_SigmaZ2", prefix).c_str(), &SigmaZ2);
        tree->Branch(std::format("{}_SigmaXPx", prefix).c_str(), &SigmaXPx);
        tree->Branch(std::format("{}_SigmaYPx", prefix).c_str(), &SigmaYPx);
        tree->Branch(std::format("{}_SigmaZPx", prefix).c_str(), &SigmaZPx);
        tree->Branch(std::format("{}_SigmaPx2", prefix).c_str(), &SigmaPx2);
        tree->Branch(std::format("{}_SigmaXPy", prefix).c_str(), &SigmaXPy);
        tree->Branch(std::format("{}_SigmaYPy", prefix).c_str(), &SigmaYPy);
        tree->Branch(std::format("{}_SigmaZPy", prefix).c_str(), &SigmaZPy);
        tree->Branch(std::format("{}_SigmaPxPy", prefix).c_str(), &SigmaPxPy);
        tree->Branch(std::format("{}_SigmaPy2", prefix).c_str(), &SigmaPy2);
        tree->Branch(std::format("{}_SigmaXPz", prefix).c_str(), &SigmaXPz);
        tree->Branch(std::format("{}_SigmaYPz", prefix).c_str(), &SigmaYPz);
        tree->Branch(std::format("{}_SigmaZPz", prefix).c_str(), &SigmaZPz);
        tree->Branch(std::format("{}_SigmaPxPz", prefix).c_str(), &SigmaPxPz);
        tree->Branch(std::format("{}_SigmaPyPz", prefix).c_str(), &SigmaPyPz);
        tree->Branch(std::format("{}_SigmaPz2", prefix).c_str(), &SigmaPz2);
    }
    void ReadBranches_CovMatrices_NoE(TTree* tree, std::string_view prefix) {
        Utils::ReadBranch(tree, std::format("{}_SigmaX2", prefix), &SigmaX2);
        Utils::ReadBranch(tree, std::format("{}_SigmaXY", prefix), &SigmaXY);
        Utils::ReadBranch(tree, std::format("{}_SigmaY2", prefix), &SigmaY2);
        Utils::ReadBranch(tree, std::format("{}_SigmaXZ", prefix), &SigmaXZ);
        Utils::ReadBranch(tree, std::format("{}_SigmaYZ", prefix), &SigmaYZ);
        Utils::ReadBranch(tree, std::format("{}_SigmaZ2", prefix), &SigmaZ2);
        Utils::ReadBranch(tree, std::format("{}_SigmaXPx", prefix), &SigmaXPx);
        Utils::ReadBranch(tree, std::format("{}_SigmaYPx", prefix), &SigmaYPx);
        Utils::ReadBranch(tree, std::format("{}_SigmaZPx", prefix), &SigmaZPx);
        Utils::ReadBranch(tree, std::format("{}_SigmaPx2", prefix), &SigmaPx2);
        Utils::ReadBranch(tree, std::format("{}_SigmaXPy", prefix), &SigmaXPy);
        Utils::ReadBranch(tree, std::format("{}_SigmaYPy", prefix), &SigmaYPy);
        Utils::ReadBranch(tree, std::format("{}_SigmaZPy", prefix), &SigmaZPy);
        Utils::ReadBranch(tree, std::format("{}_SigmaPxPy", prefix), &SigmaPxPy);
        Utils::ReadBranch(tree, std::format("{}_SigmaPy2", prefix), &SigmaPy2);
        Utils::ReadBranch(tree, std::format("{}_SigmaXPz", prefix), &SigmaXPz);
        Utils::ReadBranch(tree, std::format("{}_SigmaYPz", prefix), &SigmaYPz);
        Utils::ReadBranch(tree, std::format("{}_SigmaZPz", prefix), &SigmaZPz);
        Utils::ReadBranch(tree, std::format("{}_SigmaPxPz", prefix), &SigmaPxPz);
        Utils::ReadBranch(tree, std::format("{}_SigmaPyPz", prefix), &SigmaPyPz);
        Utils::ReadBranch(tree, std::format("{}_SigmaPz2", prefix), &SigmaPz2);
    }
};

struct alignas(T2S_SIMD_ALIGN) CovMatrices : SOV::CovMatrices_NoE {
    std::vector<float>* SigmaXE{nullptr};
    std::vector<float>* SigmaYE{nullptr};
    std::vector<float>* SigmaZE{nullptr};
    std::vector<float>* SigmaPxE{nullptr};
    std::vector<float>* SigmaPyE{nullptr};
    std::vector<float>* SigmaPzE{nullptr};
    std::vector<float>* SigmaE2{nullptr};

    void Clear_CovMatrices() {
        Clear_CovMatrices_NoE();
        SigmaXE->clear();
        SigmaYE->clear();
        SigmaZE->clear();
        SigmaPxE->clear();
        SigmaPyE->clear();
        SigmaPzE->clear();
        SigmaE2->clear();
    }
    void CreateBranches_CovMatrices(TTree* tree, std::string_view prefix) {
        CreateBranches_CovMatrices_NoE(tree, prefix);
        tree->Branch(std::format("{}_SigmaXE", prefix).c_str(), &SigmaXE);
        tree->Branch(std::format("{}_SigmaYE", prefix).c_str(), &SigmaYE);
        tree->Branch(std::format("{}_SigmaZE", prefix).c_str(), &SigmaZE);
        tree->Branch(std::format("{}_SigmaPxE", prefix).c_str(), &SigmaPxE);
        tree->Branch(std::format("{}_SigmaPyE", prefix).c_str(), &SigmaPyE);
        tree->Branch(std::format("{}_SigmaPzE", prefix).c_str(), &SigmaPzE);
        tree->Branch(std::format("{}_SigmaE2", prefix).c_str(), &SigmaE2);
    }
    void ReadBranches_CovMatrices(TTree* tree, std::string_view prefix) {
        ReadBranches_CovMatrices_NoE(tree, prefix);
        Utils::ReadBranch(tree, std::format("{}_SigmaXE", prefix), &SigmaXE);
        Utils::ReadBranch(tree, std::format("{}_SigmaYE", prefix), &SigmaYE);
        Utils::ReadBranch(tree, std::format("{}_SigmaZE", prefix), &SigmaZE);
        Utils::ReadBranch(tree, std::format("{}_SigmaPxE", prefix), &SigmaPxE);
        Utils::ReadBranch(tree, std::format("{}_SigmaPyE", prefix), &SigmaPyE);
        Utils::ReadBranch(tree, std::format("{}_SigmaPzE", prefix), &SigmaPzE);
        Utils::ReadBranch(tree, std::format("{}_SigmaE2", prefix), &SigmaE2);
    }
};

// MC Information //

// `Entry` + `PdgCode` + `ReactionID` + `IsTrue` + `IsSignal` + `IsSecondary`.
struct alignas(T2S_SIMD_ALIGN) MCInfo {

    std::vector<int>* Entry{nullptr};
    std::vector<int>* PdgCode{nullptr};
    std::vector<int>* ReactionID{nullptr};
    std::vector<char>* IsTrue{nullptr};
    std::vector<char>* IsSignal{nullptr};
    std::vector<char>* IsSecondary{nullptr};

    void Clear_MCInfo() {
        Entry->clear();
        PdgCode->clear();
        ReactionID->clear();
        IsTrue->clear();
        IsSignal->clear();
        IsSecondary->clear();
    }
    void CreateBranches_MCInfo(TTree* tree, std::string_view acronym = "") {
        tree->Branch(std::format("MC_{}_Entry", acronym).c_str(), &Entry);
        tree->Branch(std::format("MC_{}_PdgCode", acronym).c_str(), &PdgCode);
        tree->Branch(std::format("MC_{}_ReactionID", acronym).c_str(), &ReactionID);
        tree->Branch(std::format("MC_{}_IsTrue", acronym).c_str(), &IsTrue);
        tree->Branch(std::format("MC_{}_IsSignal", acronym).c_str(), &IsSignal);
        tree->Branch(std::format("MC_{}_IsSecondary", acronym).c_str(), &IsSecondary);
    }
    void ReadBranches_MCInfo(TTree* tree, std::string_view acronym = "") {
        Utils::ReadBranch(tree, std::format("MC_{}_Entry", acronym), &Entry);
        Utils::ReadBranch(tree, std::format("MC_{}_PdgCode", acronym), &PdgCode);
        Utils::ReadBranch(tree, std::format("MC_{}_ReactionID", acronym), &ReactionID);
        Utils::ReadBranch(tree, std::format("MC_{}_IsTrue", acronym), &IsTrue);
        Utils::ReadBranch(tree, std::format("MC_{}_IsSignal", acronym), &IsSignal);
        Utils::ReadBranch(tree, std::format("MC_{}_IsSecondary", acronym), &IsSecondary);
    }
};

// `SOV::MCInfo` + `SOV::PxPyPz`.
struct alignas(T2S_SIMD_ALIGN) MCInfo_PxPyPz : SOV::MCInfo, SOV::PxPyPz {
    void Clear_MCInfo_PxPyPz() {
        Clear_MCInfo();
        Clear_PxPyPz();
    }
    void CreateBranches_MCInfo_PxPyPz(TTree* tree, std::string_view acronym = "") {
        CreateBranches_MCInfo(tree, acronym);
        CreateBranches_PxPyPz(tree, std::format("MC_{}", acronym));
    }
    void ReadBranches_MCInfo_PxPyPz(TTree* tree, std::string_view acronym = "") {
        ReadBranches_MCInfo(tree, acronym);
        ReadBranches_PxPyPz(tree, std::format("MC_{}", acronym));
    }
};

// `SOV::MCInfo` + `SOV::States` + `Mother_Entry` + `Mother_PdgCode`.
struct alignas(T2S_SIMD_ALIGN) MCInfo_States_Mother : SOV::MCInfo, SOV::States {

    std::vector<int>* Mother_Entry{nullptr};
    std::vector<int>* Mother_PdgCode{nullptr};

    void Clear_MCInfo_States_Mother() {
        Clear_MCInfo();
        Clear_States();
        Mother_Entry->clear();
        Mother_PdgCode->clear();
    }
    void CreateBranches_MCInfo_States_Mother(TTree* tree, std::string_view acronym = "") {
        CreateBranches_MCInfo(tree, std::format("{}", acronym));
        CreateBranches_States(tree, std::format("MC_{}", acronym));
        tree->Branch(std::format("MC_{}_Mother_Entry", acronym).c_str(), &Mother_Entry);
        tree->Branch(std::format("MC_{}_Mother_PdgCode", acronym).c_str(), &Mother_PdgCode);
    }
    void ReadBranches_MCInfo_States_Mother(TTree* tree, std::string_view acronym = "") {
        ReadBranches_MCInfo(tree, std::format("{}", acronym));
        ReadBranches_States(tree, std::format("MC_{}", acronym));
        Utils::ReadBranch(tree, std::format("MC_{}_Mother_Entry", acronym), &Mother_Entry);
        Utils::ReadBranch(tree, std::format("MC_{}_Mother_PdgCode", acronym), &Mother_PdgCode);
    }
};

}  // namespace Tree2Secondaries::DF::SOV
