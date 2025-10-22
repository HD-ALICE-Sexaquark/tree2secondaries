#ifndef T2S_DF_HXX
#define T2S_DF_HXX

#include <format>
#include <string_view>

#include <TTree.h>

#include "App/Utilities.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::DF {

namespace SOV {
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
    void CreateBranches_CovMatrices_NoE(TTree* tree, std::string_view suffix = "") {
        tree->Branch(std::format("{}_SigmaX2", suffix).c_str(), &SigmaX2);
        tree->Branch(std::format("{}_SigmaXY", suffix).c_str(), &SigmaXY);
        tree->Branch(std::format("{}_SigmaY2", suffix).c_str(), &SigmaY2);
        tree->Branch(std::format("{}_SigmaXZ", suffix).c_str(), &SigmaXZ);
        tree->Branch(std::format("{}_SigmaYZ", suffix).c_str(), &SigmaYZ);
        tree->Branch(std::format("{}_SigmaZ2", suffix).c_str(), &SigmaZ2);
        tree->Branch(std::format("{}_SigmaXPx", suffix).c_str(), &SigmaXPx);
        tree->Branch(std::format("{}_SigmaYPx", suffix).c_str(), &SigmaYPx);
        tree->Branch(std::format("{}_SigmaZPx", suffix).c_str(), &SigmaZPx);
        tree->Branch(std::format("{}_SigmaPx2", suffix).c_str(), &SigmaPx2);
        tree->Branch(std::format("{}_SigmaXPy", suffix).c_str(), &SigmaXPy);
        tree->Branch(std::format("{}_SigmaYPy", suffix).c_str(), &SigmaYPy);
        tree->Branch(std::format("{}_SigmaZPy", suffix).c_str(), &SigmaZPy);
        tree->Branch(std::format("{}_SigmaPxPy", suffix).c_str(), &SigmaPxPy);
        tree->Branch(std::format("{}_SigmaPy2", suffix).c_str(), &SigmaPy2);
        tree->Branch(std::format("{}_SigmaXPz", suffix).c_str(), &SigmaXPz);
        tree->Branch(std::format("{}_SigmaYPz", suffix).c_str(), &SigmaYPz);
        tree->Branch(std::format("{}_SigmaZPz", suffix).c_str(), &SigmaZPz);
        tree->Branch(std::format("{}_SigmaPxPz", suffix).c_str(), &SigmaPxPz);
        tree->Branch(std::format("{}_SigmaPyPz", suffix).c_str(), &SigmaPyPz);
        tree->Branch(std::format("{}_SigmaPz2", suffix).c_str(), &SigmaPz2);
    }
    void ReadBranches_CovMatrices_NoE(TTree* tree, std::string_view suffix = "") {
        Utils::ReadBranch(tree, std::format("{}_SigmaX2", suffix), &SigmaX2);
        Utils::ReadBranch(tree, std::format("{}_SigmaXY", suffix), &SigmaXY);
        Utils::ReadBranch(tree, std::format("{}_SigmaY2", suffix), &SigmaY2);
        Utils::ReadBranch(tree, std::format("{}_SigmaXZ", suffix), &SigmaXZ);
        Utils::ReadBranch(tree, std::format("{}_SigmaYZ", suffix), &SigmaYZ);
        Utils::ReadBranch(tree, std::format("{}_SigmaZ2", suffix), &SigmaZ2);
        Utils::ReadBranch(tree, std::format("{}_SigmaXPx", suffix), &SigmaXPx);
        Utils::ReadBranch(tree, std::format("{}_SigmaYPx", suffix), &SigmaYPx);
        Utils::ReadBranch(tree, std::format("{}_SigmaZPx", suffix), &SigmaZPx);
        Utils::ReadBranch(tree, std::format("{}_SigmaPx2", suffix), &SigmaPx2);
        Utils::ReadBranch(tree, std::format("{}_SigmaXPy", suffix), &SigmaXPy);
        Utils::ReadBranch(tree, std::format("{}_SigmaYPy", suffix), &SigmaYPy);
        Utils::ReadBranch(tree, std::format("{}_SigmaZPy", suffix), &SigmaZPy);
        Utils::ReadBranch(tree, std::format("{}_SigmaPxPy", suffix), &SigmaPxPy);
        Utils::ReadBranch(tree, std::format("{}_SigmaPy2", suffix), &SigmaPy2);
        Utils::ReadBranch(tree, std::format("{}_SigmaXPz", suffix), &SigmaXPz);
        Utils::ReadBranch(tree, std::format("{}_SigmaYPz", suffix), &SigmaYPz);
        Utils::ReadBranch(tree, std::format("{}_SigmaZPz", suffix), &SigmaZPz);
        Utils::ReadBranch(tree, std::format("{}_SigmaPxPz", suffix), &SigmaPxPz);
        Utils::ReadBranch(tree, std::format("{}_SigmaPyPz", suffix), &SigmaPyPz);
        Utils::ReadBranch(tree, std::format("{}_SigmaPz2", suffix), &SigmaPz2);
    }
};

struct alignas(T2S_SIMD_ALIGN) CovMatrices : CovMatrices_NoE {
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
    void CreateBranches_CovMatrices(TTree* tree, std::string_view suffix = "") {
        CreateBranches_CovMatrices_NoE(tree, suffix);
        tree->Branch(std::format("{}_SigmaXE", suffix).c_str(), &SigmaXE);
        tree->Branch(std::format("{}_SigmaYE", suffix).c_str(), &SigmaYE);
        tree->Branch(std::format("{}_SigmaZE", suffix).c_str(), &SigmaZE);
        tree->Branch(std::format("{}_SigmaPxE", suffix).c_str(), &SigmaPxE);
        tree->Branch(std::format("{}_SigmaPyE", suffix).c_str(), &SigmaPyE);
        tree->Branch(std::format("{}_SigmaPzE", suffix).c_str(), &SigmaPzE);
        tree->Branch(std::format("{}_SigmaE2", suffix).c_str(), &SigmaE2);
    }
    void ReadBranches_CovMatrices(TTree* tree, std::string_view suffix = "") {
        ReadBranches_CovMatrices_NoE(tree, suffix);
        Utils::ReadBranch(tree, std::format("{}_SigmaXE", suffix), &SigmaXE);
        Utils::ReadBranch(tree, std::format("{}_SigmaYE", suffix), &SigmaYE);
        Utils::ReadBranch(tree, std::format("{}_SigmaZE", suffix), &SigmaZE);
        Utils::ReadBranch(tree, std::format("{}_SigmaPxE", suffix), &SigmaPxE);
        Utils::ReadBranch(tree, std::format("{}_SigmaPyE", suffix), &SigmaPyE);
        Utils::ReadBranch(tree, std::format("{}_SigmaPzE", suffix), &SigmaPzE);
        Utils::ReadBranch(tree, std::format("{}_SigmaE2", suffix), &SigmaE2);
    }
};

struct alignas(T2S_SIMD_ALIGN) States_NoE {
    std::vector<float>* X{nullptr};
    std::vector<float>* Y{nullptr};
    std::vector<float>* Z{nullptr};
    std::vector<float>* Px{nullptr};
    std::vector<float>* Py{nullptr};
    std::vector<float>* Pz{nullptr};

    void Clear_States_NoE() {
        X->clear();
        Y->clear();
        Z->clear();
        Px->clear();
        Py->clear();
        Pz->clear();
    }
    void CreateBranches_States_NoE(TTree* tree, std::string_view suffix = "") {
        tree->Branch(std::format("{}_X", suffix).c_str(), &X);
        tree->Branch(std::format("{}_Y", suffix).c_str(), &Y);
        tree->Branch(std::format("{}_Z", suffix).c_str(), &Z);
        tree->Branch(std::format("{}_Px", suffix).c_str(), &Px);
        tree->Branch(std::format("{}_Py", suffix).c_str(), &Py);
        tree->Branch(std::format("{}_Pz", suffix).c_str(), &Pz);
    }
    void ReadBranches_States_NoE(TTree* tree, std::string_view suffix = "") {
        Utils::ReadBranch(tree, std::format("{}_X", suffix), &X);
        Utils::ReadBranch(tree, std::format("{}_Y", suffix), &Y);
        Utils::ReadBranch(tree, std::format("{}_Z", suffix), &Z);
        Utils::ReadBranch(tree, std::format("{}_Px", suffix), &Px);
        Utils::ReadBranch(tree, std::format("{}_Py", suffix), &Py);
        Utils::ReadBranch(tree, std::format("{}_Pz", suffix), &Pz);
    }
};

struct alignas(T2S_SIMD_ALIGN) States : States_NoE {
    std::vector<float>* Energy{nullptr};

    void Clear_States() {
        Clear_States_NoE();
        Energy->clear();
    }
    void CreateBranches_States(TTree* tree, std::string_view suffix = "") {
        CreateBranches_States_NoE(tree, suffix);
        tree->Branch(std::format("{}_E", suffix).c_str(), &Energy);
    }
    void ReadBranches_States(TTree* tree, std::string_view suffix = "") {
        ReadBranches_States_NoE(tree, suffix);
        Utils::ReadBranch(tree, std::format("{}_E", suffix), &Energy);
    }
};

struct alignas(T2S_SIMD_ALIGN) MCInfo : States {
    std::vector<int>* Entry{nullptr};
    std::vector<int>* PdgCode{nullptr};
    std::vector<int>* Mother_Entry{nullptr};
    std::vector<int>* Mother_PdgCode{nullptr};
    std::vector<int>* ReactionID{nullptr};
    std::vector<char>* IsTrue{nullptr};
    std::vector<char>* IsSignal{nullptr};
    std::vector<char>* IsSecondary{nullptr};

    void Clear_MCInfo() {
        Clear_States();
        Entry->clear();
        PdgCode->clear();
        Mother_Entry->clear();
        Mother_PdgCode->clear();
        ReactionID->clear();
        IsTrue->clear();
        IsSignal->clear();
        IsSecondary->clear();
    }
    void CreateBranches_MCInfo(TTree* tree, std::string_view acronym = "") {
        CreateBranches_States(tree, std::format("MC_{}", acronym));
        tree->Branch(std::format("MC_{}_Entry", acronym).c_str(), &Entry);
        tree->Branch(std::format("MC_{}_PdgCode", acronym).c_str(), &PdgCode);
        tree->Branch(std::format("MC_{}_Mother_Entry", acronym).c_str(), &Mother_Entry);
        tree->Branch(std::format("MC_{}_Mother_PdgCode", acronym).c_str(), &Mother_PdgCode);
        tree->Branch(std::format("MC_{}_ReactionID", acronym).c_str(), &ReactionID);
        tree->Branch(std::format("MC_{}_IsTrue", acronym).c_str(), &IsTrue);
        tree->Branch(std::format("MC_{}_IsSignal", acronym).c_str(), &IsSignal);
        tree->Branch(std::format("MC_{}_IsSecondary", acronym).c_str(), &IsSecondary);
    }
    void ReadBranches_MCInfo(TTree* tree, std::string_view acronym = "") {
        ReadBranches_States(tree, std::format("MC_{}", acronym));
        Utils::ReadBranch(tree, std::format("MC_{}_Entry", acronym), &Entry);
        Utils::ReadBranch(tree, std::format("MC_{}_PdgCode", acronym), &PdgCode);
        Utils::ReadBranch(tree, std::format("MC_{}_Mother_Entry", acronym), &Mother_Entry);
        Utils::ReadBranch(tree, std::format("MC_{}_Mother_PdgCode", acronym), &Mother_PdgCode);
        Utils::ReadBranch(tree, std::format("MC_{}_ReactionID", acronym), &ReactionID);
        Utils::ReadBranch(tree, std::format("MC_{}_ IsTrue", acronym), &IsTrue);
        Utils::ReadBranch(tree, std::format("MC_{}_ IsSignal", acronym), &IsSignal);
        Utils::ReadBranch(tree, std::format("MC_{}_ IsSecondary", acronym), &IsSecondary);
    }
};

struct alignas(T2S_SIMD_ALIGN) MCInfo_Reduced {
    std::vector<float>* Px{nullptr};
    std::vector<float>* Py{nullptr};
    std::vector<float>* Pz{nullptr};
    std::vector<int>* Entry{nullptr};
    std::vector<int>* PdgCode{nullptr};
    std::vector<int>* ReactionID{nullptr};
    std::vector<char>* IsTrue{nullptr};
    std::vector<char>* IsSignal{nullptr};
    std::vector<char>* IsSecondary{nullptr};

    void Clear_MCInfo_Reduced() {
        Px->clear();
        Py->clear();
        Pz->clear();
        Entry->clear();
        PdgCode->clear();
        ReactionID->clear();
        IsTrue->clear();
        IsSignal->clear();
        IsSecondary->clear();
    }
    void CreateBranches_MCInfo_Reduced(TTree* tree, std::string_view acronym = "") {
        tree->Branch(std::format("MC_{}_Px", acronym).c_str(), &Px);
        tree->Branch(std::format("MC_{}_Py", acronym).c_str(), &Py);
        tree->Branch(std::format("MC_{}_Pz", acronym).c_str(), &Pz);
        tree->Branch(std::format("MC_{}_Entry", acronym).c_str(), &Entry);
        tree->Branch(std::format("MC_{}_PdgCode", acronym).c_str(), &PdgCode);
        tree->Branch(std::format("MC_{}_ReactionID", acronym).c_str(), &ReactionID);
        tree->Branch(std::format("MC_{}_IsTrue", acronym).c_str(), &IsTrue);
        tree->Branch(std::format("MC_{}_IsSignal", acronym).c_str(), &IsSignal);
        tree->Branch(std::format("MC_{}_IsSecondary", acronym).c_str(), &IsSecondary);
    }
    void ReadBranches_MCInfo_Reduced(TTree* tree, std::string_view acronym = "") {
        Utils::ReadBranch(tree, std::format("MC_{}_Px", acronym), &Px);
        Utils::ReadBranch(tree, std::format("MC_{}_Py", acronym), &Py);
        Utils::ReadBranch(tree, std::format("MC_{}_Pz", acronym), &Pz);
        Utils::ReadBranch(tree, std::format("MC_{}_Entry", acronym), &Entry);
        Utils::ReadBranch(tree, std::format("MC_{}_PdgCode", acronym), &PdgCode);
        Utils::ReadBranch(tree, std::format("MC_{}_ReactionID", acronym), &ReactionID);
        Utils::ReadBranch(tree, std::format("MC_{}_ IsTrue", acronym), &IsTrue);
        Utils::ReadBranch(tree, std::format("MC_{}_ IsSignal", acronym), &IsSignal);
        Utils::ReadBranch(tree, std::format("MC_{}_ IsSecondary", acronym), &IsSecondary);
    }
};

}  // namespace SOV

namespace Flat {
struct alignas(T2S_SIMD_ALIGN) Coordinates {
    float X{};
    float Y{};
    float Z{};
    void Fill_Coordinates(float xx, float yy, float zz) {
        X = xx;
        Y = yy;
        Z = zz;
    }
    void CreateBranches_Coordinates(TTree* tree, std::string_view suffix = "") {
        tree->Branch(std::format("{}_X", suffix).c_str(), &X);
        tree->Branch(std::format("{}_Y", suffix).c_str(), &Y);
        tree->Branch(std::format("{}_Z", suffix).c_str(), &Z);
    }
};
struct alignas(T2S_SIMD_ALIGN) LorentzVector {
    float Px{};
    float Py{};
    float Pz{};
    float E{};
    void Fill_LV(float px, float py, float pz, float ee) {
        Px = px;
        Py = py;
        Pz = pz;
        E = ee;
    }
    void CreateBranches_LV(TTree* tree, std::string_view suffix = "") {
        tree->Branch(std::format("{}_Px", suffix).c_str(), &Px);
        tree->Branch(std::format("{}_Py", suffix).c_str(), &Py);
        tree->Branch(std::format("{}_Pz", suffix).c_str(), &Pz);
        tree->Branch(std::format("{}_E", suffix).c_str(), &E);
    }
};
struct alignas(T2S_SIMD_ALIGN) State : Coordinates, LorentzVector {
    void Fill_State(float xx, float yy, float zz, float px, float py, float pz, float ee) {
        Fill_Coordinates(xx, yy, zz);
        Fill_LV(px, py, pz, ee);
    }
    void CreateBranches_State(TTree* tree, std::string_view suffix = "") {
        CreateBranches_Coordinates(tree, suffix);
        CreateBranches_LV(tree, suffix);
    }
};
}  // namespace Flat

}  // namespace Tree2Secondaries::DF

#endif  // T2S_DF_HXX
