#pragma once

#include <format>
#include <string_view>

#include <TTree.h>

#include "App/Utilities.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::DF::Flat {

// General //

struct alignas(T2S_SIMD_ALIGN) Coordinates {
    float X{};
    float Y{};
    float Z{};

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
    void Fill_Coordinates(float xx, float yy, float zz) {
        X = xx;
        Y = yy;
        Z = zz;
    }
};

struct alignas(T2S_SIMD_ALIGN) PxPyPz {
    float Px{};
    float Py{};
    float Pz{};

    void CreateBranches_PxPyPz(TTree* tree, std::string_view prefix) {
        tree->Branch(std::format("{}_Px", prefix).c_str(), &Px);
        tree->Branch(std::format("{}_Py", prefix).c_str(), &Py);
        tree->Branch(std::format("{}_Pz", prefix).c_str(), &Pz);
    }
    void Fill_PxPyPz(float px, float py, float pz) {
        Px = px;
        Py = py;
        Pz = pz;
    }
};

struct alignas(T2S_SIMD_ALIGN) LorentzVector : Flat::PxPyPz {
    float E{};

    void CreateBranches_LorentzVector(TTree* tree, std::string_view prefix) {
        CreateBranches_PxPyPz(tree, prefix);
        tree->Branch(std::format("{}_E", prefix).c_str(), &E);
    }
    void Fill_LorentzVector(float px, float py, float pz, float ee) {
        Fill_PxPyPz(px, py, pz);
        E = ee;
    }
};

struct alignas(T2S_SIMD_ALIGN) State : Flat::Coordinates, Flat::LorentzVector {
    void CreateBranches_State(TTree* tree, std::string_view prefix) {
        CreateBranches_Coordinates(tree, prefix, "");
        CreateBranches_LorentzVector(tree, prefix);
    }
    void Fill_State(float xx, float yy, float zz, float px, float py, float pz, float ee) {
        Fill_Coordinates(xx, yy, zz);
        Fill_LorentzVector(px, py, pz, ee);
    }
};

// `Flat::State` + `Index`.
struct alignas(T2S_SIMD_ALIGN) Track : Flat::State {
    int Index{};

    void CreateBranches_Track(TTree* tree, std::string_view prefix) {
        CreateBranches_State(tree, prefix);
        tree->Branch(std::format("{}_Index", prefix).c_str(), &Index);
    }
};

// `Flat::State` + `Neg_atV0` (`Flat::Track`) + `Pos_atV0` (`Flat::Track`) + `Entry`.
struct alignas(T2S_SIMD_ALIGN) V0 : Flat::State {
    Flat::Track Neg_atV0{};
    Flat::Track Pos_atV0{};
    int Entry{};

    void CreateBranches_V0(TTree* tree, std::string_view prefix) {
        CreateBranches_State(tree, prefix);
        Neg_atV0.CreateBranches_Track(tree, std::format("{}_Neg", prefix));
        Pos_atV0.CreateBranches_Track(tree, std::format("{}_Pos", prefix));
        tree->Branch(std::format("{}_Entry", prefix).c_str(), &Entry);
    }
};

// MC Information //

// `Entry` + `PdgCode` + `ReactionID` + `IsTrue` + `IsSignal` + `IsSecondary`.
// NOTE: member functions will add the `MC_` prefix.
struct alignas(T2S_SIMD_ALIGN) MCInfo {
    int Entry{};
    int PdgCode{};
    int ReactionID{};
    bool IsTrue{};
    bool IsSignal{};
    bool IsSecondary{};

    void CreateBranches_MCInfo(TTree* tree, std::string_view acronym = "") {
        tree->Branch(std::format("MC_{}_Entry", acronym).c_str(), &Entry);
        tree->Branch(std::format("MC_{}_PdgCode", acronym).c_str(), &PdgCode);
        tree->Branch(std::format("MC_{}_ReactionID", acronym).c_str(), &ReactionID);
        tree->Branch(std::format("MC_{}_IsTrue", acronym).c_str(), &IsTrue);
        tree->Branch(std::format("MC_{}_IsSignal", acronym).c_str(), &IsSignal);
        tree->Branch(std::format("MC_{}_IsSecondary", acronym).c_str(), &IsSecondary);
    }
};

// `Flat::MCInfo` + `Flat::PxPyPz`.
// NOTE: member functions will add the `MC_` prefix.
struct alignas(T2S_SIMD_ALIGN) MCInfo_PxPyPz : Flat::MCInfo, Flat::PxPyPz {
    void CreateBranches_MCInfo_PxPyPz(TTree* tree, std::string_view acronym = "") {
        CreateBranches_MCInfo(tree, acronym);
        CreateBranches_PxPyPz(tree, std::format("MC_{}", acronym));
    }
};

// `Flat::MCInfo` + `Flat::LorentzVector` + `Mother_Entry` + `Mother_PdgCode`.
// NOTE: member functions will add the `MC_` prefix.
struct alignas(T2S_SIMD_ALIGN) MCInfo_LV_Mother : Flat::MCInfo, Flat::LorentzVector {
    int Mother_Entry{};
    int Mother_PdgCode{};

    void CreateBranches_MCInfo_LV_Mother(TTree* tree, std::string_view acronym = "") {
        CreateBranches_MCInfo(tree, acronym);
        CreateBranches_LorentzVector(tree, std::format("MC_{}", acronym));
        tree->Branch(std::format("MC_{}_Mother_Entry", acronym).c_str(), &Mother_Entry);
        tree->Branch(std::format("MC_{}_Mother_PdgCode", acronym).c_str(), &Mother_PdgCode);
    }
};

// `Flat::MCInfo_LV_Mother` + `GrandMother_Entry` + `GrandMother_PdgCode`.
// NOTE: member functions will add the `MC_` prefix.
struct alignas(T2S_SIMD_ALIGN) MCInfo_Bachelor : Flat::MCInfo_LV_Mother {
    int GrandMother_Entry{};
    int GrandMother_PdgCode{};

    void CreateBranches_MCInfo_Bachelor(TTree* tree, std::string_view acronym = "") {  //
        CreateBranches_MCInfo_LV_Mother(tree, acronym);
        tree->Branch(std::format("MC_{}_GrandMother_Entry", acronym).c_str(), &GrandMother_Entry);
        tree->Branch(std::format("MC_{}_GrandMother_PdgCode", acronym).c_str(), &GrandMother_PdgCode);
    }
};

// `Flat::MCInfo_LV_Mother` + `Neg` (`Flat::MCInfo_PxPyPz`) + `Pos` (`Flat::MCInfo_PxPyPz`) + `AtDecay` (`Flat::Coordinates`) + `IsHybrid`.
// NOTE: member functions will add the `MC_` prefix.
struct alignas(T2S_SIMD_ALIGN) MCInfo_V0 : Flat::MCInfo_LV_Mother {
    Flat::MCInfo_PxPyPz Neg;
    Flat::MCInfo_PxPyPz Pos;
    Flat::Coordinates AtDecay{};
    bool IsHybrid{};

    void CreateBranches_MCInfo_V0(TTree* tree, std::string_view acronym = "") {
        CreateBranches_MCInfo_LV_Mother(tree, acronym);
        Neg.CreateBranches_MCInfo_PxPyPz(tree, acronym);
        Pos.CreateBranches_MCInfo_PxPyPz(tree, acronym);
        AtDecay.CreateBranches_Coordinates(tree, std::format("MC_{}", acronym), "");
        tree->Branch(std::format("MC_{}_IsHybrid", acronym).c_str(), &IsHybrid);
    }
};

}  // namespace Tree2Secondaries::DF::Flat
