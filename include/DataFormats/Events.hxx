#pragma once

#include <vector>

#include <TTree.h>

#include "App/Utilities.hxx"
#include "DataFormats/Flat.hxx"
#include "DataFormats/StructsOfVectors.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::DF::Events {

// `PV` (`Flat::Coordinates`) + `MC_PV` (`Flat::Coordinates`) + `RunNumber` + `DirNumber` + `DirNumberB` + `EventNumber` + `Centrality` +
// `MagneticField`.
struct alignas(T2S_SIMD_ALIGN) Event {
    Flat::Coordinates PV{};
    Flat::Coordinates MC_PV{};  // NOTE: only read when analyzing MC
    unsigned int RunNumber{0};
    unsigned int DirNumber{0};
    unsigned int DirNumberB{0};  // NOTE: only read when analyzing RD
    unsigned int EventNumber{0};
    float Centrality{0.};
    float MagneticField{0.};

    void CreateBranches_Event(TTree *tree, bool is_mc) {
        PV.CreateBranches_Coordinates(tree, "PV", "v");
        if (is_mc) MC_PV.CreateBranches_Coordinates(tree, "MC_PV", "v");
        tree->Branch("RunNumber", &RunNumber);
        tree->Branch("DirNumber", &DirNumber);
        if (!is_mc) tree->Branch("DirNumberB", &DirNumberB);
        tree->Branch("EventNumber", &EventNumber);
        tree->Branch("Centrality", &Centrality);
        tree->Branch("MagneticField", &MagneticField);
    }
    void ReadBranches_Event(TTree *tree, bool is_mc) {
        PV.ReadBranches_Coordinates(tree, "PV", "v");
        if (is_mc) MC_PV.ReadBranches_Coordinates(tree, "MC_PV", "v");
        Utils::ReadBranch(tree, "RunNumber", &RunNumber);
        Utils::ReadBranch(tree, "DirNumber", &DirNumber);
        if (!is_mc) Utils::ReadBranch(tree, "DirNumberB", &DirNumberB);
        Utils::ReadBranch(tree, "EventNumber", &EventNumber);
        Utils::ReadBranch(tree, "Centrality", &Centrality);
        Utils::ReadBranch(tree, "MagneticField", &MagneticField);
    }
};

// `SOV::States` + `PdgCode` + `MotherEntry` + `Status` + `Generator` + `IsPrimary` + `IsSecFromMat` + `IsSecFromWeak`.
struct alignas(T2S_SIMD_ALIGN) MCParticles : SOV::States {
    std::vector<int> *PdgCode{nullptr};
    std::vector<int> *MotherEntry{nullptr};
    std::vector<int> *Status{nullptr};
    std::vector<int> *Generator{nullptr};
    std::vector<char> *IsPrimary{nullptr};
    std::vector<char> *IsSecFromMat{nullptr};
    std::vector<char> *IsSecFromWeak{nullptr};

    void ReadBranches_MCParticles(TTree *tree) {
        ReadBranches_States(tree, "MC");
        Utils::ReadBranch(tree, "MC_PdgCode", &PdgCode);
        Utils::ReadBranch(tree, "MC_Mother_McEntry", &MotherEntry);
        Utils::ReadBranch(tree, "MC_Status", &Status);
        Utils::ReadBranch(tree, "MC_Generator", &Generator);
        Utils::ReadBranch(tree, "MC_IsPrimary", &IsPrimary);
        Utils::ReadBranch(tree, "MC_IsSecFromMat", &IsSecFromMat);
        Utils::ReadBranch(tree, "MC_IsSecFromWeak", &IsSecFromWeak);
    }
};

// `SOV::States_NoE` + `SOV::CovMatrices_NoE` + `Charge` + `DCAxy` + `DCAz` + `TPCSignal` + `NSigmaPion` + `NSigmaKaon` + `NSigmaProton` + `McEntry`.
struct alignas(T2S_SIMD_ALIGN) Tracks : SOV::States_NoE, SOV::CovMatrices_NoE {
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

    void ReadBranches_MCParticles(TTree *tree, bool is_mc) {
        ReadBranches_States_NoE(tree, "Track");
        ReadBranches_CovMatrices_NoE(tree, "Track");
        Utils::ReadBranch(tree, "Track_Charge", &Charge);
        Utils::ReadBranch(tree, "Track_DCAxy", &DCAxy);
        Utils::ReadBranch(tree, "Track_DCAz", &DCAz);
        // -- pid info
        Utils::ReadBranch(tree, "Track_TPCSignal", &TPCSignal);
        Utils::ReadBranch(tree, "Track_NSigmaPion", &NSigmaPion);
        Utils::ReadBranch(tree, "Track_NSigmaKaon", &NSigmaKaon);
        Utils::ReadBranch(tree, "Track_NSigmaProton", &NSigmaProton);
        // -- mc info
        if (is_mc) Utils::ReadBranch(tree, "Track_McEntry", &McEntry);
    }
};

}  // namespace Tree2Secondaries::DF::Events
