#pragma once

#include <TTree.h>

#include "DataFormats/Flat.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::DF::Found {

// `Flat::State` + `PV` + `RunNumber` + `DirNumber` + `DirNumberB` + `EventNumber` + `MagneticField` + `Chi2NDF` + `E_MinusNucleon` + `AntiChannel`
struct alignas(T2S_SIMD_ALIGN) Sexaquark : Flat::State {
    // event properties
    Flat::Coordinates PV{};
    unsigned int RunNumber{};
    unsigned int DirNumber{};
    unsigned int DirNumberB{};  // NOTE: not used when analyzing MC
    unsigned int EventNumber{};
    float MagneticField{};
    // fit info
    float Chi2NDF{};
    // extra info
    float E_MinusNucleon{};
    bool AntiChannel{};

    void CreateBranches_Sexaquark(TTree* tree, bool is_mc) {
        CreateBranches_State(tree, "Sexa");
        // event properties
        PV.CreateBranches_Coordinates(tree, "PV", "v");
        tree->Branch("RunNumber", &RunNumber);
        tree->Branch("DirNumber", &DirNumber);
        if (!is_mc) tree->Branch("DirNumberB", &DirNumberB);
        tree->Branch("EventNumber", &EventNumber);
        tree->Branch("MagneticField", &MagneticField);
        // fit info
        tree->Branch("Chi2NDF", &Chi2NDF);
        // extra info
        tree->Branch("E_MinusNucleon", &E_MinusNucleon);
        tree->Branch("AntiChannel", &AntiChannel);
    }
};

// `Before` (`Flat::LorentzVector`) + `After` (`Flat::LorentzVector`) + `Nucleon` (`Flat::LorentzVector`) + `SV` (`Flat::Coordinates`) + `PV`
// (`Flat::Coordinates`) + `ReactionID`  + `IsSignal`  + `IsHybrid`
struct alignas(T2S_SIMD_ALIGN) MC_Sexaquark {
    Flat::LorentzVector Before{};
    Flat::LorentzVector After{};
    Flat::LorentzVector Nucleon{};
    // secondary vertex
    Flat::Coordinates SV{};
    // event properties
    Flat::Coordinates PV{};
    // reaction id + flags
    int ReactionID{};
    bool IsSignal{};
    bool IsHybrid{};

    void CreateBranches_MC_Sexaquark(TTree* tree) {
        Before.CreateBranches_LorentzVector(tree, "MC_Before");
        After.CreateBranches_LorentzVector(tree, "MC_After");
        Nucleon.CreateBranches_LorentzVector(tree, "MC_Nucleon");
        // secondary vertex
        SV.CreateBranches_Coordinates(tree, "MC_SV", "");
        // event properties
        PV.CreateBranches_Coordinates(tree, "MC_PV", "v");
        // reaction id + flags
        tree->Branch("ReactionID", &ReactionID);
        tree->Branch("IsSignal", &IsSignal);
        tree->Branch("IsHybrid", &IsHybrid);
    }
};

}  // namespace Tree2Secondaries::DF::Found
