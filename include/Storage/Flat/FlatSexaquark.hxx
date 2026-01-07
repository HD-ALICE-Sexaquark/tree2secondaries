#pragma once

#include <TTree.h>

#include "Math/Constants.hxx"
#include "Storage/Flat/BaseFlat.hxx"
#include "Storage/Flat/FlatEvent.hxx"

namespace Tree2Secondaries::Storage::Flat {

// `Flat::State` +
// `Event` (`Flat::Event`) +
// `Chi2NDF` + `E_MinusNucleon` + `AntiChannel`.
struct alignas(T2S_SIMD_ALIGN) Sexaquark : Flat::State {
    Flat::Event Event{};
    // fit info
    float Chi2NDF{};
    // extra info
    float E_MinusNucleon{};
    bool AntiChannel{};

    void CreateBranches_FlatSexaquark(TTree* tree, bool is_mc) {
        CreateBranches_FlatState(tree, "Sexa");
        Event.CreateBranches_FlatEvent(tree, is_mc);
        // fit info
        tree->Branch("Chi2NDF", &Chi2NDF);
        // extra info
        tree->Branch("E_MinusNucleon", &E_MinusNucleon);
        tree->Branch("AntiChannel", &AntiChannel);
    }
};

// `Before` (`Flat::LV`) + `After` (`Flat::LV`) + `Nucleon` (`Flat::LV`) +
// `PV` (`Flat::Coordinates`) + `SV` (`Flat::Coordinates`) +
// `ReactionID`  + `IsSignal`  + `IsHybrid`.
struct alignas(T2S_SIMD_ALIGN) MC_Sexaquark {
    Flat::LorentzVector Before{};
    Flat::LorentzVector After{};
    Flat::LorentzVector Nucleon{};
    Flat::Coordinates PV{};
    Flat::Coordinates SV{};
    int ReactionID{};
    bool IsSignal{};
    bool IsHybrid{};

    void CreateBranches_FlatMC_Sexaquark(TTree* tree) {
        Before.CreateBranches_FlatLV(tree, "MC_Before");
        After.CreateBranches_FlatLV(tree, "MC_After");
        Nucleon.CreateBranches_FlatLV(tree, "MC_Nucleon");
        // secondary vertex
        PV.CreateBranches_FlatCoordinates(tree, "MC_PV", "v");
        SV.CreateBranches_FlatCoordinates(tree, "MC_SV", "");
        // reaction id + flags
        tree->Branch("ReactionID", &ReactionID);
        tree->Branch("IsSignal", &IsSignal);
        tree->Branch("IsHybrid", &IsHybrid);
    }
};

}  // namespace Tree2Secondaries::Storage::Flat
