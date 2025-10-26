#ifndef T2S_DF_SEXAQUARK_HXX
#define T2S_DF_SEXAQUARK_HXX

#include <TTree.h>

#include "App/Utilities.hxx"
#include "DataFormats/DataFormats.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::DF {

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
        CreateBranches_State(tree);
        // event properties
        PV.CreateBranches_Coordinates(tree, "PV", "v");
        Utils::CreateBranch(tree, "RunNumber", &RunNumber);
        Utils::CreateBranch(tree, "DirNumber", &DirNumber);
        if (!is_mc) Utils::CreateBranch(tree, "DirNumberB", &DirNumberB);
        Utils::CreateBranch(tree, "EventNumber", &EventNumber);
        Utils::CreateBranch(tree, "MagneticField", &MagneticField);
        // fit info
        Utils::CreateBranch(tree, "Chi2NDF", &Chi2NDF);
        // extra info
        Utils::CreateBranch(tree, "E_MinusNucleon", &E_MinusNucleon);
        Utils::CreateBranch(tree, "AntiChannel", &AntiChannel);
    }
};

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
        SV.CreateBranches_Coordinates(tree, "MC_SV");
        // event properties
        PV.CreateBranches_Coordinates(tree, "MC_PV", "v");
        // reaction id + flags
        Utils::CreateBranch(tree, "ReactionID", &ReactionID);
        Utils::CreateBranch(tree, "IsSignal", &IsSignal);
        Utils::CreateBranch(tree, "IsHybrid", &IsHybrid);
    }
};

}  // namespace Tree2Secondaries::DF

#endif  // T2S_DF_SEXAQUARK_HXX
