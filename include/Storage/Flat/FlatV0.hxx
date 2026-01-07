#pragma once

#include "Math/Constants.hxx"
#include "Storage/Flat/BaseFlat.hxx"
#include "Storage/Flat/FlatTrack.hxx"

namespace Tree2Secondaries::Storage::Flat {

// `Flat::State` +
// `Neg` (`Flat::Track`) + `Pos` (`Flat::Track`) +
// `Neg_atV0` (`Flat::State_NoE`) + `Pos_atV0` (`Flat::State_NoE`) +
// `Chi2NDF`.
struct alignas(T2S_SIMD_ALIGN) V0 : Flat::State {
    Flat::Track Neg{};
    Flat::Track Pos{};
    Flat::State_NoE Neg_atV0{};
    Flat::State_NoE Pos_atV0{};
    float Chi2NDF{};

    void CreateBranches_FlatV0(TTree* tree, std::string_view prefix) {
        CreateBranches_FlatState(tree, prefix);
        Neg.CreateBranches_FlatTrack(tree, std::format("{}_Neg", prefix));
        Pos.CreateBranches_FlatTrack(tree, std::format("{}_Pos", prefix));
        Neg_atV0.CreateBranches_FlatState_NoE(tree, std::format("{}_Neg_atV0", prefix));
        Pos_atV0.CreateBranches_FlatState_NoE(tree, std::format("{}_Pos_atV0", prefix));
        tree->Branch(std::format("{}_Chi2NDF", prefix).c_str(), &Chi2NDF);
    }
};

// `Flat::MC` + `Flat::LV` +
// `Neg` (`Flat::MC`) + `Pos` (`Flat::MC`) +
// `Mother` (`Flat::MC_Id`) +
// `Neg_Momentum` (`Flat::PxPyPz`) + `Pos_Momentum` (`Flat::PxPyPz`) +
// `AtDecay` (`Flat::Coordinates`) +
// `IsHybrid`.
struct alignas(T2S_SIMD_ALIGN) MC_V0 : Flat::MC, Flat::LorentzVector {
    Flat::MC Neg{};
    Flat::MC Pos{};
    Flat::MC_Id Mother{};
    Flat::PxPyPz Neg_Momentum{};
    Flat::PxPyPz Pos_Momentum{};
    Flat::Coordinates AtDecay{};
    bool IsHybrid{};

    // NOTE: member functions will add the `MC_` prefix.
    void CreateBranches_FlatMC_V0(TTree* tree, std::string_view acronym = "") {
        CreateBranches_FlatMC(tree, acronym);
        CreateBranches_FlatLV(tree, std::format("MC_{}", acronym));
        Neg.CreateBranches_FlatMC(tree, std::format("{}_Neg", acronym));
        Pos.CreateBranches_FlatMC(tree, std::format("{}_Pos", acronym));
        Mother.CreateBranches_FlatMC_Id(tree, std::format("{}_Mother", acronym));
        Neg_Momentum.CreateBranches_FlatPxPyPz(tree, std::format("MC_{}_Neg", acronym));
        Pos_Momentum.CreateBranches_FlatPxPyPz(tree, std::format("MC_{}_Pos", acronym));
        AtDecay.CreateBranches_FlatCoordinates(tree, std::format("MC_{}_atDecay", acronym), "");
        tree->Branch(std::format("MC_{}_IsHybrid", acronym).c_str(), &IsHybrid);
    }
};

}  // namespace Tree2Secondaries::Storage::Flat
