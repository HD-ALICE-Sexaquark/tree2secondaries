#pragma once

#include "Math/Constants.hxx"
#include "Schema/SchemaFlatPrimitives.hxx"
#include "Schema/SchemaFlatTrack.hxx"

class TTree;

namespace Tree2Secondaries::Schema {

struct alignas(T2S_SIMD_ALIGN) FlatV0 {
    // Storage //
    void CreateBranches(TTree* t, std::string_view acronym) {
        decay.CreateBranches(t, std::format("{}_Decay", acronym));
        lv.CreateBranches(t, acronym);
        Utils::CreateBranch(t, std::format("{}_Chi2NDF", acronym), &chi2ndf);
        neg.CreateBranches(t, std::format("{}_Neg", acronym));
        neg_pca_v0.CreateBranches(t, std::format("{}_Neg_PCAwrtV0", acronym));
        pos.CreateBranches(t, std::format("{}_Pos", acronym));
        pos_pca_v0.CreateBranches(t, std::format("{}_Pos_PCAwrtV0", acronym));
    }

    // Member Variables //
    Schema::FlatTrack neg{};
    Schema::FlatTrack pos{};
    Schema::State6 neg_pca_v0{};
    Schema::State6 pos_pca_v0{};
    Schema::Mom4 lv{};
    Schema::Coordinates decay{};
    float chi2ndf{Const::DummyFloat};
};

struct alignas(T2S_SIMD_ALIGN) FlatMcV0 {
    // Storage //
    void CreateBranches(TTree* t, std::string_view acronym) {
        Utils::CreateBranch(t, std::format("MC_{}_McEntry", acronym), &mc_entry);
        Utils::CreateBranch(t, std::format("MC_{}_PdgCode", acronym), &pdg_code);
        origin.CreateBranches(t, std::format("MC_{}_Origin", acronym));
        decay.CreateBranches(t, std::format("MC_{}_Decay", acronym));
        lv.CreateBranches(t, std::format("MC_{}", acronym));
        Utils::CreateBranch(t, std::format("MC_{}_IsTrue", acronym), &is_true);
        Utils::CreateBranch(t, std::format("MC_{}_IsSecondary", acronym), &is_secondary);
        Utils::CreateBranch(t, std::format("MC_{}_IsSignal", acronym), &is_signal);
        Utils::CreateBranch(t, std::format("MC_{}_IsHybrid", acronym), &is_hybrid);
        Utils::CreateBranch(t, std::format("MC_{}_ReactionID", acronym), &reaction_id);
        Utils::CreateBranch(t, std::format("MC_{}_Mother_McEntry", acronym), &mother_mc_entry);
        Utils::CreateBranch(t, std::format("MC_{}_Mother_PdgCode", acronym), &mother_pdg_code);
    }

    // Member Variables //
    Schema::Mom4 lv{};
    Schema::Coordinates origin{};
    Schema::Coordinates decay{};
    int mc_entry{Const::DummyInt};  // NOTE: negative when not found
    int pdg_code{Const::DummyInt};
    int reaction_id{Const::DummyInt};
    int mother_mc_entry{Const::DummyInt};
    int mother_pdg_code{Const::DummyInt};
    bool is_true{false};
    bool is_secondary{false};
    bool is_signal{false};
    bool is_hybrid{false};
};

}  // namespace Tree2Secondaries::Schema
