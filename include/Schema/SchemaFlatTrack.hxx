#pragma once

#include "Math/Constants.hxx"
#include "Schema/SchemaFlatPrimitives.hxx"

class TTree;

namespace Tree2Secondaries::Schema {

struct alignas(T2S_SIMD_ALIGN) FlatTrack {
    // Storage //
    void CreateBranches(TTree* t, std::string_view prefix) {
        Utils::CreateBranch(t, std::format("{}_EsdEntry", prefix), &esd_entry);
        Utils::CreateBranch(t, std::format("{}_Charge", prefix), &charge);
        state.CreateBranches(t, prefix);
        Utils::CreateBranch(t, std::format("{}_PreDCAxy", prefix), &pre_dca_xy);
        Utils::CreateBranch(t, std::format("{}_PreDCAz", prefix), &pre_dca_z);
        Utils::CreateBranch(t, std::format("{}_TPC_Signal", prefix), &tpc_signal);
        Utils::CreateBranch(t, std::format("{}_NSigmaPion", prefix), &n_sigma_pion);
        Utils::CreateBranch(t, std::format("{}_NSigmaKaon", prefix), &n_sigma_kaon);
        Utils::CreateBranch(t, std::format("{}_NSigmaProton", prefix), &n_sigma_proton);
    }

    // Member Variables //
    Schema::State6 state{};
    unsigned int esd_entry{0};
    char charge{0};
    float pre_dca_xy{Const::DummyFloat};
    float pre_dca_z{Const::DummyFloat};
    float tpc_signal{Const::DummyFloat};
    float n_sigma_pion{Const::DummyFloat};
    float n_sigma_kaon{Const::DummyFloat};
    float n_sigma_proton{Const::DummyFloat};
};

struct alignas(T2S_SIMD_ALIGN) FlatMcTrack {
    // Storage //
    void CreateBranches(TTree* t, bool ascendants_info, std::string_view acronym) {
        Utils::CreateBranch(t, std::format("MC_{}_McEntry", acronym), &mc_entry);
        Utils::CreateBranch(t, std::format("MC_{}_PdgCode", acronym), &pdg_code);
        if (ascendants_info) origin.CreateBranches(t, std::format("MC_{}_Origin", acronym));
        lv.CreateBranches(t, std::format("MC_{}", acronym));
        Utils::CreateBranch(t, std::format("MC_{}_IsTrue", acronym), &is_true);
        Utils::CreateBranch(t, std::format("MC_{}_IsSecondary", acronym), &is_secondary);
        Utils::CreateBranch(t, std::format("MC_{}_IsSignal", acronym), &is_signal);
        Utils::CreateBranch(t, std::format("MC_{}_ReactionID", acronym), &reaction_id);
        if (ascendants_info) {
            Utils::CreateBranch(t, std::format("MC_{}_Mother_McEntry", acronym), &mother_mc_entry);
            Utils::CreateBranch(t, std::format("MC_{}_Mother_PdgCode", acronym), &mother_pdg_code);
            Utils::CreateBranch(t, std::format("MC_{}_GM_McEntry", acronym), &gm_mc_entry);
            Utils::CreateBranch(t, std::format("MC_{}_GM_PdgCode", acronym), &gm_pdg_code);
        }
    }

    // Member Variables //
    Schema::Mom4 lv{};
    Schema::Coordinates origin{};
    int mc_entry{Const::DummyInt};  // NOTE: negative when not found
    int pdg_code{Const::DummyInt};
    int reaction_id{Const::DummyInt};
    bool is_true{false};
    bool is_signal{false};
    bool is_secondary{false};
    int mother_mc_entry{Const::DummyInt};
    int mother_pdg_code{Const::DummyInt};
    int gm_mc_entry{Const::DummyInt};
    int gm_pdg_code{Const::DummyInt};
};

}  // namespace Tree2Secondaries::Schema
