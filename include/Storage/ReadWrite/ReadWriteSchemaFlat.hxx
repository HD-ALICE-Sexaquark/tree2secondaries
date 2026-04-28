#pragma once

#include <format>
#include <string_view>

#include <TTree.h>

#include "Storage/ReadWrite/ReadWritePrimitives.hxx"
#include "Storage/Schema/SchemaFlat.hxx"

namespace Tree2Secondaries::IO {

// Event //

// NOTE: used in `Packager` and `Finder`.
inline void ReadBranches(TTree* t, Schema::Flat::Event& event, bool is_mc) {
    Utils::ReadBranch(t, "RunNumber", &event.run_number);
    Utils::ReadBranch(t, "DirNumber", &event.dir_number);
    if (!is_mc) Utils::ReadBranch(t, "DirNumberB", &event.dir_number_b);
    Utils::ReadBranch(t, "EventNumber", &event.event_number);
    Utils::ReadBranch(t, "Centrality", &event.centrality);
    Utils::ReadBranch(t, "MagneticField", &event.magnetic_field);
    ReadBranches(t, event.pv, "PV");
    if (is_mc) ReadBranches(t, event.mc_pv, "MC_PV");
}

// NOTE: used in `Packager`.
inline void CreateBranches(TTree* t, Schema::Flat::Event& event, bool is_mc) {
    Utils::CreateBranch(t, "RunNumber", &event.run_number);
    Utils::CreateBranch(t, "DirNumber", &event.dir_number);
    if (!is_mc) Utils::CreateBranch(t, "DirNumberB", &event.dir_number_b);
    Utils::CreateBranch(t, "EventNumber", &event.event_number);
    Utils::CreateBranch(t, "Centrality", &event.centrality);
    Utils::CreateBranch(t, "MagneticField", &event.magnetic_field);
    CreateBranches(t, event.pv, "PV");
    if (is_mc) CreateBranches(t, event.mc_pv, "MC_PV");
}

// Injected //

inline void CreateBranches(TTree* t, Schema::Flat::Injected& inj) {
    Utils::CreateBranch(t, "RunNumber", &inj.run_number);
    Utils::CreateBranch(t, "DirNumber", &inj.dir_number);
    Utils::CreateBranch(t, "EventNumber", &inj.event_number);
    Utils::CreateBranch(t, "ReactionID", &inj.reaction_id);
    CreateBranches(t, inj.sv, "SV");
    CreateBranches(t, inj.lv, "Sexa");
    CreateBranches(t, inj.lv_nucleon, "Nucleon");
}

// Track //

inline void CreateBranches(TTree* t, Schema::Flat::Track& track, std::string_view prefix) {
    Utils::CreateBranch(t, std::format("{}_EsdEntry", prefix), &track.esd_entry);
    Utils::CreateBranch(t, std::format("{}_Charge", prefix), &track.charge);
    CreateBranches(t, track.state, prefix);
    Utils::CreateBranch(t, std::format("{}_PreDCAxy", prefix), &track.pre_dca_xy);
    Utils::CreateBranch(t, std::format("{}_PreDCAz", prefix), &track.pre_dca_z);
    Utils::CreateBranch(t, std::format("{}_TPC_Signal", prefix), &track.tpc_signal);
    Utils::CreateBranch(t, std::format("{}_NSigmaPion", prefix), &track.n_sigma_pion);
    Utils::CreateBranch(t, std::format("{}_NSigmaKaon", prefix), &track.n_sigma_kaon);
    Utils::CreateBranch(t, std::format("{}_NSigmaProton", prefix), &track.n_sigma_proton);
}

// MC_Track //

inline void CreateBranches(TTree* t, Schema::Flat::MC_Track& mc_track, bool ascendants_info, std::string_view acronym) {
    Utils::CreateBranch(t, std::format("MC_{}_McEntry", acronym), &mc_track.mc_entry);
    Utils::CreateBranch(t, std::format("MC_{}_PdgCode", acronym), &mc_track.pdg_code);
    if (ascendants_info) CreateBranches(t, mc_track.origin, std::format("MC_{}_Origin", acronym));
    CreateBranches(t, mc_track.lv, std::format("MC_{}", acronym));
    Utils::CreateBranch(t, std::format("MC_{}_IsTrue", acronym), &mc_track.is_true);
    Utils::CreateBranch(t, std::format("MC_{}_IsSecondary", acronym), &mc_track.is_secondary);
    Utils::CreateBranch(t, std::format("MC_{}_IsSignal", acronym), &mc_track.is_signal);
    Utils::CreateBranch(t, std::format("MC_{}_ReactionID", acronym), &mc_track.reaction_id);
    if (ascendants_info) {
        Utils::CreateBranch(t, std::format("MC_{}_Mother_McEntry", acronym), &mc_track.mother_mc_entry);
        Utils::CreateBranch(t, std::format("MC_{}_Mother_PdgCode", acronym), &mc_track.mother_pdg_code);
        Utils::CreateBranch(t, std::format("MC_{}_GM_McEntry", acronym), &mc_track.gm_mc_entry);
        Utils::CreateBranch(t, std::format("MC_{}_GM_PdgCode", acronym), &mc_track.gm_pdg_code);
    }
}

// V0 //

inline void CreateBranches(TTree* t, Schema::Flat::V0& v0, std::string_view acronym) {
    CreateBranches(t, v0.decay, std::format("{}_Decay", acronym));
    CreateBranches(t, v0.lv, acronym);
    Utils::CreateBranch(t, std::format("{}_Chi2NDF", acronym), &v0.chi2ndf);
    CreateBranches(t, v0.neg, std::format("{}_Neg", acronym));
    CreateBranches(t, v0.neg_pca_v0, std::format("{}_Neg_PCAwrtV0", acronym));
    CreateBranches(t, v0.pos, std::format("{}_Pos", acronym));
    CreateBranches(t, v0.pos_pca_v0, std::format("{}_Pos_PCAwrtV0", acronym));
}

// MC_V0 //

inline void CreateBranches(TTree* t, Schema::Flat::MC_V0& mc_v0, std::string_view acronym) {
    Utils::CreateBranch(t, std::format("MC_{}_McEntry", acronym), &mc_v0.mc_entry);
    Utils::CreateBranch(t, std::format("MC_{}_PdgCode", acronym), &mc_v0.pdg_code);
    CreateBranches(t, mc_v0.origin, std::format("MC_{}_Origin", acronym));
    CreateBranches(t, mc_v0.decay, std::format("MC_{}_Decay", acronym));
    CreateBranches(t, mc_v0.lv, std::format("MC_{}", acronym));
    Utils::CreateBranch(t, std::format("MC_{}_IsTrue", acronym), &mc_v0.is_true);
    Utils::CreateBranch(t, std::format("MC_{}_IsSecondary", acronym), &mc_v0.is_secondary);
    Utils::CreateBranch(t, std::format("MC_{}_IsSignal", acronym), &mc_v0.is_signal);
    Utils::CreateBranch(t, std::format("MC_{}_IsHybrid", acronym), &mc_v0.is_hybrid);
    Utils::CreateBranch(t, std::format("MC_{}_ReactionID", acronym), &mc_v0.reaction_id);
    Utils::CreateBranch(t, std::format("MC_{}_Mother_McEntry", acronym), &mc_v0.mother_mc_entry);
    Utils::CreateBranch(t, std::format("MC_{}_Mother_PdgCode", acronym), &mc_v0.mother_pdg_code);
    CreateBranches(t, mc_v0.neg, false, std::format("{}_Neg", acronym));
    CreateBranches(t, mc_v0.pos, false, std::format("{}_Pos", acronym));
}

// ChannelA //

inline void CreateBranches(TTree* t, Schema::Flat::ChannelA& sexa, bool is_mc) {
    CreateBranches(t, sexa.event, is_mc);
    CreateBranches(t, sexa.sv, "SV");
    CreateBranches(t, sexa.lv, "Sexa");
    Utils::CreateBranch(t, "Sexa_Chi2NDF", &sexa.chi2ndf);
    Utils::CreateBranch(t, "Sexa_E_MinusNucleon", &sexa.e_minus_nucleon);
    Utils::CreateBranch(t, "Sexa_IsAntiChannel", &sexa.is_anti);
    //
    CreateBranches(t, sexa.v0a, "V0A");
    CreateBranches(t, sexa.v0b, "V0B");
    CreateBranches(t, sexa.v0a_pca_sv, "V0A_PCAwrtSV");
    CreateBranches(t, sexa.v0b_pca_sv, "V0B_PCAwrtSV");
}

// MC_ChannelA //

inline void CreateBranches(TTree* t, Schema::Flat::MC_ChannelA& mc_sexa) {
    CreateBranches(t, mc_sexa.before, "MC_Before");
    CreateBranches(t, mc_sexa.after, "MC_After");
    CreateBranches(t, mc_sexa.nucleon, "MC_Nucleon");
    CreateBranches(t, mc_sexa.sv, "MC_SV");
    Utils::CreateBranch(t, "MC_ReactionID", &mc_sexa.reaction_id);
    Utils::CreateBranch(t, "MC_IsSignal", &mc_sexa.is_signal);
    Utils::CreateBranch(t, "MC_IsHybrid", &mc_sexa.is_hybrid);
    //
    CreateBranches(t, mc_sexa.v0a, "V0A");
    CreateBranches(t, mc_sexa.v0b, "V0B");
}

// ChannelD //

inline void CreateBranches(TTree* t, Schema::Flat::ChannelD& sexa, bool is_mc) {
    CreateBranches(t, sexa.event, is_mc);
    CreateBranches(t, sexa.sv, "SV");
    CreateBranches(t, sexa.lv, "Sexa");
    Utils::CreateBranch(t, "Sexa_Chi2NDF", &sexa.chi2ndf);
    Utils::CreateBranch(t, "Sexa_E_MinusNucleon", &sexa.e_minus_nucleon);
    Utils::CreateBranch(t, "Sexa_IsAntiChannel", &sexa.is_anti);
    //
    CreateBranches(t, sexa.v0, "V0");
    CreateBranches(t, sexa.kaon, "Kaon");
    CreateBranches(t, sexa.v0_pca_sv, "V0_PCAwrtSV");
    CreateBranches(t, sexa.kaon_pca_sv, "Kaon_PCAwrtSV");
}

// MC_ChannelD //

inline void CreateBranches(TTree* t, Schema::Flat::MC_ChannelD& mc_sexa) {
    CreateBranches(t, mc_sexa.before, "MC_Before");
    CreateBranches(t, mc_sexa.after, "MC_After");
    CreateBranches(t, mc_sexa.nucleon, "MC_Nucleon");
    CreateBranches(t, mc_sexa.sv, "MC_SV");
    Utils::CreateBranch(t, "MC_ReactionID", &mc_sexa.reaction_id);
    Utils::CreateBranch(t, "MC_IsSignal", &mc_sexa.is_signal);
    Utils::CreateBranch(t, "MC_IsHybrid", &mc_sexa.is_hybrid);
    //
    CreateBranches(t, mc_sexa.v0, "V0");
    CreateBranches(t, mc_sexa.kaon, true, "Kaon");
}

}  // namespace Tree2Secondaries::IO
