#pragma once

#include <format>
#include <string_view>

#include <TTree.h>

#include "App/Utilities.hxx"
#include "Storage/ReadWrite/ReadWritePrimitives.hxx"
#include "Storage/Schema/SchemaVector.hxx"

namespace Tree2Secondaries::IO {

// MCParticles //

inline void ReadBranches(TTree* t, Schema::Vector::MCParticles& mc) {
    Utils::ReadBranch(t, "MC_PdgCode", &mc.pdg_code);
    Utils::ReadBranch(t, "MC_Charge", &mc.charge);
    Utils::ReadBranch(t, "MC_Mother_McEntry", &mc.mother_mc_entry);
    Utils::ReadBranch(t, "MC_N_Daughters", &mc.n_daughters);
    Utils::ReadBranch(t, "MC_FirstDau_McEntry", &mc.firstdau_mc_entry);
    Utils::ReadBranch(t, "MC_LastDau_McEntry", &mc.lastdau_mc_entry);
    ReadBranches(t, mc.origin, "MC_Origin");
    ReadBranches(t, mc.lv, "MC");
    Utils::ReadBranch(t, "MC_Status", &mc.status);
    Utils::ReadBranch(t, "MC_Generator", &mc.generator);
    Utils::ReadBranch(t, "MC_IsPhysPrimary", &mc.is_phys_primary);
    Utils::ReadBranch(t, "MC_IsSecFromMat", &mc.is_sec_from_mat);
    Utils::ReadBranch(t, "MC_IsSecFromWeak", &mc.is_sec_from_weak);
}

// Injected //

inline void ReadBranches(TTree* t, Schema::Vector::Injected& inj, bool include_sv) {
    Utils::ReadBranch(t, "ReactionID", &inj.reaction_id);
    ReadBranches(t, inj.mom, "Sexa");
    if (include_sv) ReadBranches(t, inj.sv, "SV");
    ReadBranches(t, inj.mom_nucleon, "Nucleon");
}

inline void CreateBranches(TTree* t, Schema::Vector::Injected& inj, bool include_sv) {
    Utils::CreateBranch(t, "ReactionID", &inj.reaction_id);
    CreateBranches(t, inj.mom, "Sexa");
    if (include_sv) CreateBranches(t, inj.sv, "SV");
    CreateBranches(t, inj.mom_nucleon, "Nucleon");
}

inline void ClearBranches(Schema::Vector::Injected& inj, bool include_sv) {
    inj.reaction_id->clear();
    ClearBranches(inj.mom);
    if (include_sv) ClearBranches(inj.sv);
    ClearBranches(inj.mom_nucleon);
}

// Tracks //

inline void ReadBranches(TTree* t, Schema::Vector::Tracks& track, bool include_mc_entry, std::string_view prefix = "Track") {
    Utils::ReadBranch(t, std::format("{}_EsdEntry", prefix), &track.esd_entry);
    Utils::ReadBranch(t, std::format("{}_Charge", prefix), &track.charge);
    ReadBranches(t, track.state, prefix);
    Utils::ReadBranch(t, std::format("{}_PreDCAxy", prefix), &track.pre_dca_xy);
    Utils::ReadBranch(t, std::format("{}_PreDCAz", prefix), &track.pre_dca_z);
    Utils::ReadBranch(t, std::format("{}_TPC_Signal", prefix), &track.tpc_signal);
    Utils::ReadBranch(t, std::format("{}_NSigmaPion", prefix), &track.n_sigma_pion);
    Utils::ReadBranch(t, std::format("{}_NSigmaKaon", prefix), &track.n_sigma_kaon);
    Utils::ReadBranch(t, std::format("{}_NSigmaProton", prefix), &track.n_sigma_proton);
    ReadBranches(t, track.cov, prefix);
    if (include_mc_entry) Utils::ReadBranch(t, std::format("{}_McEntry", prefix), &track.mc_entry);
}

inline void CreateBranches(TTree* t, Schema::Vector::Tracks& track, bool include_mc_entry, std::string_view prefix = "Track") {
    Utils::CreateBranch(t, std::format("{}_EsdEntry", prefix), &track.esd_entry);
    Utils::CreateBranch(t, std::format("{}_Charge", prefix), &track.charge);
    CreateBranches(t, track.state, prefix);
    Utils::CreateBranch(t, std::format("{}_PreDCAxy", prefix), &track.pre_dca_xy);
    Utils::CreateBranch(t, std::format("{}_PreDCAz", prefix), &track.pre_dca_z);
    Utils::CreateBranch(t, std::format("{}_TPC_Signal", prefix), &track.tpc_signal);
    Utils::CreateBranch(t, std::format("{}_NSigmaPion", prefix), &track.n_sigma_pion);
    Utils::CreateBranch(t, std::format("{}_NSigmaKaon", prefix), &track.n_sigma_kaon);
    Utils::CreateBranch(t, std::format("{}_NSigmaProton", prefix), &track.n_sigma_proton);
    CreateBranches(t, track.cov, prefix);
    if (include_mc_entry) Utils::CreateBranch(t, std::format("{}_McEntry", prefix), &track.mc_entry);
}

inline void ClearBranches(Schema::Vector::Tracks& track, bool include_mc_entry) {
    track.esd_entry->clear();
    track.charge->clear();
    ClearBranches(track.state);
    track.pre_dca_xy->clear();
    track.pre_dca_z->clear();
    track.tpc_signal->clear();
    track.n_sigma_pion->clear();
    track.n_sigma_kaon->clear();
    track.n_sigma_proton->clear();
    ClearBranches(track.cov);
    if (include_mc_entry) track.mc_entry->clear();
}

// MC_Tracks //

inline void ReadBranches(TTree* t, Schema::Vector::MC_Tracks& mc_track, bool ascendants_info, std::string_view acronym) {
    Utils::ReadBranch(t, std::format("MC_{}_McEntry", acronym), &mc_track.mc_entry);
    Utils::ReadBranch(t, std::format("MC_{}_PdgCode", acronym), &mc_track.pdg_code);
    ReadBranches(t, mc_track.lv, std::format("MC_{}", acronym));
    Utils::ReadBranch(t, std::format("MC_{}_IsTrue", acronym), &mc_track.is_true);
    Utils::ReadBranch(t, std::format("MC_{}_IsSecondary", acronym), &mc_track.is_secondary);
    Utils::ReadBranch(t, std::format("MC_{}_IsSignal", acronym), &mc_track.is_signal);
    Utils::ReadBranch(t, std::format("MC_{}_ReactionID", acronym), &mc_track.reaction_id);
    if (!ascendants_info) return;
    ReadBranches(t, mc_track.origin, std::format("MC_{}_Origin", acronym));
    Utils::ReadBranch(t, std::format("MC_{}_Mother_McEntry", acronym), &mc_track.mother_mc_entry);
    Utils::ReadBranch(t, std::format("MC_{}_Mother_PdgCode", acronym), &mc_track.mother_pdg_code);
    Utils::ReadBranch(t, std::format("MC_{}_GM_McEntry", acronym), &mc_track.gm_mc_entry);
    Utils::ReadBranch(t, std::format("MC_{}_GM_PdgCode", acronym), &mc_track.gm_pdg_code);
}

inline void CreateBranches(TTree* t, Schema::Vector::MC_Tracks& mc_track, bool ascendants_info, std::string_view acronym) {
    Utils::CreateBranch(t, std::format("MC_{}_McEntry", acronym), &mc_track.mc_entry);
    Utils::CreateBranch(t, std::format("MC_{}_PdgCode", acronym), &mc_track.pdg_code);
    CreateBranches(t, mc_track.lv, std::format("MC_{}", acronym));
    Utils::CreateBranch(t, std::format("MC_{}_IsTrue", acronym), &mc_track.is_true);
    Utils::CreateBranch(t, std::format("MC_{}_IsSecondary", acronym), &mc_track.is_secondary);
    Utils::CreateBranch(t, std::format("MC_{}_IsSignal", acronym), &mc_track.is_signal);
    Utils::CreateBranch(t, std::format("MC_{}_ReactionID", acronym), &mc_track.reaction_id);
    if (!ascendants_info) return;
    CreateBranches(t, mc_track.origin, std::format("MC_{}_Origin", acronym));
    Utils::CreateBranch(t, std::format("MC_{}_Mother_McEntry", acronym), &mc_track.mother_mc_entry);
    Utils::CreateBranch(t, std::format("MC_{}_Mother_PdgCode", acronym), &mc_track.mother_pdg_code);
    Utils::CreateBranch(t, std::format("MC_{}_GM_McEntry", acronym), &mc_track.gm_mc_entry);
    Utils::CreateBranch(t, std::format("MC_{}_GM_PdgCode", acronym), &mc_track.gm_pdg_code);
}

inline void ClearBranches(Schema::Vector::MC_Tracks& mc_track, bool ascendants_info) {
    mc_track.mc_entry->clear();
    mc_track.pdg_code->clear();
    ClearBranches(mc_track.lv);
    mc_track.is_true->clear();
    mc_track.is_secondary->clear();
    mc_track.is_signal->clear();
    mc_track.reaction_id->clear();
    if (!ascendants_info) return;
    ClearBranches(mc_track.origin);
    mc_track.mother_mc_entry->clear();
    mc_track.mother_pdg_code->clear();
    mc_track.gm_mc_entry->clear();
    mc_track.gm_pdg_code->clear();
}

// V0s //

inline void ReadBranches(TTree* t, Schema::Vector::V0s& v0, std::string_view acronym) {
    ReadBranches(t, v0.decay, std::format("{}_Decay", acronym));
    ReadBranches(t, v0.lv, acronym);
    Utils::ReadBranch(t, std::format("{}_Chi2NDF", acronym), &v0.chi2ndf);
    ReadBranches(t, v0.cov, acronym);
    ReadBranches(t, v0.neg, false, std::format("{}_Neg", acronym));
    ReadBranches(t, v0.neg_pca_v0, std::format("{}_Neg_PCAwrtV0", acronym));
    ReadBranches(t, v0.pos, false, std::format("{}_Pos", acronym));
    ReadBranches(t, v0.pos_pca_v0, std::format("{}_Pos_PCAwrtV0", acronym));
}

inline void CreateBranches(TTree* t, Schema::Vector::V0s& v0, std::string_view acronym) {
    CreateBranches(t, v0.decay, std::format("{}_Decay", acronym));
    CreateBranches(t, v0.lv, acronym);
    Utils::CreateBranch(t, std::format("{}_Chi2NDF", acronym), &v0.chi2ndf);
    CreateBranches(t, v0.cov, acronym);
    CreateBranches(t, v0.neg, false, std::format("{}_Neg", acronym));
    CreateBranches(t, v0.neg_pca_v0, std::format("{}_Neg_PCAwrtV0", acronym));
    CreateBranches(t, v0.pos, false, std::format("{}_Pos", acronym));
    CreateBranches(t, v0.pos_pca_v0, std::format("{}_Pos_PCAwrtV0", acronym));
}

inline void ClearBranches(Schema::Vector::V0s& v0) {
    ClearBranches(v0.decay);
    ClearBranches(v0.lv);
    v0.chi2ndf->clear();
    ClearBranches(v0.cov);
    ClearBranches(v0.neg, false);
    ClearBranches(v0.neg_pca_v0);
    ClearBranches(v0.pos, false);
    ClearBranches(v0.pos_pca_v0);
}

// MC_V0s //

inline void ReadBranches(TTree* t, Schema::Vector::MC_V0s& mc_v0, std::string_view acronym) {
    Utils::ReadBranch(t, std::format("MC_{}_McEntry", acronym), &mc_v0.mc_entry);
    Utils::ReadBranch(t, std::format("MC_{}_PdgCode", acronym), &mc_v0.pdg_code);
    ReadBranches(t, mc_v0.origin, std::format("MC_{}_Origin", acronym));
    ReadBranches(t, mc_v0.decay, std::format("MC_{}_Decay", acronym));
    ReadBranches(t, mc_v0.lv, std::format("MC_{}", acronym));
    Utils::ReadBranch(t, std::format("MC_{}_IsTrue", acronym), &mc_v0.is_true);
    Utils::ReadBranch(t, std::format("MC_{}_IsSecondary", acronym), &mc_v0.is_secondary);
    Utils::ReadBranch(t, std::format("MC_{}_IsSignal", acronym), &mc_v0.is_signal);
    Utils::ReadBranch(t, std::format("MC_{}_IsHybrid", acronym), &mc_v0.is_hybrid);
    Utils::ReadBranch(t, std::format("MC_{}_ReactionID", acronym), &mc_v0.reaction_id);
    Utils::ReadBranch(t, std::format("MC_{}_Mother_McEntry", acronym), &mc_v0.mother_mc_entry);
    Utils::ReadBranch(t, std::format("MC_{}_Mother_PdgCode", acronym), &mc_v0.mother_pdg_code);
    ReadBranches(t, mc_v0.neg, false, std::format("{}_Neg", acronym));
    ReadBranches(t, mc_v0.pos, false, std::format("{}_Pos", acronym));
}

inline void CreateBranches(TTree* t, Schema::Vector::MC_V0s& mc_v0, std::string_view acronym) {
    Utils::CreateBranch(t, std::format("MC_{}_McEntry", acronym), &mc_v0.mc_entry);
    Utils::CreateBranch(t, std::format("MC_{}_PdgCode", acronym), &mc_v0.pdg_code);
    CreateBranches(t, mc_v0.origin, std::format("MC_{}_Origin", acronym));
    CreateBranches(t, mc_v0.decay, std::format("MC_{}_Decay", acronym));
    CreateBranches(t, mc_v0.lv, std::format("MC_{}", acronym));
    Utils::CreateBranch(t, std::format("MC_{}_IsTrue", acronym), &mc_v0.is_true);
    Utils::CreateBranch(t, std::format("MC_{}_IsSecondary", acronym), &mc_v0.is_secondary);
    Utils::CreateBranch(t, std::format("MC_{}_IsHybrid", acronym), &mc_v0.is_hybrid);
    Utils::CreateBranch(t, std::format("MC_{}_IsSignal", acronym), &mc_v0.is_signal);
    Utils::CreateBranch(t, std::format("MC_{}_ReactionID", acronym), &mc_v0.reaction_id);
    Utils::CreateBranch(t, std::format("MC_{}_Mother_McEntry", acronym), &mc_v0.mother_mc_entry);
    Utils::CreateBranch(t, std::format("MC_{}_Mother_PdgCode", acronym), &mc_v0.mother_pdg_code);
    CreateBranches(t, mc_v0.neg, false, std::format("{}_Neg", acronym));
    CreateBranches(t, mc_v0.pos, false, std::format("{}_Pos", acronym));
}

inline void ClearBranches(Schema::Vector::MC_V0s& mc_v0) {
    mc_v0.mc_entry->clear();
    mc_v0.pdg_code->clear();
    ClearBranches(mc_v0.origin);
    ClearBranches(mc_v0.decay);
    ClearBranches(mc_v0.lv);
    mc_v0.is_true->clear();
    mc_v0.is_secondary->clear();
    mc_v0.is_signal->clear();
    mc_v0.is_hybrid->clear();
    mc_v0.reaction_id->clear();
    mc_v0.mother_mc_entry->clear();
    mc_v0.mother_pdg_code->clear();
    ClearBranches(mc_v0.neg, false);
    ClearBranches(mc_v0.pos, false);
}

}  // namespace Tree2Secondaries::IO
