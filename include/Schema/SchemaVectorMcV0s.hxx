#pragma once

#include <format>
#include <vector>

#include "App/Utilities.hxx"
#include "Math/Constants.hxx"
#include "Schema/SchemaVectorMcTracks.hxx"
#include "Schema/SchemaVectorPrimitives.hxx"

class TTree;

namespace Tree2Secondaries::Schema {

struct VecMcV0s {

    // ## Storage ## //

    void ReadBranches(TTree* t, std::string_view acronym) {
        Utils::ReadBranch(t, std::format("MC_{}_McEntry", acronym), &mc_entry);
        Utils::ReadBranch(t, std::format("MC_{}_PdgCode", acronym), &pdg_code);
        origin.ReadBranches(t, std::format("MC_{}_Origin", acronym));
        decay.ReadBranches(t, std::format("MC_{}_Decay", acronym));
        lv.ReadBranches(t, std::format("MC_{}", acronym));
        Utils::ReadBranch(t, std::format("MC_{}_IsTrue", acronym), &is_true);
        Utils::ReadBranch(t, std::format("MC_{}_IsSecondary", acronym), &is_secondary);
        Utils::ReadBranch(t, std::format("MC_{}_IsSignal", acronym), &is_signal);
        Utils::ReadBranch(t, std::format("MC_{}_IsHybrid", acronym), &is_hybrid);
        Utils::ReadBranch(t, std::format("MC_{}_ReactionID", acronym), &reaction_id);
        Utils::ReadBranch(t, std::format("MC_{}_Mother_McEntry", acronym), &mother_mc_entry);
        Utils::ReadBranch(t, std::format("MC_{}_Mother_PdgCode", acronym), &mother_pdg_code);
    }

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

    void PushDummy() {
        lv.PushDummy();
        origin.PushDummy();
        decay.PushDummy();
        mc_entry->push_back(Const::DummyInt);
        pdg_code->push_back(Const::DummyInt);
        is_true->push_back(static_cast<char>(false));
        is_secondary->push_back(static_cast<char>(false));
        is_signal->push_back(static_cast<char>(false));
        is_hybrid->push_back(static_cast<char>(false));
        reaction_id->push_back(Const::DummyInt);
        mother_mc_entry->push_back(Const::DummyInt);
        mother_pdg_code->push_back(Const::DummyInt);
    }

    void ClearBranches() {
        lv.ClearBranches();
        origin.ClearBranches();
        decay.ClearBranches();
        mc_entry->clear();
        pdg_code->clear();
        is_true->clear();
        is_secondary->clear();
        is_signal->clear();
        is_hybrid->clear();
        reaction_id->clear();
        mother_mc_entry->clear();
        mother_pdg_code->clear();
    }

    // ## Member Variables ## //

    Schema::VecMom4 lv{};
    Schema::VecCoordinates origin{};
    Schema::VecCoordinates decay{};
    std::vector<int>* mc_entry{nullptr};
    std::vector<int>* pdg_code{nullptr};
    std::vector<char>* is_true{nullptr};
    std::vector<char>* is_secondary{nullptr};
    std::vector<char>* is_signal{nullptr};
    std::vector<char>* is_hybrid{nullptr};
    std::vector<int>* reaction_id{nullptr};
    std::vector<int>* mother_mc_entry{nullptr};
    std::vector<int>* mother_pdg_code{nullptr};
};

}  // namespace Tree2Secondaries::Schema
