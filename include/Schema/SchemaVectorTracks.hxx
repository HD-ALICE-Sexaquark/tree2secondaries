#pragma once

#include <format>
#include <vector>

#include "App/Utilities.hxx"
#include "Math/Constants.hxx"
#include "Schema/SchemaVectorPrimitives.hxx"

class TTree;

namespace Tree2Secondaries::Schema {

// NOTE:
// - read by `Packager` (with `include_mc_entry=IsMC()`).
// - written by `Packager` and read by `Finder` (with `include_mc_entry=false`).
struct VecTracks {

    // ## Storage ## //

    void ReadBranches(TTree* t, bool include_mc_entry, std::string_view prefix = "Track") {
        Utils::ReadBranch(t, std::format("{}_EsdEntry", prefix), &esd_entry);
        Utils::ReadBranch(t, std::format("{}_Charge", prefix), &charge);
        state.ReadBranches(t, prefix);
        Utils::ReadBranch(t, std::format("{}_PreDCAxy", prefix), &pre_dca_xy);
        Utils::ReadBranch(t, std::format("{}_PreDCAz", prefix), &pre_dca_z);
        Utils::ReadBranch(t, std::format("{}_TPC_Signal", prefix), &tpc_signal);
        Utils::ReadBranch(t, std::format("{}_NSigmaPion", prefix), &n_sigma_pion);
        Utils::ReadBranch(t, std::format("{}_NSigmaKaon", prefix), &n_sigma_kaon);
        Utils::ReadBranch(t, std::format("{}_NSigmaProton", prefix), &n_sigma_proton);
        cov.ReadBranches(t, prefix);
        if (include_mc_entry) Utils::ReadBranch(t, std::format("{}_McEntry", prefix), &mc_entry);
    }

    void CreateBranches(TTree* t, bool include_mc_entry, std::string_view prefix = "Track") {
        Utils::CreateBranch(t, std::format("{}_EsdEntry", prefix), &esd_entry);
        Utils::CreateBranch(t, std::format("{}_Charge", prefix), &charge);
        state.CreateBranches(t, prefix);
        Utils::CreateBranch(t, std::format("{}_PreDCAxy", prefix), &pre_dca_xy);
        Utils::CreateBranch(t, std::format("{}_PreDCAz", prefix), &pre_dca_z);
        Utils::CreateBranch(t, std::format("{}_TPC_Signal", prefix), &tpc_signal);
        Utils::CreateBranch(t, std::format("{}_NSigmaPion", prefix), &n_sigma_pion);
        Utils::CreateBranch(t, std::format("{}_NSigmaKaon", prefix), &n_sigma_kaon);
        Utils::CreateBranch(t, std::format("{}_NSigmaProton", prefix), &n_sigma_proton);
        cov.CreateBranches(t, prefix);
        if (include_mc_entry) Utils::CreateBranch(t, std::format("{}_McEntry", prefix), &mc_entry);
    }

    void ClearBranches(bool include_mc_entry) {
        state.ClearBranches();
        cov.ClearBranches();
        esd_entry->clear();
        charge->clear();
        pre_dca_xy->clear();
        pre_dca_z->clear();
        tpc_signal->clear();
        n_sigma_pion->clear();
        n_sigma_kaon->clear();
        n_sigma_proton->clear();
        if (include_mc_entry) mc_entry->clear();
    }

    // ## Member Variables ## //

    Schema::VecState6 state{};
    Schema::VecCovMatrix<6> cov{};
    std::vector<unsigned int>* esd_entry{nullptr};
    std::vector<char>* charge{nullptr};
    std::vector<float>* pre_dca_xy{nullptr};
    std::vector<float>* pre_dca_z{nullptr};
    std::vector<float>* tpc_signal{nullptr};
    std::vector<float>* n_sigma_pion{nullptr};
    std::vector<float>* n_sigma_kaon{nullptr};
    std::vector<float>* n_sigma_proton{nullptr};
    std::vector<int>* mc_entry{nullptr};  // MC only.
};

}  // namespace Tree2Secondaries::Schema
