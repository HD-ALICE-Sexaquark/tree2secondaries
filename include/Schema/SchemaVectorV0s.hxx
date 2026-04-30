#pragma once

#include <format>
#include <vector>

#include "App/Utilities.hxx"
#include "Schema/SchemaVectorPrimitives.hxx"
#include "Schema/SchemaVectorTracks.hxx"

class TTree;

namespace Tree2Secondaries::Schema {

struct VecV0s {

    // ## Storage ## //

    void ReadBranches(TTree* t, std::string_view acronym) {
        decay.ReadBranches(t, std::format("{}_Decay", acronym));
        lv.ReadBranches(t, acronym);
        Utils::ReadBranch(t, std::format("{}_Chi2NDF", acronym), &chi2ndf);
        cov.ReadBranches(t, acronym);
        neg.ReadBranches(t, false, std::format("{}_Neg", acronym));
        neg_pca_v0.ReadBranches(t, std::format("{}_Neg_PCAwrtV0", acronym));
        pos.ReadBranches(t, false, std::format("{}_Pos", acronym));
        pos_pca_v0.ReadBranches(t, std::format("{}_Pos_PCAwrtV0", acronym));
    }

    void CreateBranches(TTree* t, std::string_view acronym) {
        decay.CreateBranches(t, std::format("{}_Decay", acronym));
        lv.CreateBranches(t, acronym);
        Utils::CreateBranch(t, std::format("{}_Chi2NDF", acronym), &chi2ndf);
        cov.CreateBranches(t, acronym);
        neg.CreateBranches(t, false, std::format("{}_Neg", acronym));
        neg_pca_v0.CreateBranches(t, std::format("{}_Neg_PCAwrtV0", acronym));
        pos.CreateBranches(t, false, std::format("{}_Pos", acronym));
        pos_pca_v0.CreateBranches(t, std::format("{}_Pos_PCAwrtV0", acronym));
    }

    void ClearBranches() {
        decay.ClearBranches();
        lv.ClearBranches();
        chi2ndf->clear();
        cov.ClearBranches();
        neg.ClearBranches(false);
        neg_pca_v0.ClearBranches();
        pos.ClearBranches(false);
        pos_pca_v0.ClearBranches();
    }

    // ## Member Variables ## //

    Schema::VecTracks neg{};
    Schema::VecTracks pos{};
    Schema::VecState6 neg_pca_v0{};
    Schema::VecState6 pos_pca_v0{};
    Schema::VecMom4 lv{};
    Schema::VecCoordinates decay{};
    Schema::VecCovMatrix<7> cov{};
    std::vector<float>* chi2ndf{nullptr};
};

}  // namespace Tree2Secondaries::Schema
