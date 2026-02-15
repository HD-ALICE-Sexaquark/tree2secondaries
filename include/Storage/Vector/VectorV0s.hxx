#pragma once

#include "Storage/Vector/BaseVector.hxx"
#include "Storage/Vector/VectorTracks.hxx"

namespace Tree2Secondaries::Storage::Vector {

// `Vector::States` + `Vector::CovMatrices` +
// `Neg` (`Vector::Tracks`) + `Pos` (`Vector::Tracks`) + `Neg_atPCA` (`States_NoE`) +
// `Pos_atPCA` (`States_NoE`) + `Chi2NDF`.
struct V0s : Vector::States, Vector::CovMatrices {
    Vector::Tracks Neg{};
    Vector::Tracks Pos{};
    Vector::States_NoE Neg_atPCA{};
    Vector::States_NoE Pos_atPCA{};
    std::vector<float>* Chi2NDF{nullptr};

    void Clear_VectorV0s() {
        Clear_VectorStates();
        Clear_VectorCovMatrices();
        Neg.Clear_VectorTracks(false, true);
        Pos.Clear_VectorTracks(false, true);
        Neg_atPCA.Clear_VectorStates_NoE();
        Pos_atPCA.Clear_VectorStates_NoE();
        Chi2NDF->clear();
    }
    void CreateBranches_VectorV0s(TTree* tree, std::string_view prefix) {
        const std::string neg_prefix = std::format("{}_Neg", prefix);
        const std::string pos_prefix = std::format("{}_Pos", prefix);
        CreateBranches_VectorStates(tree, prefix);
        CreateBranches_VectorCovMatrices(tree, prefix);
        Neg.CreateBranches_VectorTracks(tree, false, true, neg_prefix);
        Pos.CreateBranches_VectorTracks(tree, false, true, pos_prefix);
        Neg_atPCA.CreateBranches_VectorStates_NoE(tree, std::format("{}_atPCA", neg_prefix));
        Pos_atPCA.CreateBranches_VectorStates_NoE(tree, std::format("{}_atPCA", pos_prefix));
        tree->Branch(std::format("{}_Chi2NDF", prefix).c_str(), &Chi2NDF);
    }
    void ReadBranches_VectorV0s(TTree* tree, std::string_view prefix) {
        const std::string neg_prefix = std::format("{}_Neg", prefix);
        const std::string pos_prefix = std::format("{}_Pos", prefix);
        ReadBranches_VectorStates(tree, prefix);
        ReadBranches_VectorCovMatrices(tree, prefix);
        Neg.ReadBranches_VectorTracks(tree, false, true, neg_prefix);
        Pos.ReadBranches_VectorTracks(tree, false, true, pos_prefix);
        Neg_atPCA.ReadBranches_VectorStates_NoE(tree, std::format("{}_atPCA", neg_prefix));
        Pos_atPCA.ReadBranches_VectorStates_NoE(tree, std::format("{}_atPCA", pos_prefix));
        Utils::ReadBranch(tree, std::format("{}_Chi2NDF", prefix), &Chi2NDF);
    }
};

// `Vector::MC` + `Vector::States` +
// `Neg` (`Vector::MC`) + `Pos` (`Vector::MC`) +
// `Mother` (`Vector::MC_Id`) +
// `Neg_Momentum` (`Vector::PxPyPz`) + `Pos_Momentum` (`Vector::PxPyPz`) +
// `AtDecay` (`Vector::Coordinates`) +
// `IsHybrid`.
struct MC_V0s : Vector::MC, Vector::States {
    Vector::MC Neg{};
    Vector::MC Pos{};
    Vector::MC_Id Mother{};
    Vector::PxPyPz Neg_Momentum{};
    Vector::PxPyPz Pos_Momentum{};
    Vector::Coordinates AtDecay{};
    std::vector<char>* IsHybrid{nullptr};  // NOTE: `char` instead of `bool`, because ROOT

    void Clear_VectorMC_V0s() {
        Clear_VectorMC();
        Clear_VectorStates();
        Neg.Clear_VectorMC();
        Pos.Clear_VectorMC();
        Mother.Clear_VectorMC_Id();
        Neg_Momentum.Clear_VectorPxPyPz();
        Pos_Momentum.Clear_VectorPxPyPz();
        AtDecay.Clear_VectorCoordinates();
        IsHybrid->clear();
    }
    void CreateBranches_VectorMC_V0s(TTree* tree, std::string_view acronym = "") {
        CreateBranches_VectorMC(tree, acronym);
        CreateBranches_VectorStates(tree, std::format("MC_{}", acronym));
        Neg.CreateBranches_VectorMC(tree, std::format("{}_Neg", acronym));
        Pos.CreateBranches_VectorMC(tree, std::format("{}_Pos", acronym));
        Mother.CreateBranches_VectorMC_Id(tree, std::format("{}_Mother", acronym));
        Neg_Momentum.CreateBranches_VectorPxPyPz(tree, std::format("MC_{}_Neg", acronym));
        Pos_Momentum.CreateBranches_VectorPxPyPz(tree, std::format("MC_{}_Pos", acronym));
        AtDecay.CreateBranches_VectorCoordinates(tree, std::format("MC_{}_atDecay", acronym), "");
        tree->Branch(std::format("MC_{}_IsHybrid", acronym).c_str(), &IsHybrid);
    }
    void ReadBranches_VectorMC_V0s(TTree* tree, std::string_view acronym = "") {
        ReadBranches_VectorMC(tree, acronym);
        ReadBranches_VectorStates(tree, std::format("MC_{}", acronym));
        Neg.ReadBranches_VectorMC(tree, std::format("{}_Neg", acronym));
        Pos.ReadBranches_VectorMC(tree, std::format("{}_Pos", acronym));
        Mother.ReadBranches_VectorMC_Id(tree, std::format("{}_Mother", acronym));
        Neg_Momentum.ReadBranches_VectorPxPyPz(tree, std::format("MC_{}_Neg", acronym));
        Pos_Momentum.ReadBranches_VectorPxPyPz(tree, std::format("MC_{}_Pos", acronym));
        AtDecay.ReadBranches_VectorCoordinates(tree, std::format("MC_{}_atDecay", acronym), "");
        Utils::ReadBranch(tree, std::format("MC_{}_IsHybrid", acronym), &IsHybrid);
    }
};

}  // namespace Tree2Secondaries::Storage::Vector
