#pragma once

#include <vector>

#include "App/Utilities.hxx"
#include "Schema/SchemaVectorPrimitives.hxx"

class TTree;

namespace Tree2Secondaries::Schema {

// Handles MC particle info, with some cached info from related particles (mother and daughters).
// NOTE: only read by `Packager`.
struct VecMcParticles {

    // ## Storage ## //

    void ReadBranches(TTree* t) {
        Utils::ReadBranch(t, "MC_PdgCode", &pdg_code);
        Utils::ReadBranch(t, "MC_Charge", &charge);
        origin.ReadBranches(t, "MC_Origin");
        lv.ReadBranches(t, "MC");
        Utils::ReadBranch(t, "MC_Mother_McEntry", &mother_mc_entry);
        Utils::ReadBranch(t, "MC_N_Daughters", &n_daughters);
        Utils::ReadBranch(t, "MC_FirstDau_McEntry", &firstdau_mc_entry);
        Utils::ReadBranch(t, "MC_LastDau_McEntry", &lastdau_mc_entry);
        Utils::ReadBranch(t, "MC_Status", &status);
        Utils::ReadBranch(t, "MC_Generator", &generator);
        Utils::ReadBranch(t, "MC_IsPhysPrimary", &is_phys_primary);
        Utils::ReadBranch(t, "MC_IsSecFromMat", &is_sec_from_mat);
        Utils::ReadBranch(t, "MC_IsSecFromWeak", &is_sec_from_weak);
    }

    // ## Member Variables ## //

    Schema::VecMom4 lv{};
    Schema::VecCoordinates origin{};
    std::vector<int>* pdg_code{nullptr};
    std::vector<char>* charge{nullptr};
    std::vector<int>* mother_mc_entry{nullptr};
    std::vector<unsigned int>* n_daughters{nullptr};
    std::vector<int>* firstdau_mc_entry{nullptr};
    std::vector<int>* lastdau_mc_entry{nullptr};
    std::vector<unsigned int>* status{nullptr};
    std::vector<char>* generator{nullptr};
    std::vector<char>* is_phys_primary{nullptr};
    std::vector<char>* is_sec_from_mat{nullptr};
    std::vector<char>* is_sec_from_weak{nullptr};
};

}  // namespace Tree2Secondaries::Schema
