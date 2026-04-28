#pragma once

#include <vector>

#include "Storage/BaseStorage.hxx"

namespace Tree2Secondaries::Schema::Vector {

// NOTE: read by `Packager`.
struct MCParticles {
    Storage::VecMom4 lv{};
    Storage::VecCoordinates origin{};
    std::vector<int>* pdg_code{nullptr};
    std::vector<char>* charge{nullptr};
    std::vector<int>* mother_mc_entry{nullptr};
    std::vector<int>* n_daughters{nullptr};
    std::vector<int>* firstdau_mc_entry{nullptr};
    std::vector<int>* lastdau_mc_entry{nullptr};
    std::vector<int>* status{nullptr};
    std::vector<int>* generator{nullptr};
    std::vector<char>* is_phys_primary{nullptr};
    std::vector<char>* is_sec_from_mat{nullptr};
    std::vector<char>* is_sec_from_weak{nullptr};
};

// NOTE:
// - read by `Packager` (with `include_sv=false`).
// - written by `Packager` and read by `Finder` (with `include_sv=true`).
struct Injected {
    Storage::VecMom3 mom{};
    Storage::VecCoordinates sv{};
    Storage::VecMom3 mom_nucleon{};
    std::vector<int>* reaction_id{nullptr};
};

// NOTE:
// - read by `Packager` (with `include_mc_entry=IsMC()`).
// - written by `Packager` and read by `Finder` (with `include_mc_entry=false`).
struct Tracks {
    Storage::VecState6 state{};
    Storage::VecCovMatrix<6> cov{};
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

// NOTE:
// - read by `Packager` (with `include_mc_entry=IsMC()`).
// - written by `Packager` and read by `Finder` (with `include_mc_entry=false`).
struct MC_Tracks {
    Storage::VecMom4 lv{};
    Storage::VecCoordinates origin{};
    std::vector<int>* mc_entry{nullptr};
    std::vector<int>* pdg_code{nullptr};
    std::vector<int>* reaction_id{nullptr};
    std::vector<char>* is_true{nullptr};
    std::vector<char>* is_secondary{nullptr};
    std::vector<char>* is_signal{nullptr};
    std::vector<int>* mother_mc_entry{nullptr};
    std::vector<int>* mother_pdg_code{nullptr};
    std::vector<int>* gm_mc_entry{nullptr};
    std::vector<int>* gm_pdg_code{nullptr};
};

struct V0s {
    Vector::Tracks neg{};
    Vector::Tracks pos{};
    Storage::VecState6 neg_pca_v0{};
    Storage::VecState6 pos_pca_v0{};
    Storage::VecMom4 lv{};
    Storage::VecCoordinates decay{};
    Storage::VecCovMatrix<7> cov{};
    std::vector<float>* chi2ndf{nullptr};
};

struct MC_V0s {
    Vector::MC_Tracks neg{};
    Vector::MC_Tracks pos{};
    Storage::VecMom4 lv{};
    Storage::VecCoordinates origin{};
    Storage::VecCoordinates decay{};
    std::vector<int>* mc_entry{nullptr};
    std::vector<int>* pdg_code{nullptr};
    std::vector<int>* reaction_id{nullptr};
    std::vector<char>* is_true{nullptr};
    std::vector<char>* is_signal{nullptr};
    std::vector<char>* is_secondary{nullptr};
    std::vector<char>* is_hybrid{nullptr};
    std::vector<int>* mother_mc_entry{nullptr};
    std::vector<int>* mother_pdg_code{nullptr};
};

}  // namespace Tree2Secondaries::Schema::Vector
