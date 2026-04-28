#pragma once

#include "Math/Constants.hxx"
#include "Storage/BaseStorage.hxx"

namespace Tree2Secondaries::Schema::Flat {

// Event //

// NOTE: read and written by `Packager` and `Finder`.
struct alignas(T2S_SIMD_ALIGN) Event {
    Storage::Coordinates pv{};
    Storage::Coordinates mc_pv{};  // NOTE: only read when analyzing MC
    unsigned int run_number{0};
    unsigned int dir_number{0};
    unsigned int dir_number_b{0};  // NOTE: only read when analyzing RD
    unsigned int event_number{0};
    float centrality{Const::DummyFloat};
    float magnetic_field{Const::DummyFloat};
};

// Injected //

// NOTE: written by `Finder` to another TTree.
struct alignas(T2S_SIMD_ALIGN) Injected {
    Storage::Mom4 lv{};
    Storage::Mom4 lv_nucleon{};
    Storage::Coordinates sv{};
    unsigned int run_number{0};
    unsigned int dir_number{0};
    unsigned int event_number{0};
    unsigned int reaction_id{0};
};

// Track //

struct alignas(T2S_SIMD_ALIGN) Track {
    Storage::State6 state{};
    unsigned int esd_entry{0};
    char charge{0};
    float pre_dca_xy{Const::DummyFloat};
    float pre_dca_z{Const::DummyFloat};
    float tpc_signal{Const::DummyFloat};
    float n_sigma_pion{Const::DummyFloat};
    float n_sigma_kaon{Const::DummyFloat};
    float n_sigma_proton{Const::DummyFloat};
};

struct alignas(T2S_SIMD_ALIGN) MC_Track {
    Storage::Mom4 lv{};
    Storage::Coordinates origin{};
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

// V0 //

struct alignas(T2S_SIMD_ALIGN) V0 {
    Flat::Track neg{};
    Flat::Track pos{};
    Storage::State6 neg_pca_v0{};
    Storage::State6 pos_pca_v0{};
    Storage::Mom4 lv{};
    Storage::Coordinates decay{};
    float chi2ndf{Const::DummyFloat};
};

struct alignas(T2S_SIMD_ALIGN) MC_V0 {
    Flat::MC_Track neg{};
    Flat::MC_Track pos{};
    Storage::Mom4 lv{};
    Storage::Coordinates origin{};
    Storage::Coordinates decay{};
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

// Channel A //

struct alignas(T2S_SIMD_ALIGN) ChannelA {
    Flat::Event event{};
    Storage::Coordinates sv{};
    Storage::Mom4 lv{};
    float chi2ndf{Const::DummyFloat};
    float e_minus_nucleon{Const::DummyFloat};
    bool is_anti{false};
    //
    Flat::V0 v0a{};
    Flat::V0 v0b{};
    Storage::State6 v0a_pca_sv{};
    Storage::State6 v0b_pca_sv{};
};

struct alignas(T2S_SIMD_ALIGN) MC_ChannelA {
    Storage::Mom4 before{};
    Storage::Mom4 after{};
    Storage::Mom4 nucleon{};
    Storage::Coordinates sv{};
    int reaction_id{Const::DummyInt};  // NOTE: negative when not found
    bool is_signal{false};
    bool is_hybrid{false};
    //
    Flat::MC_V0 v0a{};
    Flat::MC_V0 v0b{};
};

// Channel D //

struct alignas(T2S_SIMD_ALIGN) ChannelD {
    Flat::Event event{};
    Storage::Coordinates sv{};
    Storage::Mom4 lv{};
    float chi2ndf{Const::DummyFloat};
    float e_minus_nucleon{Const::DummyFloat};
    bool is_anti{false};
    //
    Flat::Track kaon{};
    Flat::V0 v0{};
    Storage::State6 kaon_pca_sv{};
    Storage::State6 v0_pca_sv{};
};

struct alignas(T2S_SIMD_ALIGN) MC_ChannelD {
    Storage::Mom4 before{};
    Storage::Mom4 after{};
    Storage::Mom4 nucleon{};
    Storage::Coordinates sv{};
    int reaction_id{Const::DummyInt};  // NOTE: negative when not found
    bool is_signal{false};
    bool is_hybrid{false};
    //
    Flat::MC_V0 v0{};
    Flat::MC_Track kaon{};
};

}  // namespace Tree2Secondaries::Schema::Flat
