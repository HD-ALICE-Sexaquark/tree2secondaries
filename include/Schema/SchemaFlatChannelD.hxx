#pragma once

#include "Math/Constants.hxx"
#include "Schema/SchemaFlatEvent.hxx"
#include "Schema/SchemaFlatPrimitives.hxx"
#include "Schema/SchemaFlatTrack.hxx"
#include "Schema/SchemaFlatV0.hxx"

class TTree;

namespace Tree2Secondaries::Schema {

struct alignas(T2S_SIMD_ALIGN) FlatChannelD {

    // ## Storage ## //

    void CreateBranches(TTree* t, bool is_mc) {
        event.CreateBranches(t, is_mc);
        sv.CreateBranches(t, "SV");
        lv.CreateBranches(t, "Sexa");
        Utils::CreateBranch(t, "Sexa_Chi2NDF", &chi2ndf);
        Utils::CreateBranch(t, "Sexa_E_MinusNucleon", &e_minus_nucleon);
        Utils::CreateBranch(t, "Sexa_IsAntiChannel", &is_anti);
        v0.CreateBranches(t, "V0");
        kaon.CreateBranches(t, "Kaon");
        v0_pca_sv.CreateBranches(t, "V0_PCAwrtSV");
        kaon_pca_sv.CreateBranches(t, "Kaon_PCAwrtSV");
    }

    // ## Member Variables ## //

    Schema::FlatEvent event{};
    Schema::Coordinates sv{};
    Schema::Mom4 lv{};
    float chi2ndf{Const::DummyFloat};
    float e_minus_nucleon{Const::DummyFloat};
    bool is_anti{false};
    Schema::FlatTrack kaon{};
    Schema::FlatV0 v0{};
    Schema::State6 kaon_pca_sv{};
    Schema::State6 v0_pca_sv{};
};

struct alignas(T2S_SIMD_ALIGN) FlatMcChannelD {

    // ## Storage ## //

    void CreateBranches(TTree* t) {
        before.CreateBranches(t, "MC_Before");
        after.CreateBranches(t, "MC_After");
        nucleon.CreateBranches(t, "MC_Nucleon");
        sv.CreateBranches(t, "MC_SV");
        Utils::CreateBranch(t, "MC_ReactionID", &reaction_id);
        Utils::CreateBranch(t, "MC_IsSignal", &is_signal);
        Utils::CreateBranch(t, "MC_IsHybrid", &is_hybrid);
        v0.CreateBranches(t, "V0");
        kaon.CreateBranches(t, true, "Kaon");
    }

    // ## Member Variables ## //

    Schema::Mom4 before{};
    Schema::Mom4 after{};
    Schema::Mom4 nucleon{};
    Schema::Coordinates sv{};
    int reaction_id{Const::DummyInt};  // NOTE: negative when not found
    bool is_signal{false};
    bool is_hybrid{false};
    Schema::FlatMcV0 v0{};
    Schema::FlatMcTrack kaon{};
};

}  // namespace Tree2Secondaries::Schema
