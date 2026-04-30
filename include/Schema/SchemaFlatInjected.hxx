#pragma once

#include "Math/Constants.hxx"
#include "Schema/SchemaFlatPrimitives.hxx"

class TTree;

namespace Tree2Secondaries::Schema {

// NOTE: written by `Finder` to another TTree.
struct alignas(T2S_SIMD_ALIGN) FlatInjected {

    // ## Storage ## //

    void CreateBranches(TTree* t) {
        Utils::CreateBranch(t, "RunNumber", &run_number);
        Utils::CreateBranch(t, "DirNumber", &dir_number);
        Utils::CreateBranch(t, "EventNumber", &event_number);
        Utils::CreateBranch(t, "ReactionID", &reaction_id);
        sv.CreateBranches(t, "SV");
        lv.CreateBranches(t, "Sexa");
        lv_nucleon.CreateBranches(t, "Nucleon");
    }

    // ## Member Variables ## //

    Schema::Mom4 lv{};
    Schema::Mom4 lv_nucleon{};
    Schema::Coordinates sv{};
    unsigned int run_number{0};
    unsigned int dir_number{0};
    unsigned int event_number{0};
    unsigned int reaction_id{0};
};

}  // namespace Tree2Secondaries::Schema
