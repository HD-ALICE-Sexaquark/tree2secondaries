#pragma once

#include <vector>

#include "App/Utilities.hxx"
#include "Schema/SchemaVectorPrimitives.hxx"

class TTree;

namespace Tree2Secondaries::Schema {

// NOTE:
// - read by `Packager` (with `include_sv=false`).
// - written by `Packager` and read by `Finder` (with `include_sv=true`).
struct VecInjected {

    // ## Storage ## //

    void ReadBranches(TTree* t, bool include_sv) {
        Utils::ReadBranch(t, "ReactionID", &reaction_id);
        mom.ReadBranches(t, "Sexa");
        if (include_sv) sv.ReadBranches(t, "SV");
        mom_nucleon.ReadBranches(t, "Nucleon");
    }

    void CreateBranches(TTree* t, bool include_sv) {
        Utils::CreateBranch(t, "ReactionID", &reaction_id);
        mom.CreateBranches(t, "Sexa");
        if (include_sv) sv.CreateBranches(t, "SV");
        mom_nucleon.CreateBranches(t, "Nucleon");
    }

    void ClearBranches(bool include_sv) {
        reaction_id->clear();
        mom.ClearBranches();
        if (include_sv) sv.ClearBranches();
        mom_nucleon.ClearBranches();
    }

    // ## Member Variables ## //

    Schema::VecMom3 mom{};
    Schema::VecCoordinates sv{};
    Schema::VecMom3 mom_nucleon{};
    std::vector<unsigned int>* reaction_id{nullptr};
};

}  // namespace Tree2Secondaries::Schema
