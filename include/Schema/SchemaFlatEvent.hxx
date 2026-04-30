#pragma once

#include "Math/Constants.hxx"
#include "Schema/SchemaFlatPrimitives.hxx"

class TTree;

namespace Tree2Secondaries::Schema {

// NOTE: read and written by `Packager` and `Finder`.
struct alignas(T2S_SIMD_ALIGN) FlatEvent {

    // ## Getters ## //

    [[nodiscard]] unsigned int RunNumber() const { return run_number; }
    [[nodiscard]] unsigned int DirNumber() const { return dir_number; }
    [[nodiscard]] unsigned int DirNumberB() const { return dir_number_b; }
    [[nodiscard]] unsigned int EventNumber() const { return event_number; }
    [[nodiscard]] float Centrality() const { return centrality; }
    [[nodiscard]] float MagneticField() const { return magnetic_field; }

    // ## Storage ## //

    // NOTE: used in `Packager` and `Finder`.
    void ReadBranches(TTree* t, bool is_mc) {
        Utils::ReadBranch(t, "RunNumber", &run_number);
        Utils::ReadBranch(t, "DirNumber", &dir_number);
        if (!is_mc) Utils::ReadBranch(t, "DirNumberB", &dir_number_b);
        Utils::ReadBranch(t, "EventNumber", &event_number);
        Utils::ReadBranch(t, "Centrality", &centrality);
        Utils::ReadBranch(t, "MagneticField", &magnetic_field);
        pv.ReadBranches(t, "PV");
        if (is_mc) mc_pv.ReadBranches(t, "MC_PV");
    }

    // NOTE: used in `Packager`.
    void CreateBranches(TTree* t, bool is_mc) {
        Utils::CreateBranch(t, "RunNumber", &run_number);
        Utils::CreateBranch(t, "DirNumber", &dir_number);
        if (!is_mc) Utils::CreateBranch(t, "DirNumberB", &dir_number_b);
        Utils::CreateBranch(t, "EventNumber", &event_number);
        Utils::CreateBranch(t, "Centrality", &centrality);
        Utils::CreateBranch(t, "MagneticField", &magnetic_field);
        pv.CreateBranches(t, "PV");
        if (is_mc) mc_pv.CreateBranches(t, "MC_PV");
    }

    // ## Member Variables ## //

    Schema::Coordinates pv{};
    Schema::Coordinates mc_pv{};  // NOTE: only read when analyzing MC
    unsigned int run_number{0};
    unsigned int dir_number{0};
    unsigned int dir_number_b{0};  // NOTE: only read when analyzing RD
    unsigned int event_number{0};
    float centrality{Const::DummyFloat};
    float magnetic_field{Const::DummyFloat};
};

}  // namespace Tree2Secondaries::Schema
