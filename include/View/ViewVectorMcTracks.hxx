#pragma once

#include <cstddef>

#include "Schema/SchemaVectorMcTracks.hxx"

namespace Tree2Secondaries::View {

struct VecMcTracks {

    VecMcTracks() = delete;
    VecMcTracks(const Schema::VecMcTracks* df, std::size_t idx) : source{df}, index{idx} {}

    // ## Getters ## //

    [[nodiscard]] int McEntry() const { return (*source->mc_entry)[index]; }
    [[nodiscard]] int PdgCode() const { return (*source->pdg_code)[index]; }
    [[nodiscard]] int ReactionID() const { return (*source->reaction_id)[index]; }
    [[nodiscard]] bool IsTrue() const { return static_cast<bool>((*source->is_true)[index]); }
    [[nodiscard]] bool IsSecondary() const { return static_cast<bool>((*source->is_secondary)[index]); }
    [[nodiscard]] bool IsSignal() const { return static_cast<bool>((*source->is_signal)[index]); }

    [[nodiscard]] float Origin_X() const { return (*source->origin.x)[index]; }
    [[nodiscard]] float Origin_Y() const { return (*source->origin.y)[index]; }
    [[nodiscard]] float Origin_Z() const { return (*source->origin.z)[index]; }
    [[nodiscard]] float Px() const { return (*source->lv.px)[index]; }
    [[nodiscard]] float Py() const { return (*source->lv.py)[index]; }
    [[nodiscard]] float Pz() const { return (*source->lv.pz)[index]; }
    [[nodiscard]] float Energy() const { return (*source->lv.energy)[index]; }

    [[nodiscard]] int Mother_McEntry() const { return (*source->mother_mc_entry)[index]; }
    [[nodiscard]] int Mother_PdgCode() const { return (*source->mother_pdg_code)[index]; }

    [[nodiscard]] int GrandMother_McEntry() const { return (*source->gm_mc_entry)[index]; }
    [[nodiscard]] int GrandMother_PdgCode() const { return (*source->gm_pdg_code)[index]; }

    // ## Member Variables ## //

    const Schema::VecMcTracks* source{};
    std::size_t index{};
};

}  // namespace Tree2Secondaries::View
