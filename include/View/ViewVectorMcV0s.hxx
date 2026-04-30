#pragma once

#include "Schema/SchemaVectorMcV0s.hxx"
#include "View/ViewVectorMcTracks.hxx"

namespace Tree2Secondaries::View {

struct VecMcV0s {

    VecMcV0s() = delete;
    VecMcV0s(const Schema::VecMcV0s* df, std::size_t idx) : source{df}, index{idx} {}

    // ## Getters ## //

    [[nodiscard]] int McEntry() const { return (*source->mc_entry)[index]; }
    [[nodiscard]] int PdgCode() const { return (*source->pdg_code)[index]; }
    [[nodiscard]] bool IsTrue() const { return static_cast<bool>((*source->is_true)[index]); }
    [[nodiscard]] bool IsSignal() const { return static_cast<bool>((*source->is_signal)[index]); }
    [[nodiscard]] bool IsSecondary() const { return static_cast<bool>((*source->is_secondary)[index]); }
    [[nodiscard]] bool IsHybrid() const { return static_cast<bool>((*source->is_hybrid)[index]); }
    [[nodiscard]] int ReactionID() const { return (*source->reaction_id)[index]; }

    [[nodiscard]] float Origin_X() const { return (*source->origin.x)[index]; }
    [[nodiscard]] float Origin_Y() const { return (*source->origin.y)[index]; }
    [[nodiscard]] float Origin_Z() const { return (*source->origin.z)[index]; }
    [[nodiscard]] float Decay_X() const { return (*source->decay.x)[index]; }
    [[nodiscard]] float Decay_Y() const { return (*source->decay.y)[index]; }
    [[nodiscard]] float Decay_Z() const { return (*source->decay.z)[index]; }
    [[nodiscard]] float Px() const { return (*source->lv.px)[index]; }
    [[nodiscard]] float Py() const { return (*source->lv.py)[index]; }
    [[nodiscard]] float Pz() const { return (*source->lv.pz)[index]; }
    [[nodiscard]] float Energy() const { return (*source->lv.energy)[index]; }

    [[nodiscard]] int Mother_McEntry() const { return (*source->mother_mc_entry)[index]; }
    [[nodiscard]] int Mother_PdgCode() const { return (*source->mother_pdg_code)[index]; }

    // ## Member Variables ## //

    const Schema::VecMcV0s* source{};
    std::size_t index{};
};

}  // namespace Tree2Secondaries::View
