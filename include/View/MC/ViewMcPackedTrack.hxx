#pragma once

#include "Storage/Schema/SchemaVector.hxx"
#include "View/BaseView.hxx"

namespace Tree2Secondaries::View::MC {

// Slight transformations to `View::MC::Particle`.
struct PackedTrack : View::Base<Schema::Vector::MC_Tracks, int> {

    PackedTrack() = delete;
    PackedTrack(const Schema::Vector::MC_Tracks* df, int entry)  //
        : View::Base<Schema::Vector::MC_Tracks, int>{df, entry} {}

    [[nodiscard]] int McEntry() const { return (*Source->mc_entry)[Entry]; }
    [[nodiscard]] int PdgCode() const { return (*Source->pdg_code)[Entry]; }
    [[nodiscard]] int ReactionID() const { return (*Source->reaction_id)[Entry]; }
    [[nodiscard]] bool IsTrue() const { return static_cast<bool>((*Source->is_true)[Entry]); }
    [[nodiscard]] bool IsSignal() const { return static_cast<bool>((*Source->is_signal)[Entry]); }
    [[nodiscard]] bool IsSecondary() const { return static_cast<bool>((*Source->is_secondary)[Entry]); }

    [[nodiscard]] float Origin_X() const { return (*Source->origin.x)[Entry]; }
    [[nodiscard]] float Origin_Y() const { return (*Source->origin.y)[Entry]; }
    [[nodiscard]] float Origin_Z() const { return (*Source->origin.z)[Entry]; }
    [[nodiscard]] float Px() const { return (*Source->lv.px)[Entry]; }
    [[nodiscard]] float Py() const { return (*Source->lv.py)[Entry]; }
    [[nodiscard]] float Pz() const { return (*Source->lv.pz)[Entry]; }
    [[nodiscard]] float Energy() const { return (*Source->lv.energy)[Entry]; }

    [[nodiscard]] int Mother_McEntry() const { return (*Source->mother_mc_entry)[Entry]; }
    [[nodiscard]] int Mother_PdgCode() const { return (*Source->mother_pdg_code)[Entry]; }

    [[nodiscard]] int GrandMother_McEntry() const { return (*Source->gm_mc_entry)[Entry]; }
    [[nodiscard]] int GrandMother_PdgCode() const { return (*Source->gm_pdg_code)[Entry]; }
};

}  // namespace Tree2Secondaries::View::MC
