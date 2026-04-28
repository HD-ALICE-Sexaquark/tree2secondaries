#pragma once

#include "Storage/Schema/SchemaVector.hxx"
#include "View/BaseView.hxx"
#include "View/MC/ViewMcPackedTrack.hxx"

namespace Tree2Secondaries::View::MC {

struct PackedV0 : View::Base<Schema::Vector::MC_V0s, int> {

    PackedV0() = delete;
    PackedV0(const Schema::Vector::MC_V0s* df, int entry)  //
        : View::Base<Schema::Vector::MC_V0s, int>{df, entry} {}

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

    [[nodiscard]] float Decay_X() const { return (*Source->decay.x)[Entry]; }
    [[nodiscard]] float Decay_Y() const { return (*Source->decay.y)[Entry]; }
    [[nodiscard]] float Decay_Z() const { return (*Source->decay.z)[Entry]; }

    [[nodiscard]] int Neg_Entry() const { return (*Source->neg.mc_entry)[Entry]; }
    [[nodiscard]] int Neg_PdgCode() const { return (*Source->neg.pdg_code)[Entry]; }
    [[nodiscard]] int Neg_ReactionID() const { return (*Source->neg.reaction_id)[Entry]; }
    [[nodiscard]] bool Neg_IsTrue() const { return static_cast<bool>((*Source->neg.is_true)[Entry]); }
    [[nodiscard]] bool Neg_IsSignal() const { return static_cast<bool>((*Source->neg.is_signal)[Entry]); }
    [[nodiscard]] bool Neg_IsSecondary() const { return static_cast<bool>((*Source->neg.is_secondary)[Entry]); }

    [[nodiscard]] int Pos_Entry() const { return (*Source->pos.mc_entry)[Entry]; }
    [[nodiscard]] int Pos_PdgCode() const { return (*Source->pos.pdg_code)[Entry]; }
    [[nodiscard]] int Pos_ReactionID() const { return (*Source->pos.reaction_id)[Entry]; }
    [[nodiscard]] bool Pos_IsTrue() const { return static_cast<bool>((*Source->pos.is_true)[Entry]); }
    [[nodiscard]] bool Pos_IsSignal() const { return static_cast<bool>((*Source->pos.is_signal)[Entry]); }
    [[nodiscard]] bool Pos_IsSecondary() const { return static_cast<bool>((*Source->pos.is_secondary)[Entry]); }

    [[nodiscard]] int Mother_McEntry() const { return (*Source->mother_mc_entry)[Entry]; }
    [[nodiscard]] int Mother_PdgCode() const { return (*Source->mother_pdg_code)[Entry]; }

    [[nodiscard]] float Neg_Px() const { return (*Source->neg.lv.px)[Entry]; }
    [[nodiscard]] float Neg_Py() const { return (*Source->neg.lv.py)[Entry]; }
    [[nodiscard]] float Neg_Pz() const { return (*Source->neg.lv.pz)[Entry]; }

    [[nodiscard]] float Pos_Px() const { return (*Source->pos.lv.px)[Entry]; }
    [[nodiscard]] float Pos_Py() const { return (*Source->pos.lv.py)[Entry]; }
    [[nodiscard]] float Pos_Pz() const { return (*Source->pos.lv.pz)[Entry]; }

    [[nodiscard]] bool IsHybrid() const { return static_cast<bool>((*Source->is_hybrid)[Entry]); }
};

}  // namespace Tree2Secondaries::View::MC
