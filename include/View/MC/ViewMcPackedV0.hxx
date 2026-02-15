#pragma once

#include <cstdlib>

#include "Storage/Vector/VectorV0s.hxx"
#include "View/BaseView.hxx"

namespace Tree2Secondaries::View::MC {

struct PackedV0 : View::Base<Storage::Vector::MC_V0s, size_t> {

    PackedV0() = delete;
    PackedV0(const Storage::Vector::MC_V0s* df, size_t entry)  //
        : View::Base<Storage::Vector::MC_V0s, size_t>{df, entry} {}

    [[nodiscard]] int McEntry() const { return (*Source->McEntry)[Entry]; }
    [[nodiscard]] int PdgCode() const { return (*Source->PdgCode)[Entry]; }
    [[nodiscard]] int ReactionID() const { return (*Source->ReactionID)[Entry]; }
    [[nodiscard]] bool IsTrue() const { return static_cast<bool>((*Source->IsTrue)[Entry]); }
    [[nodiscard]] bool IsSignal() const { return static_cast<bool>((*Source->IsSignal)[Entry]); }
    [[nodiscard]] bool IsSecondary() const { return static_cast<bool>((*Source->IsSecondary)[Entry]); }

    [[nodiscard]] float X() const { return (*Source->X)[Entry]; }
    [[nodiscard]] float Y() const { return (*Source->Y)[Entry]; }
    [[nodiscard]] float Z() const { return (*Source->Z)[Entry]; }
    [[nodiscard]] float Px() const { return (*Source->Px)[Entry]; }
    [[nodiscard]] float Py() const { return (*Source->Py)[Entry]; }
    [[nodiscard]] float Pz() const { return (*Source->Pz)[Entry]; }
    [[nodiscard]] float Energy() const { return (*Source->Energy)[Entry]; }

    [[nodiscard]] int Neg_Entry() const { return (*Source->Neg.McEntry)[Entry]; }
    [[nodiscard]] int Neg_PdgCode() const { return (*Source->Neg.PdgCode)[Entry]; }
    [[nodiscard]] int Neg_ReactionID() const { return (*Source->Neg.ReactionID)[Entry]; }
    [[nodiscard]] bool Neg_IsTrue() const { return static_cast<bool>((*Source->Neg.IsTrue)[Entry]); }
    [[nodiscard]] bool Neg_IsSignal() const { return static_cast<bool>((*Source->Neg.IsSignal)[Entry]); }
    [[nodiscard]] bool Neg_IsSecondary() const { return static_cast<bool>((*Source->Neg.IsSecondary)[Entry]); }

    [[nodiscard]] int Pos_Entry() const { return (*Source->Pos.McEntry)[Entry]; }
    [[nodiscard]] int Pos_PdgCode() const { return (*Source->Pos.PdgCode)[Entry]; }
    [[nodiscard]] int Pos_ReactionID() const { return (*Source->Pos.ReactionID)[Entry]; }
    [[nodiscard]] bool Pos_IsTrue() const { return static_cast<bool>((*Source->Pos.IsTrue)[Entry]); }
    [[nodiscard]] bool Pos_IsSignal() const { return static_cast<bool>((*Source->Pos.IsSignal)[Entry]); }
    [[nodiscard]] bool Pos_IsSecondary() const { return static_cast<bool>((*Source->Pos.IsSecondary)[Entry]); }

    [[nodiscard]] int Mother_Entry() const { return (*Source->Mother.McEntry)[Entry]; }
    [[nodiscard]] int Mother_PdgCode() const { return (*Source->Mother.PdgCode)[Entry]; }

    [[nodiscard]] float Neg_Px() const { return (*Source->Neg_Momentum.Px)[Entry]; }
    [[nodiscard]] float Neg_Py() const { return (*Source->Neg_Momentum.Py)[Entry]; }
    [[nodiscard]] float Neg_Pz() const { return (*Source->Neg_Momentum.Pz)[Entry]; }

    [[nodiscard]] float Pos_Px() const { return (*Source->Pos_Momentum.Px)[Entry]; }
    [[nodiscard]] float Pos_Py() const { return (*Source->Pos_Momentum.Py)[Entry]; }
    [[nodiscard]] float Pos_Pz() const { return (*Source->Pos_Momentum.Pz)[Entry]; }

    [[nodiscard]] float Decay_X() const { return (*Source->AtDecay.X)[Entry]; }
    [[nodiscard]] float Decay_Y() const { return (*Source->AtDecay.Y)[Entry]; }
    [[nodiscard]] float Decay_Z() const { return (*Source->AtDecay.Z)[Entry]; }

    [[nodiscard]] bool IsHybrid() const { return static_cast<bool>((*Source->IsHybrid)[Entry]); }
};

}  // namespace Tree2Secondaries::View::MC
