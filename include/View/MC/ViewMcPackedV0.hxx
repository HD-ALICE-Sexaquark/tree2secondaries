#pragma once

#include <cmath>

#include "Storage/Vector/VectorV0s.hxx"
#include "View/BaseView.hxx"

namespace Tree2Secondaries::View::MC {

struct PackedV0 : View::Base<Storage::Vector::MC_V0s> {

    PackedV0() = delete;
    PackedV0(const Storage::Vector::MC_V0s* df, int entry)  //
        : View::Base<Storage::Vector::MC_V0s>{.Source = df, .Entry = entry} {}

    [[nodiscard]] int McEntry() const { return Source->McEntry->at(Entry); };
    [[nodiscard]] int PdgCode() const { return Source->PdgCode->at(Entry); };
    [[nodiscard]] int ReactionID() const { return Source->ReactionID->at(Entry); };
    [[nodiscard]] bool IsTrue() const { return static_cast<bool>(Source->IsTrue->at(Entry)); };
    [[nodiscard]] bool IsSignal() const { return static_cast<bool>(Source->IsSignal->at(Entry)); };
    [[nodiscard]] bool IsSecondary() const { return static_cast<bool>(Source->IsSecondary->at(Entry)); };

    [[nodiscard]] float X() const { return Source->X->at(Entry); };
    [[nodiscard]] float Y() const { return Source->Y->at(Entry); };
    [[nodiscard]] float Z() const { return Source->Z->at(Entry); };
    [[nodiscard]] float Px() const { return Source->Px->at(Entry); };
    [[nodiscard]] float Py() const { return Source->Py->at(Entry); };
    [[nodiscard]] float Pz() const { return Source->Pz->at(Entry); };
    [[nodiscard]] float Energy() const { return Source->Energy->at(Entry); };

    [[nodiscard]] int Neg_Entry() const { return Source->Neg.McEntry->at(Entry); };
    [[nodiscard]] int Neg_PdgCode() const { return Source->Neg.PdgCode->at(Entry); };
    [[nodiscard]] int Neg_ReactionID() const { return Source->Neg.ReactionID->at(Entry); };
    [[nodiscard]] bool Neg_IsTrue() const { return static_cast<bool>(Source->Neg.IsTrue->at(Entry)); };
    [[nodiscard]] bool Neg_IsSignal() const { return static_cast<bool>(Source->Neg.IsSignal->at(Entry)); };
    [[nodiscard]] bool Neg_IsSecondary() const { return static_cast<bool>(Source->Neg.IsSecondary->at(Entry)); };

    [[nodiscard]] int Pos_Entry() const { return Source->Pos.McEntry->at(Entry); };
    [[nodiscard]] int Pos_PdgCode() const { return Source->Pos.PdgCode->at(Entry); };
    [[nodiscard]] int Pos_ReactionID() const { return Source->Pos.ReactionID->at(Entry); };
    [[nodiscard]] bool Pos_IsTrue() const { return static_cast<bool>(Source->Pos.IsTrue->at(Entry)); };
    [[nodiscard]] bool Pos_IsSignal() const { return static_cast<bool>(Source->Pos.IsSignal->at(Entry)); };
    [[nodiscard]] bool Pos_IsSecondary() const { return static_cast<bool>(Source->Pos.IsSecondary->at(Entry)); };

    [[nodiscard]] int Mother_Entry() const { return Source->Mother.McEntry->at(Entry); };
    [[nodiscard]] int Mother_PdgCode() const { return Source->Mother.PdgCode->at(Entry); };

    [[nodiscard]] float Neg_Px() const { return Source->Neg_Momentum.Px->at(Entry); };
    [[nodiscard]] float Neg_Py() const { return Source->Neg_Momentum.Py->at(Entry); };
    [[nodiscard]] float Neg_Pz() const { return Source->Neg_Momentum.Pz->at(Entry); };

    [[nodiscard]] float Pos_Px() const { return Source->Pos_Momentum.Px->at(Entry); };
    [[nodiscard]] float Pos_Py() const { return Source->Pos_Momentum.Py->at(Entry); };
    [[nodiscard]] float Pos_Pz() const { return Source->Pos_Momentum.Pz->at(Entry); };

    [[nodiscard]] float Decay_X() const { return Source->AtDecay.X->at(Entry); };
    [[nodiscard]] float Decay_Y() const { return Source->AtDecay.Y->at(Entry); };
    [[nodiscard]] float Decay_Z() const { return Source->AtDecay.Z->at(Entry); };

    [[nodiscard]] bool IsHybrid() const { return static_cast<bool>(Source->IsHybrid->at(Entry)); };
};

}  // namespace Tree2Secondaries::View::MC
