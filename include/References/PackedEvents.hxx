#pragma once

#include <cmath>

#include "DataFormats/PackedEvents.hxx"
#include "DataFormats/StructsOfVectors.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::Ref {

// Mirror of SOV //

struct alignas(T2S_SIMD_ALIGN) PackedBachelor {

    [[nodiscard]] int Entry() const { return source->Entry->at(entry); };
    [[nodiscard]] int PdgCode() const { return source->PdgCode->at(entry); };
    [[nodiscard]] int ReactionID() const { return source->ReactionID->at(entry); };
    [[nodiscard]] bool IsTrue() const { return source->IsTrue->at(entry); };
    [[nodiscard]] bool IsSignal() const { return source->IsSignal->at(entry); };
    [[nodiscard]] bool IsSecondary() const { return source->IsSecondary->at(entry); };

    [[nodiscard]] float Px() const { return source->Px->at(entry); };
    [[nodiscard]] float Py() const { return source->Py->at(entry); };
    [[nodiscard]] float Pz() const { return source->Pz->at(entry); };
    [[nodiscard]] float Energy() const { return source->Energy->at(entry); };

    [[nodiscard]] int Mother_Entry() const { return source->Mother_Entry->at(entry); };
    [[nodiscard]] int Mother_PdgCode() const { return source->Mother_PdgCode->at(entry); };

    [[nodiscard]] int GrandMother_Entry() const { return source->GrandMother_Entry->at(entry); };
    [[nodiscard]] int GrandMother_PdgCode() const { return source->GrandMother_PdgCode->at(entry); };

    const DF::Packed::LinkedTracks* source{};
    int entry{};
};

struct alignas(T2S_SIMD_ALIGN) PackedV0 {

    [[nodiscard]] int Entry() const { return source->Entry->at(entry); };
    [[nodiscard]] int PdgCode() const { return source->PdgCode->at(entry); };
    [[nodiscard]] int ReactionID() const { return source->ReactionID->at(entry); };
    [[nodiscard]] bool IsTrue() const { return source->IsTrue->at(entry); };
    [[nodiscard]] bool IsSignal() const { return source->IsSignal->at(entry); };
    [[nodiscard]] bool IsSecondary() const { return source->IsSecondary->at(entry); };
    [[nodiscard]] float Px() const { return source->Px->at(entry); };
    [[nodiscard]] float Py() const { return source->Py->at(entry); };
    [[nodiscard]] float Pz() const { return source->Pz->at(entry); };
    [[nodiscard]] float Energy() const { return source->Energy->at(entry); };
    [[nodiscard]] int Mother_Entry() const { return source->Mother_Entry->at(entry); };
    [[nodiscard]] int Mother_PdgCode() const { return source->Mother_PdgCode->at(entry); };

    [[nodiscard]] int Neg_Entry() const { return source->Neg.Entry->at(entry); };
    [[nodiscard]] int Neg_PdgCode() const { return source->Neg.PdgCode->at(entry); };
    [[nodiscard]] int Neg_ReactionID() const { return source->Neg.ReactionID->at(entry); };
    [[nodiscard]] bool Neg_IsTrue() const { return source->Neg.IsTrue->at(entry); };
    [[nodiscard]] bool Neg_IsSignal() const { return source->Neg.IsSignal->at(entry); };
    [[nodiscard]] bool Neg_IsSecondary() const { return source->Neg.IsSecondary->at(entry); };
    [[nodiscard]] float Neg_Px() const { return source->Neg.Px->at(entry); };
    [[nodiscard]] float Neg_Py() const { return source->Neg.Py->at(entry); };
    [[nodiscard]] float Neg_Pz() const { return source->Neg.Pz->at(entry); };

    [[nodiscard]] int Pos_Entry() const { return source->Pos.Entry->at(entry); };
    [[nodiscard]] int Pos_PdgCode() const { return source->Pos.PdgCode->at(entry); };
    [[nodiscard]] int Pos_ReactionID() const { return source->Pos.ReactionID->at(entry); };
    [[nodiscard]] bool Pos_IsTrue() const { return source->Pos.IsTrue->at(entry); };
    [[nodiscard]] bool Pos_IsSignal() const { return source->Pos.IsSignal->at(entry); };
    [[nodiscard]] bool Pos_IsSecondary() const { return source->Pos.IsSecondary->at(entry); };
    [[nodiscard]] float Pos_Px() const { return source->Pos.Px->at(entry); };
    [[nodiscard]] float Pos_Py() const { return source->Pos.Py->at(entry); };
    [[nodiscard]] float Pos_Pz() const { return source->Pos.Pz->at(entry); };

    [[nodiscard]] float DecayX() const { return source->AtDecay.X->at(entry); };
    [[nodiscard]] float DecayY() const { return source->AtDecay.Y->at(entry); };
    [[nodiscard]] float DecayZ() const { return source->AtDecay.Z->at(entry); };

    [[nodiscard]] int IsHybrid() const { return source->IsHybrid->at(entry); };

    const DF::Packed::LinkedV0s* source{};
    int entry{};
};

}  // namespace Tree2Secondaries::Ref
