#pragma once

#include <cmath>

#include "DataFormats/Injected.hxx"
#include "DataFormats/PackedEvents.hxx"
#include "DataFormats/StructsOfVectors.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::Ref {

// Mirror of SOV //

struct PackedBachelor {

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

struct PackedV0 {

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

struct Injected {

    [[nodiscard]] int Entry() const { return entry; };
    [[nodiscard]] bool FoundEntry() const { return Entry() > Const::DummyInt; }

    [[nodiscard]] int ReactionID() const { return FoundEntry() ? entry + 600 : Const::DummyInt; };
    [[nodiscard]] double Mass() const { return FoundEntry() ? mass : Const::DummyDouble; };

    [[nodiscard]] float BeforePx() const { return FoundEntry() ? source->Px->at(entry) : Const::DummyFloat; };
    [[nodiscard]] float BeforePy() const { return FoundEntry() ? source->Py->at(entry) : Const::DummyFloat; };
    [[nodiscard]] float BeforePz() const { return FoundEntry() ? source->Pz->at(entry) : Const::DummyFloat; };
    [[nodiscard]] double BeforeP2() const {
        return FoundEntry() ? BeforePx() * BeforePx() + BeforePy() * BeforePy() + BeforePz() * BeforePz() : Const::DummyDouble;
    };
    [[nodiscard]] double BeforeE() const { return FoundEntry() ? std::sqrt(BeforeP2() + Mass() * Mass()) : Const::DummyDouble; };

    [[nodiscard]] float NucleonPx() const { return FoundEntry() ? source->Nucleon.Px->at(entry) : Const::DummyFloat; };
    [[nodiscard]] float NucleonPy() const { return FoundEntry() ? source->Nucleon.Py->at(entry) : Const::DummyFloat; };
    [[nodiscard]] float NucleonPz() const { return FoundEntry() ? source->Nucleon.Pz->at(entry) : Const::DummyFloat; };
    [[nodiscard]] double NucleonP2() const {
        return FoundEntry() ? NucleonPx() * NucleonPx() + NucleonPy() * NucleonPy() + NucleonPz() * NucleonPz() : Const::DummyDouble;
    };
    [[nodiscard]] double NucleonMass() const { return FoundEntry() ? Const::Particle_Mass[nucleon_pid] : Const::DummyDouble; };
    [[nodiscard]] double NucleonE() const { return FoundEntry() ? std::sqrt(NucleonP2() + NucleonMass() * NucleonMass()) : Const::DummyDouble; };

    [[nodiscard]] float SV_X() const { return FoundEntry() ? source->X->at(entry) : Const::DummyFloat; };
    [[nodiscard]] float SV_Y() const { return FoundEntry() ? source->Y->at(entry) : Const::DummyFloat; };
    [[nodiscard]] float SV_Z() const { return FoundEntry() ? source->Z->at(entry) : Const::DummyFloat; };

    const DF::SOV::Injected* source{};
    double mass{};
    EParticle nucleon_pid{};
    int entry{Const::DummyInt};
};

struct ChannelA : Ref::Injected {

    ChannelA() = delete;
    ChannelA(const Ref::Injected& inj, const Ref::PackedV0& v0a_ref, const Ref::PackedV0& v0b_ref)  //
        : Ref::Injected{inj},                                                                       //
          v0a{v0a_ref},
          v0b{v0b_ref} {
        if (DaughtersFromSameReaction()) entry = v0a.ReactionID() - 600;  // override `entry`
    }

    [[nodiscard]] bool DaughtersFromSameReaction() const {  //
        return v0a.ReactionID() == v0b.ReactionID() && v0a.ReactionID() > Const::DummyInt;
    };

    [[nodiscard]] double AsChannelA_AfterPx() const { return v0a.Px() + v0b.Px(); };
    [[nodiscard]] double AsChannelA_AfterPy() const { return v0a.Py() + v0b.Py(); };
    [[nodiscard]] double AsChannelA_AfterPz() const { return v0a.Pz() + v0b.Pz(); };
    [[nodiscard]] double AsChannelA_AfterE() const { return v0a.Energy() + v0b.Energy(); };

    [[nodiscard]] int AsChannelA_ReactionID() const { return DaughtersFromSameReaction() ? entry + 600 : Const::DummyInt; }
    [[nodiscard]] bool AsChannelA_IsSignal() const { return v0a.IsSignal() && v0b.IsSignal() && DaughtersFromSameReaction(); };
    [[nodiscard]] bool AsChannelA_IsHybrid() const {
        return !AsChannelA_IsSignal() &&
               (v0a.IsHybrid() || v0b.IsHybrid() || (v0a.IsSignal() && !v0b.IsSignal()) || (!v0a.IsSignal() && v0b.IsSignal()) ||
                (v0a.IsSignal() && v0b.IsSignal() && v0a.ReactionID() != v0b.ReactionID()));
    };

    Ref::PackedV0 v0a;
    Ref::PackedV0 v0b;
};

struct ChannelD : Ref::Injected {

    ChannelD() = delete;
    ChannelD(const Ref::Injected& inj, const Ref::PackedV0& v0_ref, const Ref::PackedBachelor& kaon_ref)  //
        : Ref::Injected{inj},                                                                             //
          v0{v0_ref},
          kaon{kaon_ref} {
        if (DaughtersFromSameReaction()) entry = v0.ReactionID() - 600;  // override `entry`
    }

    [[nodiscard]] bool DaughtersFromSameReaction() const {  //
        return v0.ReactionID() == kaon.ReactionID() && v0.ReactionID() > Const::DummyInt;
    };

    [[nodiscard]] double AsChannelD_AfterPx() const { return v0.Px() + v0.Px(); };
    [[nodiscard]] double AsChannelD_AfterPy() const { return v0.Py() + v0.Py(); };
    [[nodiscard]] double AsChannelD_AfterPz() const { return v0.Pz() + v0.Pz(); };
    [[nodiscard]] double AsChannelD_AfterE() const { return v0.Energy() + v0.Energy(); };

    [[nodiscard]] int AsChannelD_ReactionID() const { return DaughtersFromSameReaction() ? entry + 600 : Const::DummyInt; }
    [[nodiscard]] bool AsChannelD_IsSignal() const { return v0.IsSignal() && kaon.IsSignal() && DaughtersFromSameReaction(); };
    [[nodiscard]] bool AsChannelD_IsHybrid() const {
        return !AsChannelD_IsSignal() && (v0.IsHybrid() || (v0.IsSignal() && !kaon.IsSignal()) || (!v0.IsSignal() && kaon.IsSignal()) ||
                                          (v0.IsSignal() && kaon.IsSignal() && v0.ReactionID() != kaon.ReactionID()));
    };

    Ref::PackedV0 v0;
    Ref::PackedBachelor kaon;
};

}  // namespace Tree2Secondaries::Ref
