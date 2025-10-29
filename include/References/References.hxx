#pragma once

#include <cmath>

#include "DataFormats/Events.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::Ref {

struct MC_Particle {

    [[nodiscard]] int Entry() const { return entry; };
    [[nodiscard]] bool FoundEntry() const { return Entry() > Const::DummyInt; }

    [[nodiscard]] float X() const { return FoundEntry() ? source->X->at(entry) : Const::DummyFloat; };
    [[nodiscard]] float Y() const { return FoundEntry() ? source->Y->at(entry) : Const::DummyFloat; };
    [[nodiscard]] float Z() const { return FoundEntry() ? source->Z->at(entry) : Const::DummyFloat; };
    [[nodiscard]] float Px() const { return FoundEntry() ? source->Px->at(entry) : Const::DummyFloat; };
    [[nodiscard]] float Py() const { return FoundEntry() ? source->Py->at(entry) : Const::DummyFloat; };
    [[nodiscard]] float Pz() const { return FoundEntry() ? source->Pz->at(entry) : Const::DummyFloat; };
    [[nodiscard]] float Energy() const { return FoundEntry() ? source->Energy->at(entry) : Const::DummyFloat; };

    [[nodiscard]] int PdgCode() const { return FoundEntry() ? source->PdgCode->at(entry) : Const::DummyInt; };
    [[nodiscard]] int MotherEntry() const { return FoundEntry() ? source->MotherEntry->at(entry) : Const::DummyInt; };
    [[nodiscard]] int Generator() const { return FoundEntry() ? source->Generator->at(entry) : Const::DummyInt; };
    [[nodiscard]] int Status() const { return FoundEntry() ? source->Status->at(entry) : Const::DummyInt; };
    [[nodiscard]] int IsSecFromMat() const { return FoundEntry() ? source->IsSecFromMat->at(entry) : Const::DummyInt; };
    [[nodiscard]] int IsSecFromWeak() const { return FoundEntry() ? source->IsSecFromWeak->at(entry) : Const::DummyInt; };

    const DF::Events::MCParticles* source{};
    int entry{};
    EParticle hyp{};  // hypothesis
};

struct MC_Track : Ref::MC_Particle {
    [[nodiscard]] bool AsTrack_IsTrue() const {  //
        return Const::Particle_PdgCode[hyp] == PdgCode();
    };
    [[nodiscard]] bool AsTrack_IsSignal() const {  //
        return AsTrack_IsTrue() && Generator() == 2;
    };
    [[nodiscard]] bool AsTrack_IsSecondary() const {  //
        return IsSecFromMat() || IsSecFromWeak() || AsTrack_IsSignal();
    };
    [[nodiscard]] int AsTrack_ReactionID(const Ref::MC_Particle& mother) const {  //
        return AsTrack_IsSignal() ? mother.Status() : Const::DummyInt;
    };
};

struct MC_V0 : Ref::MC_Particle {

    MC_V0() = delete;
    MC_V0(const Ref::MC_Track& neg_ref, const Ref::MC_Track& pos_ref, EParticle pid)
        : Ref::MC_Particle{.source = neg_ref.source,  // NOTE: same `source`
                           .entry = Const::DummyInt,  // NOTE: will be overriden
                           .hyp = pid},
          neg{neg_ref},
          pos{pos_ref} {
        if (DaughtersFromSameMother()) entry = neg.MotherEntry();  // override `entry`
    }

    [[nodiscard]] bool AsV0_IsTrue() const {  //
        return Const::Particle_PdgCode[hyp] == PdgCode() && neg.AsTrack_IsTrue() && pos.AsTrack_IsTrue();
    };
    [[nodiscard]] bool AsV0_IsSignal() const {  //
        return AsV0_IsTrue() && Generator() == 2 && neg.AsTrack_ReactionID(*this) == pos.AsTrack_ReactionID(*this);
    };
    [[nodiscard]] bool AsV0_IsSecondary() const {  //
        return IsSecFromMat() || IsSecFromWeak() || AsV0_IsSignal();
    };
    [[nodiscard]] int AsV0_ReactionID() const {  //
        return AsV0_IsSignal() ? Status() : Const::DummyInt;
    };
    [[nodiscard]] bool AsV0_IsHybrid() const {
        return !AsV0_IsSignal() &&                                      //
               ((neg.AsTrack_IsSignal() && !pos.AsTrack_IsSignal()) ||  //
                (!neg.AsTrack_IsSignal() && pos.AsTrack_IsSignal()));
    }

    [[nodiscard]] bool DaughtersFromSameMother() const { return neg.MotherEntry() == pos.MotherEntry() && neg.MotherEntry() > Const::DummyInt; };

    [[nodiscard]] float DecayX() const { return DaughtersFromSameMother() ? neg.X() : Const::DummyFloat; }
    [[nodiscard]] float DecayY() const { return DaughtersFromSameMother() ? neg.Y() : Const::DummyFloat; }
    [[nodiscard]] float DecayZ() const { return DaughtersFromSameMother() ? neg.Z() : Const::DummyFloat; }

    Ref::MC_Track neg{};
    Ref::MC_Track pos{};
};

}  // namespace Tree2Secondaries::Ref
