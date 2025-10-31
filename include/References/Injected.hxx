#pragma once

#include <cmath>

#include "DataFormats/Injected.hxx"
#include "DataFormats/StructsOfVectors.hxx"
#include "Math/Constants.hxx"
#include "References/PackedEvents.hxx"

namespace Tree2Secondaries::Ref {

struct alignas(T2S_SIMD_ALIGN) Injected {

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

struct alignas(T2S_SIMD_ALIGN) ChannelA : Ref::Injected {

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

struct alignas(T2S_SIMD_ALIGN) ChannelD : Ref::Injected {

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
