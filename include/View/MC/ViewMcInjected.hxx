#pragma once

#include "Math/Constants.hxx"
#include "Storage/Schema/SchemaVector.hxx"
#include "View/BaseView.hxx"
#include "View/MC/ViewMcPackedTrack.hxx"
#include "View/MC/ViewMcPackedV0.hxx"

namespace Tree2Secondaries::View::MC {

// NOTE: need to be guarded with `View::IsValid()` after construction.
struct Injected : View::Base<Schema::Vector::Injected, int> {

    Injected() = delete;
    Injected(const Schema::Vector::Injected* df, int entry)  //
        : View::Base<Schema::Vector::Injected, int>{df, entry} {}

    [[nodiscard]] float Px() const { return (*Source->mom.px)[EntryAsSize()]; }
    [[nodiscard]] float Py() const { return (*Source->mom.py)[EntryAsSize()]; }
    [[nodiscard]] float Pz() const { return (*Source->mom.pz)[EntryAsSize()]; }

    [[nodiscard]] float SV_X() const { return (*Source->sv.x)[EntryAsSize()]; }
    [[nodiscard]] float SV_Y() const { return (*Source->sv.y)[EntryAsSize()]; }
    [[nodiscard]] float SV_Z() const { return (*Source->sv.z)[EntryAsSize()]; }

    [[nodiscard]] float Nucleon_Px() const { return (*Source->mom_nucleon.px)[EntryAsSize()]; }
    [[nodiscard]] float Nucleon_Py() const { return (*Source->mom_nucleon.py)[EntryAsSize()]; }
    [[nodiscard]] float Nucleon_Pz() const { return (*Source->mom_nucleon.pz)[EntryAsSize()]; }

    [[nodiscard]] int ReactionID() const { return (*Source->reaction_id)[EntryAsSize()]; }
};

// NOTE: need to be guarded with `View::IsValid()` after construction.
struct ChannelA : View::MC::Injected {

    ChannelA() = delete;
    ChannelA(const Schema::Vector::Injected* df, const View::MC::PackedV0& v0a, const View::MC::PackedV0& v0b)
        : View::MC::Injected{df, Const::DummyInt},  // overridden in definition
          V0A{v0a},
          V0B{v0b} {
        if (V0A.ReactionID() != Const::DummyInt && V0A.ReactionID() == V0B.ReactionID()) {
            Entry = V0A.ReactionID() - Const::ReactionID_Offset;
        }
    }

    View::MC::PackedV0 V0A;
    View::MC::PackedV0 V0B;
};

// NOTE: need to be guarded with `View::IsValid()` after construction.
struct ChannelD : View::MC::Injected {

    ChannelD() = delete;
    ChannelD(const Schema::Vector::Injected* df, const View::MC::PackedTrack& ka, const View::MC::PackedV0& v0)
        : View::MC::Injected{df, Const::DummyInt},  // overridden in definition
          Kaon{ka},
          V0{v0} {
        if (V0.ReactionID() != Const::DummyInt && V0.ReactionID() == Kaon.ReactionID()) {
            Entry = V0.ReactionID() - Const::ReactionID_Offset;
        }
    }

    View::MC::PackedTrack Kaon;
    View::MC::PackedV0 V0;
};

}  // namespace Tree2Secondaries::View::MC
