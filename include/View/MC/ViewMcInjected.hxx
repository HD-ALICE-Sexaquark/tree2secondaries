#pragma once

#include <cstdlib>

#include "Math/Constants.hxx"
#include "Storage/Vector/VectorInjected.hxx"
#include "Storage/Vector/VectorTracks.hxx"
#include "Storage/Vector/VectorV0s.hxx"
#include "View/BaseView.hxx"
#include "View/MC/ViewMcPackedTrack.hxx"
#include "View/MC/ViewMcPackedV0.hxx"

namespace Tree2Secondaries::View::MC {

// NOTE: need to be guarded with `View::IsValid()` after construction.
struct Injected : View::Base<Storage::Vector::Injected, int> {

    Injected() = delete;
    Injected(const Storage::Vector::Injected* df, int entry)  //
        : View::Base<Storage::Vector::Injected, int>{df, entry} {}

    [[nodiscard]] float Px() const { return (*Source->Px)[EntryAsSize()]; }
    [[nodiscard]] float Py() const { return (*Source->Py)[EntryAsSize()]; }
    [[nodiscard]] float Pz() const { return (*Source->Pz)[EntryAsSize()]; }

    [[nodiscard]] float SV_X() const { return (*Source->SV.X)[EntryAsSize()]; }
    [[nodiscard]] float SV_Y() const { return (*Source->SV.Y)[EntryAsSize()]; }
    [[nodiscard]] float SV_Z() const { return (*Source->SV.Z)[EntryAsSize()]; }

    [[nodiscard]] float Nucleon_Px() const { return (*Source->Nucleon.Px)[EntryAsSize()]; }
    [[nodiscard]] float Nucleon_Py() const { return (*Source->Nucleon.Py)[EntryAsSize()]; }
    [[nodiscard]] float Nucleon_Pz() const { return (*Source->Nucleon.Pz)[EntryAsSize()]; }
};

// NOTE: need to be guarded with `View::IsValid()` after construction.
struct ChannelA : View::MC::Injected {

    ChannelA() = delete;
    ChannelA(const Storage::Vector::Injected* df, const View::MC::PackedV0& v0a, const View::MC::PackedV0& v0b)
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
    ChannelD(const Storage::Vector::Injected* df, const View::MC::PackedTrack& ka, const View::MC::PackedV0& v0)
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
