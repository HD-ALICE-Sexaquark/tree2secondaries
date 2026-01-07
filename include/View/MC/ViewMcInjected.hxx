#pragma once

#include <cmath>

#include "Math/Constants.hxx"
#include "Storage/Vector/VectorInjected.hxx"
#include "Storage/Vector/VectorTracks.hxx"
#include "Storage/Vector/VectorV0s.hxx"
#include "View/BaseView.hxx"
#include "View/MC/ViewMcPackedTrack.hxx"
#include "View/MC/ViewMcPackedV0.hxx"

namespace Tree2Secondaries::View::MC {

struct Injected : View::Base<Storage::Vector::Injected> {

    Injected() = delete;
    Injected(const Storage::Vector::Injected* df, int entry)  //
        : View::Base<Storage::Vector::Injected>{.Source = df, .Entry = entry} {}

    [[nodiscard]] float Px() const { return IsValid() ? Source->Px->at(Entry) : Const::DummyFloat; };
    [[nodiscard]] float Py() const { return IsValid() ? Source->Py->at(Entry) : Const::DummyFloat; };
    [[nodiscard]] float Pz() const { return IsValid() ? Source->Pz->at(Entry) : Const::DummyFloat; };

    [[nodiscard]] float Pt2() const {  //
        return IsValid() ? Px() * Px() + Py() * Py() : Const::DummyFloat;
    };
    [[nodiscard]] float P2() const {  //
        return IsValid() ? Pt2() + Pz() * Pz() : Const::DummyFloat;
    };

    [[nodiscard]] float SV_X() const { return IsValid() ? Source->SV.X->at(Entry) : Const::DummyFloat; };
    [[nodiscard]] float SV_Y() const { return IsValid() ? Source->SV.Y->at(Entry) : Const::DummyFloat; };
    [[nodiscard]] float SV_Z() const { return IsValid() ? Source->SV.Z->at(Entry) : Const::DummyFloat; };

    [[nodiscard]] float Nucleon_Px() const { return IsValid() ? Source->Nucleon.Px->at(Entry) : Const::DummyFloat; };
    [[nodiscard]] float Nucleon_Py() const { return IsValid() ? Source->Nucleon.Py->at(Entry) : Const::DummyFloat; };
    [[nodiscard]] float Nucleon_Pz() const { return IsValid() ? Source->Nucleon.Pz->at(Entry) : Const::DummyFloat; };

    [[nodiscard]] float Nucleon_Pt2() const {  //
        return IsValid() ? Nucleon_Px() * Nucleon_Px() + Nucleon_Py() * Nucleon_Py() : Const::DummyFloat;
    };
    [[nodiscard]] float Nucleon_P2() const {  //
        return IsValid() ? Nucleon_Pt2() + Nucleon_Pz() * Nucleon_Pz() : Const::DummyFloat;
    };

    [[nodiscard]] int ReactionID() const { return IsValid() ? Entry + Const::ReactionID_Offset : Const::DummyInt; };
};

struct ChannelA : View::MC::Injected {

    View::MC::PackedV0 V0A;
    View::MC::PackedV0 V0B;

    ChannelA() = delete;
    ChannelA(const Storage::Vector::Injected* df, const Storage::Vector::MC_V0s* df_v0a, const Storage::Vector::MC_V0s* df_v0b, int v0a_entry,
             int v0b_entry)
        : View::MC::Injected{df, Const::DummyInt},  // overridden in definition
          V0A{df_v0a, v0a_entry},
          V0B{df_v0b, v0b_entry} {
        if (V0A.ReactionID() == V0B.ReactionID()) Entry = V0A.ReactionID() - Const::ReactionID_Offset;
    }
};

struct ChannelD : View::MC::Injected {

    View::MC::PackedV0 V0;
    View::MC::PackedTrack Kaon;

    ChannelD() = delete;
    ChannelD(const Storage::Vector::Injected* df, const Storage::Vector::MC_V0s* df_v0, const Storage::Vector::MC_Tracks* df_track, int v0_entry,
             int kaon_entry)
        : View::MC::Injected{df, Const::DummyInt},  // overridden in definition
          V0{df_v0, v0_entry},
          Kaon{df_track, kaon_entry} {
        if (V0.ReactionID() == Kaon.ReactionID()) Entry = V0.ReactionID() - Const::ReactionID_Offset;
    }
};

}  // namespace Tree2Secondaries::View::MC
