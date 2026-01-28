#pragma once

#include <cmath>

#include "Storage/Vector/VectorTracks.hxx"
#include "View/BaseView.hxx"

namespace Tree2Secondaries::View::MC {

// Slight transformations to `View::MC::Particle`.
struct PackedTrack : View::Base<Storage::Vector::MC_Tracks> {

    PackedTrack() = delete;
    PackedTrack(const Storage::Vector::MC_Tracks* df, int entry)  //
        : View::Base<Storage::Vector::MC_Tracks>{.Source = df, .Entry = entry} {}

    [[nodiscard]] int McEntry() const { return Source->McEntry->at(Entry); }
    [[nodiscard]] int PdgCode() const { return Source->PdgCode->at(Entry); }
    [[nodiscard]] int ReactionID() const { return Source->ReactionID->at(Entry); }
    [[nodiscard]] bool IsTrue() const { return Source->IsTrue->at(Entry); }
    [[nodiscard]] bool IsSignal() const { return Source->IsSignal->at(Entry); }
    [[nodiscard]] bool IsSecondary() const { return Source->IsSecondary->at(Entry); }

    [[nodiscard]] float X() const { return Source->X->at(Entry); }
    [[nodiscard]] float Y() const { return Source->Y->at(Entry); }
    [[nodiscard]] float Z() const { return Source->Z->at(Entry); }
    [[nodiscard]] float Px() const { return Source->Px->at(Entry); }
    [[nodiscard]] float Py() const { return Source->Py->at(Entry); }
    [[nodiscard]] float Pz() const { return Source->Pz->at(Entry); }
    [[nodiscard]] float Energy() const { return Source->Energy->at(Entry); }

    [[nodiscard]] int Mother_Entry() const { return Source->Mother.McEntry->at(Entry); }
    [[nodiscard]] int Mother_PdgCode() const { return Source->Mother.PdgCode->at(Entry); }

    [[nodiscard]] int GrandMother_Entry() const { return Source->GrandMother.McEntry->at(Entry); }
    [[nodiscard]] int GrandMother_PdgCode() const { return Source->GrandMother.PdgCode->at(Entry); }
};

}  // namespace Tree2Secondaries::View::MC
