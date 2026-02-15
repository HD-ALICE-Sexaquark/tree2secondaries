#pragma once

#include <cstdlib>

#include "Storage/Vector/VectorTracks.hxx"
#include "View/BaseView.hxx"

namespace Tree2Secondaries::View::MC {

// Slight transformations to `View::MC::Particle`.
struct PackedTrack : View::Base<Storage::Vector::MC_Tracks, size_t> {

    PackedTrack() = delete;
    PackedTrack(const Storage::Vector::MC_Tracks* df, size_t entry)  //
        : View::Base<Storage::Vector::MC_Tracks, size_t>{df, entry} {}

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

    [[nodiscard]] int Mother_McEntry() const { return (*Source->Mother.McEntry)[Entry]; }
    [[nodiscard]] int Mother_PdgCode() const { return (*Source->Mother.PdgCode)[Entry]; }

    [[nodiscard]] int GrandMother_McEntry() const { return (*Source->GrandMother.McEntry)[Entry]; }
    [[nodiscard]] int GrandMother_PdgCode() const { return (*Source->GrandMother.PdgCode)[Entry]; }
};

}  // namespace Tree2Secondaries::View::MC
