#pragma once

#include <cstdlib>

#include "Storage/Vector/VectorTracks.hxx"
#include "View/BaseView.hxx"

namespace Tree2Secondaries::View::Rec {

struct Track : View::Base<Storage::Vector::Tracks, size_t> {

    Track() = delete;
    Track(const Storage::Vector::Tracks* df, size_t entry)  //
        : View::Base<Storage::Vector::Tracks, size_t>{df, entry} {}

    [[nodiscard]] float X() const { return (*Source->X)[Entry]; }
    [[nodiscard]] float Y() const { return (*Source->Y)[Entry]; }
    [[nodiscard]] float Z() const { return (*Source->Z)[Entry]; }
    [[nodiscard]] float Px() const { return (*Source->Px)[Entry]; }
    [[nodiscard]] float Py() const { return (*Source->Py)[Entry]; }
    [[nodiscard]] float Pz() const { return (*Source->Pz)[Entry]; }

    [[nodiscard]] float SigmaX2() const { return (*Source->SigmaX2)[Entry]; }
    [[nodiscard]] float SigmaXY() const { return (*Source->SigmaXY)[Entry]; }
    [[nodiscard]] float SigmaY2() const { return (*Source->SigmaY2)[Entry]; }
    [[nodiscard]] float SigmaXZ() const { return (*Source->SigmaXZ)[Entry]; }
    [[nodiscard]] float SigmaYZ() const { return (*Source->SigmaYZ)[Entry]; }
    [[nodiscard]] float SigmaZ2() const { return (*Source->SigmaZ2)[Entry]; }
    [[nodiscard]] float SigmaXPx() const { return (*Source->SigmaXPx)[Entry]; }
    [[nodiscard]] float SigmaYPx() const { return (*Source->SigmaYPx)[Entry]; }
    [[nodiscard]] float SigmaZPx() const { return (*Source->SigmaZPx)[Entry]; }
    [[nodiscard]] float SigmaPx2() const { return (*Source->SigmaPx2)[Entry]; }
    [[nodiscard]] float SigmaXPy() const { return (*Source->SigmaXPy)[Entry]; }
    [[nodiscard]] float SigmaYPy() const { return (*Source->SigmaYPy)[Entry]; }
    [[nodiscard]] float SigmaZPy() const { return (*Source->SigmaZPy)[Entry]; }
    [[nodiscard]] float SigmaPxPy() const { return (*Source->SigmaPxPy)[Entry]; }
    [[nodiscard]] float SigmaPy2() const { return (*Source->SigmaPy2)[Entry]; }
    [[nodiscard]] float SigmaXPz() const { return (*Source->SigmaXPz)[Entry]; }
    [[nodiscard]] float SigmaYPz() const { return (*Source->SigmaYPz)[Entry]; }
    [[nodiscard]] float SigmaZPz() const { return (*Source->SigmaZPz)[Entry]; }
    [[nodiscard]] float SigmaPxPz() const { return (*Source->SigmaPxPz)[Entry]; }
    [[nodiscard]] float SigmaPyPz() const { return (*Source->SigmaPyPz)[Entry]; }
    [[nodiscard]] float SigmaPz2() const { return (*Source->SigmaPz2)[Entry]; }

    [[nodiscard]] int Charge() const { return (*Source->Charge)[Entry]; }
    [[nodiscard]] float DCAxy() const { return (*Source->DCAxy)[Entry]; }
    [[nodiscard]] float DCAz() const { return (*Source->DCAz)[Entry]; }
    [[nodiscard]] float TPCSignal() const { return (*Source->TPCSignal)[Entry]; }
    [[nodiscard]] float NSigmaPion() const { return (*Source->NSigmaPion)[Entry]; }
    [[nodiscard]] float NSigmaKaon() const { return (*Source->NSigmaKaon)[Entry]; }
    [[nodiscard]] float NSigmaProton() const { return (*Source->NSigmaProton)[Entry]; }

    [[nodiscard]] int McEntry() const { return (*Source->McEntry)[Entry]; }  // NOTE: only valid when analyzing MC
    [[nodiscard]] size_t Index() const { return (*Source->Index)[Entry]; }   // NOTE: only valid when reading Packed Tracks
};

}  // namespace Tree2Secondaries::View::Rec
