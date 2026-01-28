#pragma once

#include <array>
#include <cmath>

#include "Storage/Vector/VectorTracks.hxx"
#include "View/BaseView.hxx"

namespace Tree2Secondaries::View::Rec {

struct Track : View::Base<Storage::Vector::Tracks> {

    Track() = delete;
    Track(const Storage::Vector::Tracks* df, int entry)  //
        : View::Base<Storage::Vector::Tracks>{.Source = df, .Entry = entry} {}

    [[nodiscard]] float X() const { return Source->X->at(Entry); }
    [[nodiscard]] float Y() const { return Source->Y->at(Entry); }
    [[nodiscard]] float Z() const { return Source->Z->at(Entry); }
    [[nodiscard]] float Px() const { return Source->Px->at(Entry); }
    [[nodiscard]] float Py() const { return Source->Py->at(Entry); }
    [[nodiscard]] float Pz() const { return Source->Pz->at(Entry); }

    [[nodiscard]] std::array<float, 6> State_NoE() const { return {X(), Y(), Z(), Px(), Py(), Pz()}; }
    [[nodiscard]] float Pt2() const { return Px() * Px() + Py() * Py(); }
    [[nodiscard]] float Pt() const { return std::sqrt(Pt2()); }
    [[nodiscard]] float P2() const { return Pt2() + Pz() * Pz(); }
    [[nodiscard]] float P() const { return std::sqrt(P2()); }

    [[nodiscard]] float SigmaX2() const { return Source->SigmaX2->at(Entry); }
    [[nodiscard]] float SigmaXY() const { return Source->SigmaXY->at(Entry); }
    [[nodiscard]] float SigmaY2() const { return Source->SigmaY2->at(Entry); }
    [[nodiscard]] float SigmaXZ() const { return Source->SigmaXZ->at(Entry); }
    [[nodiscard]] float SigmaYZ() const { return Source->SigmaYZ->at(Entry); }
    [[nodiscard]] float SigmaZ2() const { return Source->SigmaZ2->at(Entry); }
    [[nodiscard]] float SigmaXPx() const { return Source->SigmaXPx->at(Entry); }
    [[nodiscard]] float SigmaYPx() const { return Source->SigmaYPx->at(Entry); }
    [[nodiscard]] float SigmaZPx() const { return Source->SigmaZPx->at(Entry); }
    [[nodiscard]] float SigmaPx2() const { return Source->SigmaPx2->at(Entry); }
    [[nodiscard]] float SigmaXPy() const { return Source->SigmaXPy->at(Entry); }
    [[nodiscard]] float SigmaYPy() const { return Source->SigmaYPy->at(Entry); }
    [[nodiscard]] float SigmaZPy() const { return Source->SigmaZPy->at(Entry); }
    [[nodiscard]] float SigmaPxPy() const { return Source->SigmaPxPy->at(Entry); }
    [[nodiscard]] float SigmaPy2() const { return Source->SigmaPy2->at(Entry); }
    [[nodiscard]] float SigmaXPz() const { return Source->SigmaXPz->at(Entry); }
    [[nodiscard]] float SigmaYPz() const { return Source->SigmaYPz->at(Entry); }
    [[nodiscard]] float SigmaZPz() const { return Source->SigmaZPz->at(Entry); }
    [[nodiscard]] float SigmaPxPz() const { return Source->SigmaPxPz->at(Entry); }
    [[nodiscard]] float SigmaPyPz() const { return Source->SigmaPyPz->at(Entry); }
    [[nodiscard]] float SigmaPz2() const { return Source->SigmaPz2->at(Entry); }

    [[nodiscard]] std::array<float, 21> CovMatrix_NoE() const {
        return {SigmaX2(),                                                    //
                SigmaXY(),  SigmaY2(),                                        //
                SigmaXZ(),  SigmaYZ(),  SigmaZ2(),                            //
                SigmaXPx(), SigmaYPx(), SigmaZPx(), SigmaPx2(),               //
                SigmaXPy(), SigmaYPy(), SigmaZPy(), SigmaPxPy(), SigmaPy2(),  //
                SigmaXPz(), SigmaYPz(), SigmaZPz(), SigmaPxPz(), SigmaPyPz(), SigmaPz2()};
    }

    [[nodiscard]] int Charge() const { return Source->Charge->at(Entry); }
    [[nodiscard]] float DCAxy() const { return Source->DCAxy->at(Entry); }
    [[nodiscard]] float DCAz() const { return Source->DCAz->at(Entry); }
    [[nodiscard]] float TPCSignal() const { return Source->TPCSignal->at(Entry); }
    [[nodiscard]] float NSigmaPion() const { return Source->NSigmaPion->at(Entry); }
    [[nodiscard]] float NSigmaKaon() const { return Source->NSigmaKaon->at(Entry); }
    [[nodiscard]] float NSigmaProton() const { return Source->NSigmaProton->at(Entry); }

    [[nodiscard]] int McEntry() const { return Source->McEntry->at(Entry); }  // NOTE: only valid when analyzing MC
    [[nodiscard]] int Index() const { return Source->Index->at(Entry); }      // NOTE: only valid when reading Packed Tracks
};

}  // namespace Tree2Secondaries::View::Rec
