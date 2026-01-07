#pragma once

#include <cmath>

#include <KFParticle_Math.hxx>

#include "Storage/Vector/VectorV0s.hxx"
#include "View/BaseView.hxx"
#include "View/Reconstructed/ViewTrack.hxx"

namespace Tree2Secondaries::View::Rec {

struct V0 : View::Base<Storage::Vector::V0s> {

    View::Rec::Track Neg;
    View::Rec::Track Pos;

    V0() = delete;
    V0(const Storage::Vector::V0s* df, int entry)
        : View::Base<Storage::Vector::V0s>{.Source = df, .Entry = entry},
          Neg{&df->Neg, entry},  // NOTE: yes, same `entry`.
          Pos{&df->Pos, entry}   // NOTE: yes, same `entry`.
    {}

    [[nodiscard]] float X() const { return Source->X->at(Entry); };
    [[nodiscard]] float Y() const { return Source->Y->at(Entry); };
    [[nodiscard]] float Z() const { return Source->Z->at(Entry); };
    [[nodiscard]] float Px() const { return Source->Px->at(Entry); };
    [[nodiscard]] float Py() const { return Source->Py->at(Entry); };
    [[nodiscard]] float Pz() const { return Source->Pz->at(Entry); };
    [[nodiscard]] float Energy() const { return Source->Energy->at(Entry); };

    [[nodiscard]] KF::Vector<7> State() const { return {X(), Y(), Z(), Px(), Py(), Pz(), Energy()}; };

    [[nodiscard]] float SigmaX2() const { return Source->SigmaX2->at(Entry); };
    [[nodiscard]] float SigmaXY() const { return Source->SigmaXY->at(Entry); };
    [[nodiscard]] float SigmaY2() const { return Source->SigmaY2->at(Entry); };
    [[nodiscard]] float SigmaXZ() const { return Source->SigmaXZ->at(Entry); };
    [[nodiscard]] float SigmaYZ() const { return Source->SigmaYZ->at(Entry); };
    [[nodiscard]] float SigmaZ2() const { return Source->SigmaZ2->at(Entry); };
    [[nodiscard]] float SigmaXPx() const { return Source->SigmaXPx->at(Entry); };
    [[nodiscard]] float SigmaYPx() const { return Source->SigmaYPx->at(Entry); };
    [[nodiscard]] float SigmaZPx() const { return Source->SigmaZPx->at(Entry); };
    [[nodiscard]] float SigmaPx2() const { return Source->SigmaPx2->at(Entry); };
    [[nodiscard]] float SigmaXPy() const { return Source->SigmaXPy->at(Entry); };
    [[nodiscard]] float SigmaYPy() const { return Source->SigmaYPy->at(Entry); };
    [[nodiscard]] float SigmaZPy() const { return Source->SigmaZPy->at(Entry); };
    [[nodiscard]] float SigmaPxPy() const { return Source->SigmaPxPy->at(Entry); };
    [[nodiscard]] float SigmaPy2() const { return Source->SigmaPy2->at(Entry); };
    [[nodiscard]] float SigmaXPz() const { return Source->SigmaXPz->at(Entry); };
    [[nodiscard]] float SigmaYPz() const { return Source->SigmaYPz->at(Entry); };
    [[nodiscard]] float SigmaZPz() const { return Source->SigmaZPz->at(Entry); };
    [[nodiscard]] float SigmaPxPz() const { return Source->SigmaPxPz->at(Entry); };
    [[nodiscard]] float SigmaPyPz() const { return Source->SigmaPyPz->at(Entry); };
    [[nodiscard]] float SigmaPz2() const { return Source->SigmaPz2->at(Entry); };
    [[nodiscard]] float SigmaXE() const { return Source->SigmaXE->at(Entry); };
    [[nodiscard]] float SigmaYE() const { return Source->SigmaYE->at(Entry); };
    [[nodiscard]] float SigmaZE() const { return Source->SigmaZE->at(Entry); };
    [[nodiscard]] float SigmaPxE() const { return Source->SigmaPxE->at(Entry); };
    [[nodiscard]] float SigmaPyE() const { return Source->SigmaPyE->at(Entry); };
    [[nodiscard]] float SigmaPzE() const { return Source->SigmaPzE->at(Entry); };
    [[nodiscard]] float SigmaE2() const { return Source->SigmaE2->at(Entry); };

    [[nodiscard]] KF::SymMatrix<7> CovMatrix() const {
        return {SigmaX2(),                                                                 //
                SigmaXY(),  SigmaY2(),                                                     //
                SigmaXZ(),  SigmaYZ(),  SigmaZ2(),                                         //
                SigmaXPx(), SigmaYPx(), SigmaZPx(), SigmaPx2(),                            //
                SigmaXPy(), SigmaYPy(), SigmaZPy(), SigmaPxPy(), SigmaPy2(),               //
                SigmaXPz(), SigmaYPz(), SigmaZPz(), SigmaPxPz(), SigmaPyPz(), SigmaPz2(),  //
                SigmaXE(),  SigmaYE(),  SigmaZE(),  SigmaPxE(),  SigmaPyE(),  SigmaPzE(), SigmaE2()};
    };

    [[nodiscard]] float Neg_atV0_X() const { return Source->Neg_atPCA.X->at(Entry); };
    [[nodiscard]] float Neg_atV0_Y() const { return Source->Neg_atPCA.Y->at(Entry); };
    [[nodiscard]] float Neg_atV0_Z() const { return Source->Neg_atPCA.Z->at(Entry); };
    [[nodiscard]] float Neg_atV0_Px() const { return Source->Neg_atPCA.Px->at(Entry); };
    [[nodiscard]] float Neg_atV0_Py() const { return Source->Neg_atPCA.Py->at(Entry); };
    [[nodiscard]] float Neg_atV0_Pz() const { return Source->Neg_atPCA.Pz->at(Entry); };

    [[nodiscard]] float Pos_atV0_X() const { return Source->Pos_atPCA.X->at(Entry); };
    [[nodiscard]] float Pos_atV0_Y() const { return Source->Pos_atPCA.Y->at(Entry); };
    [[nodiscard]] float Pos_atV0_Z() const { return Source->Pos_atPCA.Z->at(Entry); };
    [[nodiscard]] float Pos_atV0_Px() const { return Source->Pos_atPCA.Px->at(Entry); };
    [[nodiscard]] float Pos_atV0_Py() const { return Source->Pos_atPCA.Py->at(Entry); };
    [[nodiscard]] float Pos_atV0_Pz() const { return Source->Pos_atPCA.Pz->at(Entry); };

    [[nodiscard]] float Chi2NDF() const { return Source->Chi2NDF->at(Entry); };
};

}  // namespace Tree2Secondaries::View::Rec
