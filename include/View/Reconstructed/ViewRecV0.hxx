#pragma once

#include <cstdlib>

#include "Storage/Vector/VectorV0s.hxx"
#include "View/BaseView.hxx"
#include "View/Reconstructed/ViewRecTrack.hxx"

namespace Tree2Secondaries::View::Rec {

struct V0 : View::Base<Storage::Vector::V0s, size_t> {

    V0() = delete;
    V0(const Storage::Vector::V0s* df, size_t entry)
        : View::Base<Storage::Vector::V0s, size_t>{df, entry},
          Neg{&df->Neg, entry},  // NOTE: yes, same `entry
          Pos{&df->Pos, entry}   // NOTE: yes, same `entry
    {}

    View::Rec::Track Neg;
    View::Rec::Track Pos;

    [[nodiscard]] float X() const { return (*Source->X)[Entry]; }
    [[nodiscard]] float Y() const { return (*Source->Y)[Entry]; }
    [[nodiscard]] float Z() const { return (*Source->Z)[Entry]; }
    [[nodiscard]] float Px() const { return (*Source->Px)[Entry]; }
    [[nodiscard]] float Py() const { return (*Source->Py)[Entry]; }
    [[nodiscard]] float Pz() const { return (*Source->Pz)[Entry]; }
    [[nodiscard]] float Energy() const { return (*Source->Energy)[Entry]; }

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

    [[nodiscard]] float SigmaXE() const { return (*Source->SigmaXE)[Entry]; }
    [[nodiscard]] float SigmaYE() const { return (*Source->SigmaYE)[Entry]; }
    [[nodiscard]] float SigmaZE() const { return (*Source->SigmaZE)[Entry]; }
    [[nodiscard]] float SigmaPxE() const { return (*Source->SigmaPxE)[Entry]; }
    [[nodiscard]] float SigmaPyE() const { return (*Source->SigmaPyE)[Entry]; }
    [[nodiscard]] float SigmaPzE() const { return (*Source->SigmaPzE)[Entry]; }
    [[nodiscard]] float SigmaE2() const { return (*Source->SigmaE2)[Entry]; }

    [[nodiscard]] float Neg_atPCA_X() const { return (*Source->Neg_atPCA.X)[Entry]; }
    [[nodiscard]] float Neg_atPCA_Y() const { return (*Source->Neg_atPCA.Y)[Entry]; }
    [[nodiscard]] float Neg_atPCA_Z() const { return (*Source->Neg_atPCA.Z)[Entry]; }
    [[nodiscard]] float Neg_atPCA_Px() const { return (*Source->Neg_atPCA.Px)[Entry]; }
    [[nodiscard]] float Neg_atPCA_Py() const { return (*Source->Neg_atPCA.Py)[Entry]; }
    [[nodiscard]] float Neg_atPCA_Pz() const { return (*Source->Neg_atPCA.Pz)[Entry]; }

    [[nodiscard]] float Pos_atPCA_X() const { return (*Source->Pos_atPCA.X)[Entry]; }
    [[nodiscard]] float Pos_atPCA_Y() const { return (*Source->Pos_atPCA.Y)[Entry]; }
    [[nodiscard]] float Pos_atPCA_Z() const { return (*Source->Pos_atPCA.Z)[Entry]; }
    [[nodiscard]] float Pos_atPCA_Px() const { return (*Source->Pos_atPCA.Px)[Entry]; }
    [[nodiscard]] float Pos_atPCA_Py() const { return (*Source->Pos_atPCA.Py)[Entry]; }
    [[nodiscard]] float Pos_atPCA_Pz() const { return (*Source->Pos_atPCA.Pz)[Entry]; }

    [[nodiscard]] float Chi2NDF() const { return (*Source->Chi2NDF)[Entry]; }
};

}  // namespace Tree2Secondaries::View::Rec
