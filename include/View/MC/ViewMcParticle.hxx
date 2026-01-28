#pragma once

#include <cmath>

#include "Math/Constants.hxx"
#include "Storage/Vector/VectorMcParticles.hxx"
#include "View/BaseView.hxx"

namespace Tree2Secondaries::View::MC {

struct Particle : View::Base<Storage::Vector::MCParticles> {

    Particle() = delete;
    Particle(const Storage::Vector::MCParticles* df, int entry)  //
        : View::Base<Storage::Vector::MCParticles>{.Source = df, .Entry = entry} {}

    [[nodiscard]] float X() const { return IsValid() ? Source->X->at(Entry) : Const::DummyFloat; }
    [[nodiscard]] float Y() const { return IsValid() ? Source->Y->at(Entry) : Const::DummyFloat; }
    [[nodiscard]] float Z() const { return IsValid() ? Source->Z->at(Entry) : Const::DummyFloat; }
    [[nodiscard]] float Px() const { return IsValid() ? Source->Px->at(Entry) : Const::DummyFloat; }
    [[nodiscard]] float Py() const { return IsValid() ? Source->Py->at(Entry) : Const::DummyFloat; }
    [[nodiscard]] float Pz() const { return IsValid() ? Source->Pz->at(Entry) : Const::DummyFloat; }
    [[nodiscard]] float Energy() const { return IsValid() ? Source->Energy->at(Entry) : Const::DummyFloat; }

    [[nodiscard]] int PdgCode() const { return IsValid() ? Source->PdgCode->at(Entry) : Const::DummyInt; }
    [[nodiscard]] int MotherEntry() const { return IsValid() ? Source->MotherEntry->at(Entry) : Const::DummyInt; }
    [[nodiscard]] int Generator() const { return IsValid() ? Source->Generator->at(Entry) : Const::DummyInt; }
    [[nodiscard]] int Status() const { return IsValid() ? Source->Status->at(Entry) : Const::DummyInt; }
    [[nodiscard]] int IsSecFromMat() const { return IsValid() ? Source->IsSecFromMat->at(Entry) : Const::DummyInt; }
    [[nodiscard]] int IsSecFromWeak() const { return IsValid() ? Source->IsSecFromWeak->at(Entry) : Const::DummyInt; }
};

struct V0 : View::MC::Particle {

    View::MC::Particle Neg;
    View::MC::Particle Pos;

    V0() = delete;
    V0(const Storage::Vector::MCParticles* df, int neg_entry, int pos_entry)
        : View::MC::Particle{df, Const::DummyInt},  // overriden in definition
          Neg{df, neg_entry},
          Pos{df, pos_entry} {
        if (Neg.MotherEntry() == Pos.MotherEntry()) Entry = Neg.MotherEntry();
    }
};

}  // namespace Tree2Secondaries::View::MC
