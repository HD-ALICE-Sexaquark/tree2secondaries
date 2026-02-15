#pragma once

#include <cstdlib>

#include "Math/Constants.hxx"
#include "Storage/Vector/VectorMcParticles.hxx"
#include "View/BaseView.hxx"

namespace Tree2Secondaries::View::MC {

// NOTE: need to be guarded with `View::IsValid()` after construction.
struct Particle : View::Base<Storage::Vector::MCParticles, int> {

    Particle() = delete;
    Particle(const Storage::Vector::MCParticles* df, int entry)  //
        : View::Base<Storage::Vector::MCParticles, int>{df, entry} {}

    [[nodiscard]] float X() const { return (*Source->X)[EntryAsSize()]; }
    [[nodiscard]] float Y() const { return (*Source->Y)[EntryAsSize()]; }
    [[nodiscard]] float Z() const { return (*Source->Z)[EntryAsSize()]; }
    [[nodiscard]] float Px() const { return (*Source->Px)[EntryAsSize()]; }
    [[nodiscard]] float Py() const { return (*Source->Py)[EntryAsSize()]; }
    [[nodiscard]] float Pz() const { return (*Source->Pz)[EntryAsSize()]; }
    [[nodiscard]] float Energy() const { return (*Source->Energy)[EntryAsSize()]; }

    [[nodiscard]] int PdgCode() const { return (*Source->PdgCode)[EntryAsSize()]; }
    [[nodiscard]] int MotherMcEntry() const { return (*Source->MotherMcEntry)[EntryAsSize()]; }
    [[nodiscard]] int Status() const { return (*Source->Status)[EntryAsSize()]; }
    [[nodiscard]] int Generator() const { return (*Source->Generator)[EntryAsSize()]; }
    [[nodiscard]] bool IsPrimary() const { return (*Source->IsPrimary)[EntryAsSize()]; }
    [[nodiscard]] bool IsSecFromMat() const { return (*Source->IsSecFromMat)[EntryAsSize()]; }
    [[nodiscard]] bool IsSecFromWeak() const { return (*Source->IsSecFromWeak)[EntryAsSize()]; }
};

// NOTE: need to be guarded with `View::IsValid()` after construction.
struct V0 : View::MC::Particle {

    V0() = delete;
    V0(const View::MC::Particle& neg, const View::MC::Particle& pos)
        : View::MC::Particle{neg.Source, Const::DummyInt},  // overriden in definition
          Neg{neg},
          Pos{pos} {
        if (View::IsValid(neg) && View::IsValid(pos)) {
            if (Neg.MotherMcEntry() == Pos.MotherMcEntry()) Entry = Neg.MotherMcEntry();
        }
    }

    View::MC::Particle Neg;
    View::MC::Particle Pos;
};

}  // namespace Tree2Secondaries::View::MC
