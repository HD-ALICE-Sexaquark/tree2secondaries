#pragma once

#include "Math/Constants.hxx"
#include "Storage/Schema/SchemaVector.hxx"
#include "View/BaseView.hxx"

namespace Tree2Secondaries::View::MC {

// NOTE: needs to be guarded with `View::IsValid()` after construction.
struct Particle : View::Base<Schema::Vector::MCParticles, int> {
    Particle() = delete;
    Particle(const Schema::Vector::MCParticles* df, int entry)  //
        : View::Base<Schema::Vector::MCParticles, int>{df, entry} {}

    [[nodiscard]] float Origin_X() const { return (*Source->origin.x)[EntryAsSize()]; }
    [[nodiscard]] float Origin_Y() const { return (*Source->origin.y)[EntryAsSize()]; }
    [[nodiscard]] float Origin_Z() const { return (*Source->origin.z)[EntryAsSize()]; }
    [[nodiscard]] float Px() const { return (*Source->lv.px)[EntryAsSize()]; }
    [[nodiscard]] float Py() const { return (*Source->lv.py)[EntryAsSize()]; }
    [[nodiscard]] float Pz() const { return (*Source->lv.pz)[EntryAsSize()]; }
    [[nodiscard]] float Energy() const { return (*Source->lv.energy)[EntryAsSize()]; }

    [[nodiscard]] int PdgCode() const { return (*Source->pdg_code)[EntryAsSize()]; }
    [[nodiscard]] int Mother_McEntry() const { return (*Source->mother_mc_entry)[EntryAsSize()]; }
    [[nodiscard]] int Status() const { return (*Source->status)[EntryAsSize()]; }
    [[nodiscard]] int Generator() const { return (*Source->generator)[EntryAsSize()]; }
    [[nodiscard]] bool IsPhysPrimary() const { return (*Source->is_phys_primary)[EntryAsSize()]; }
    [[nodiscard]] bool IsSecFromMat() const { return (*Source->is_sec_from_mat)[EntryAsSize()]; }
    [[nodiscard]] bool IsSecFromWeak() const { return (*Source->is_sec_from_weak)[EntryAsSize()]; }
};

// NOTE: needs to be guarded with `View::IsValid()` after construction.
struct V0 : View::MC::Particle {
    V0() = delete;
    V0(const View::MC::Particle& neg, const View::MC::Particle& pos)
        : View::MC::Particle{neg.Source, Const::DummyInt},  // overriden in definition
          Neg{neg},
          Pos{pos} {
        if (View::IsValid(neg) && View::IsValid(pos)) {
            if (Neg.Mother_McEntry() == Pos.Mother_McEntry()) Entry = Neg.Mother_McEntry();
        }
    }

    View::MC::Particle Neg;
    View::MC::Particle Pos;
};

}  // namespace Tree2Secondaries::View::MC
