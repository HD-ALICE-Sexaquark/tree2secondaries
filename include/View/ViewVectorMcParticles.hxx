#pragma once

#include <cstddef>

#include "Schema/SchemaVectorMcParticles.hxx"

namespace Tree2Secondaries::View {

struct VecMcParticles {

    VecMcParticles() = delete;
    VecMcParticles(const Schema::VecMcParticles* df, std::size_t idx) : source{df}, index{idx} {}

    // ## Getters ## //

    [[nodiscard]] int PdgCode() const { return (*source->pdg_code)[index]; }

    template <typename T>
    [[nodiscard]] T Charge() const {
        return static_cast<T>((*source->charge)[index]);
    }

    [[nodiscard]] float Origin_X() const { return (*source->origin.x)[index]; }
    [[nodiscard]] float Origin_Y() const { return (*source->origin.y)[index]; }
    [[nodiscard]] float Origin_Z() const { return (*source->origin.z)[index]; }
    [[nodiscard]] float Px() const { return (*source->lv.px)[index]; }
    [[nodiscard]] float Py() const { return (*source->lv.py)[index]; }
    [[nodiscard]] float Pz() const { return (*source->lv.pz)[index]; }
    [[nodiscard]] float Energy() const { return (*source->lv.energy)[index]; }

    [[nodiscard]] int Mother_McEntry() const { return (*source->mother_mc_entry)[index]; }
    [[nodiscard]] unsigned int N_Daughters() const { return (*source->n_daughters)[index]; }
    [[nodiscard]] int FirstDau_McEntry() const { return (*source->firstdau_mc_entry)[index]; }
    [[nodiscard]] int LastDau_McEntry() const { return (*source->lastdau_mc_entry)[index]; }

    [[nodiscard]] unsigned int Status() const { return (*source->status)[index]; }
    [[nodiscard]] char Generator() const { return (*source->generator)[index]; }
    [[nodiscard]] bool IsPhysPrimary() const { return (*source->is_phys_primary)[index]; }
    [[nodiscard]] bool IsSecFromMat() const { return (*source->is_sec_from_mat)[index]; }
    [[nodiscard]] bool IsSecFromWeak() const { return (*source->is_sec_from_weak)[index]; }

    // ## Operations ## //

    [[nodiscard]] bool IsLogicalPrimary() const { return Mother_McEntry() == Const::DummyInt; }

    // ## Member Variables ## //

    const Schema::VecMcParticles* source{};
    std::size_t index{};
};

}  // namespace Tree2Secondaries::View
