#pragma once

#include <cmath>
#include <cstddef>

#include "Schema/SchemaVectorInjected.hxx"

namespace Tree2Secondaries::View {

struct VecInjected {

    VecInjected() = delete;
    VecInjected(const Schema::VecInjected* df, std::size_t idx) : source{df}, index{idx} {}

    // ## Getters ## //

    [[nodiscard]] float SV_X() const { return (*source->sv.x)[index]; }
    [[nodiscard]] float SV_Y() const { return (*source->sv.y)[index]; }
    [[nodiscard]] float SV_Z() const { return (*source->sv.z)[index]; }
    [[nodiscard]] float Px() const { return (*source->mom.px)[index]; }
    [[nodiscard]] float Py() const { return (*source->mom.py)[index]; }
    [[nodiscard]] float Pz() const { return (*source->mom.pz)[index]; }
    [[nodiscard]] float Nucleon_Px() const { return (*source->mom_nucleon.px)[index]; }
    [[nodiscard]] float Nucleon_Py() const { return (*source->mom_nucleon.py)[index]; }
    [[nodiscard]] float Nucleon_Pz() const { return (*source->mom_nucleon.pz)[index]; }
    [[nodiscard]] unsigned int ReactionID() const { return (*source->reaction_id)[index]; }

    // ## Operations ## //

    [[nodiscard]] float Energy(double mass) const {
        auto px = static_cast<double>(Px());
        auto py = static_cast<double>(Py());
        auto pz = static_cast<double>(Pz());
        return static_cast<float>(std::sqrt(mass * mass + px * px + py * py + pz * pz));
    }

    [[nodiscard]] float Nucleon_Energy(double n_mass) const {
        auto px = static_cast<double>(Nucleon_Px());
        auto py = static_cast<double>(Nucleon_Py());
        auto pz = static_cast<double>(Nucleon_Pz());
        return static_cast<float>(std::sqrt(n_mass * n_mass + px * px + py * py + pz * pz));
    }

    // ## Member Variables ## //

    const Schema::VecInjected* source{};
    std::size_t index{};
};

}  // namespace Tree2Secondaries::View
