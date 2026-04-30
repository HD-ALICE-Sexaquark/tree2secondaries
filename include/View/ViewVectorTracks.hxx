#pragma once

#include <array>
#include <cstddef>

#include "Schema/SchemaVectorTracks.hxx"

namespace Tree2Secondaries::View {

struct VecTracks {

    VecTracks() = delete;
    VecTracks(const Schema::VecTracks* df, std::size_t idx) : source{df}, index{idx} {}

    // ## Getters ## //

    [[nodiscard]] float X() const { return (*source->state.x)[index]; }
    [[nodiscard]] float Y() const { return (*source->state.y)[index]; }
    [[nodiscard]] float Z() const { return (*source->state.z)[index]; }
    [[nodiscard]] float Px() const { return (*source->state.px)[index]; }
    [[nodiscard]] float Py() const { return (*source->state.py)[index]; }
    [[nodiscard]] float Pz() const { return (*source->state.pz)[index]; }

    template <typename T>
    [[nodiscard]] T Charge() const {
        return static_cast<T>((*source->charge)[index]);
    }

    [[nodiscard]] float PreDCAxy() const { return (*source->pre_dca_xy)[index]; }
    [[nodiscard]] float PreDCAz() const { return (*source->pre_dca_z)[index]; }
    [[nodiscard]] float TPC_Signal() const { return (*source->tpc_signal)[index]; }
    [[nodiscard]] float NSigmaPion() const { return (*source->n_sigma_pion)[index]; }
    [[nodiscard]] float NSigmaKaon() const { return (*source->n_sigma_kaon)[index]; }
    [[nodiscard]] float NSigmaProton() const { return (*source->n_sigma_proton)[index]; }

    [[nodiscard]] unsigned int EsdEntry() const { return (*source->esd_entry)[index]; }
    [[nodiscard]] int McEntry() const { return (*source->mc_entry)[index]; }  // NOTE: only valid when analyzing MC

    // ## Operations ## //

    [[nodiscard]] double Cov(std::size_t i, std::size_t j) const {
        return static_cast<double>((*source->cov.mat)[21 * index + 6 * i + j - (i * (i + 1)) / 2]);
    }
    [[nodiscard]] std::array<float, 21> Cov() const {
        std::array<float, 21> cov{};
        for (std::size_t k = 0; k < 21; ++k) {
            cov[k] = static_cast<float>((*source->cov.mat)[21 * index + k]);
        }
        return cov;
    }

    // ## Member Variables ## //

    const Schema::VecTracks* source{};
    std::size_t index{};
};

}  // namespace Tree2Secondaries::View
