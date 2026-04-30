#pragma once

#include <cstddef>

#include "Schema/SchemaVectorV0s.hxx"
#include "View/ViewVectorTracks.hxx"

namespace Tree2Secondaries::View {

struct VecV0s {

    VecV0s() = delete;
    VecV0s(const Schema::VecV0s* df, std::size_t idx)  //
        : source{df},                                  //
          neg{&df->neg, idx},
          pos{&df->pos, idx},
          index{idx} {}

    // ## Getters ## //

    [[nodiscard]] float X() const { return (*source->decay.x)[index]; }
    [[nodiscard]] float Y() const { return (*source->decay.y)[index]; }
    [[nodiscard]] float Z() const { return (*source->decay.z)[index]; }
    [[nodiscard]] float Px() const { return (*source->lv.px)[index]; }
    [[nodiscard]] float Py() const { return (*source->lv.py)[index]; }
    [[nodiscard]] float Pz() const { return (*source->lv.pz)[index]; }
    [[nodiscard]] float Energy() const { return (*source->lv.energy)[index]; }
    [[nodiscard]] float Chi2NDF() const { return (*source->chi2ndf)[index]; }

    [[nodiscard]] float Neg_PCAwrtV0_X() const { return (*source->neg_pca_v0.x)[index]; }
    [[nodiscard]] float Neg_PCAwrtV0_Y() const { return (*source->neg_pca_v0.y)[index]; }
    [[nodiscard]] float Neg_PCAwrtV0_Z() const { return (*source->neg_pca_v0.z)[index]; }
    [[nodiscard]] float Neg_PCAwrtV0_Px() const { return (*source->neg_pca_v0.px)[index]; }
    [[nodiscard]] float Neg_PCAwrtV0_Py() const { return (*source->neg_pca_v0.py)[index]; }
    [[nodiscard]] float Neg_PCAwrtV0_Pz() const { return (*source->neg_pca_v0.pz)[index]; }

    [[nodiscard]] float Pos_PCAwrtV0_X() const { return (*source->pos_pca_v0.x)[index]; }
    [[nodiscard]] float Pos_PCAwrtV0_Y() const { return (*source->pos_pca_v0.y)[index]; }
    [[nodiscard]] float Pos_PCAwrtV0_Z() const { return (*source->pos_pca_v0.z)[index]; }
    [[nodiscard]] float Pos_PCAwrtV0_Px() const { return (*source->pos_pca_v0.px)[index]; }
    [[nodiscard]] float Pos_PCAwrtV0_Py() const { return (*source->pos_pca_v0.py)[index]; }
    [[nodiscard]] float Pos_PCAwrtV0_Pz() const { return (*source->pos_pca_v0.pz)[index]; }

    // ## Operations ## //

    [[nodiscard]] double Cov(std::size_t i, std::size_t j) const {
        return static_cast<double>((*source->cov.mat)[28 * index + 7 * i + j - (i * (i + 1)) / 2]);
    }
    [[nodiscard]] std::array<float, 28> Cov() const {
        std::array<float, 28> cov{};
        for (std::size_t k = 0; k < 28; ++k) {
            cov[k] = static_cast<float>((*source->cov.mat)[28 * index + k]);
        }
        return cov;
    }

    // ## Member Variables ## //

    const Schema::VecV0s* source{};
    View::VecTracks neg;
    View::VecTracks pos;
    std::size_t index{};
};

}  // namespace Tree2Secondaries::View
