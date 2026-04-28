#pragma once

#include <cstddef>
#include <span>

#include "Storage/Schema/SchemaVector.hxx"
#include "View/BaseView.hxx"
#include "View/Reconstructed/ViewRecTrack.hxx"

namespace Tree2Secondaries::View::Rec {

struct V0 : View::Base<Schema::Vector::V0s, unsigned int> {

    V0() = delete;
    V0(const Schema::Vector::V0s* df, unsigned int entry)
        : View::Base<Schema::Vector::V0s, unsigned int>{df, entry},
          Neg{&df->neg, entry},  // NOTE: yes, same `entry`
          Pos{&df->pos, entry}   // NOTE: yes, same `entry`
    {}

    View::Rec::Track Neg;
    View::Rec::Track Pos;

    [[nodiscard]] float X() const { return (*Source->decay.x)[Entry]; }
    [[nodiscard]] float Y() const { return (*Source->decay.y)[Entry]; }
    [[nodiscard]] float Z() const { return (*Source->decay.z)[Entry]; }
    [[nodiscard]] float Px() const { return (*Source->lv.px)[Entry]; }
    [[nodiscard]] float Py() const { return (*Source->lv.py)[Entry]; }
    [[nodiscard]] float Pz() const { return (*Source->lv.pz)[Entry]; }
    [[nodiscard]] float Energy() const { return (*Source->lv.energy)[Entry]; }

    [[nodiscard]] std::span<const float, Storage::VecCovMatrix<7>::n_elements> Cov() const {
        constexpr std::size_t N = Storage::VecCovMatrix<7>::n_elements;
        return std::span<const float, N>{Source->cov.mat->data() + N * Entry, N};
    }

    [[nodiscard]] float Neg_PCAwrtV0_X() const { return (*Source->neg_pca_v0.x)[Entry]; }
    [[nodiscard]] float Neg_PCAwrtV0_Y() const { return (*Source->neg_pca_v0.y)[Entry]; }
    [[nodiscard]] float Neg_PCAwrtV0_Z() const { return (*Source->neg_pca_v0.z)[Entry]; }
    [[nodiscard]] float Neg_PCAwrtV0_Px() const { return (*Source->neg_pca_v0.px)[Entry]; }
    [[nodiscard]] float Neg_PCAwrtV0_Py() const { return (*Source->neg_pca_v0.py)[Entry]; }
    [[nodiscard]] float Neg_PCAwrtV0_Pz() const { return (*Source->neg_pca_v0.pz)[Entry]; }

    [[nodiscard]] float Pos_PCAwrtV0_X() const { return (*Source->pos_pca_v0.x)[Entry]; }
    [[nodiscard]] float Pos_PCAwrtV0_Y() const { return (*Source->pos_pca_v0.y)[Entry]; }
    [[nodiscard]] float Pos_PCAwrtV0_Z() const { return (*Source->pos_pca_v0.z)[Entry]; }
    [[nodiscard]] float Pos_PCAwrtV0_Px() const { return (*Source->pos_pca_v0.px)[Entry]; }
    [[nodiscard]] float Pos_PCAwrtV0_Py() const { return (*Source->pos_pca_v0.py)[Entry]; }
    [[nodiscard]] float Pos_PCAwrtV0_Pz() const { return (*Source->pos_pca_v0.pz)[Entry]; }

    [[nodiscard]] float Chi2NDF() const { return (*Source->chi2ndf)[Entry]; }
};

}  // namespace Tree2Secondaries::View::Rec
