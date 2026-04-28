#pragma once

#include <cstddef>
#include <span>

#include "Storage/BaseStorage.hxx"
#include "Storage/Schema/SchemaVector.hxx"
#include "View/BaseView.hxx"

namespace Tree2Secondaries::View::Rec {

struct Track : View::Base<Schema::Vector::Tracks, unsigned int> {

    Track() = delete;
    Track(const Schema::Vector::Tracks* df, unsigned int entry)  //
        : View::Base<Schema::Vector::Tracks, unsigned int>{df, entry} {}

    [[nodiscard]] float X() const { return (*Source->state.x)[Entry]; }
    [[nodiscard]] float Y() const { return (*Source->state.y)[Entry]; }
    [[nodiscard]] float Z() const { return (*Source->state.z)[Entry]; }
    [[nodiscard]] float Px() const { return (*Source->state.px)[Entry]; }
    [[nodiscard]] float Py() const { return (*Source->state.py)[Entry]; }
    [[nodiscard]] float Pz() const { return (*Source->state.pz)[Entry]; }

    [[nodiscard]] std::span<const float, Storage::VecCovMatrix<6>::n_elements> Cov() const {
        constexpr std::size_t N = Storage::VecCovMatrix<6>::n_elements;
        return std::span<const float, N>{Source->cov.mat->data() + N * Entry, N};
    }

    void AppendCov(std::vector<float>& out) const {
        const auto cov = Cov();
        out.insert(out.end(), cov.begin(), cov.end());
    }

    template <typename T>
    [[nodiscard]] T Charge() const {
        return static_cast<T>((*Source->charge)[Entry]);
    }

    [[nodiscard]] float PreDCAxy() const { return (*Source->pre_dca_xy)[Entry]; }
    [[nodiscard]] float PreDCAz() const { return (*Source->pre_dca_z)[Entry]; }
    [[nodiscard]] float TPC_Signal() const { return (*Source->tpc_signal)[Entry]; }
    [[nodiscard]] float NSigmaPion() const { return (*Source->n_sigma_pion)[Entry]; }
    [[nodiscard]] float NSigmaKaon() const { return (*Source->n_sigma_kaon)[Entry]; }
    [[nodiscard]] float NSigmaProton() const { return (*Source->n_sigma_proton)[Entry]; }

    [[nodiscard]] unsigned int EsdEntry() const { return (*Source->esd_entry)[Entry]; }
    [[nodiscard]] int McEntry() const { return (*Source->mc_entry)[Entry]; }  // NOTE: only valid when analyzing MC
};

}  // namespace Tree2Secondaries::View::Rec
