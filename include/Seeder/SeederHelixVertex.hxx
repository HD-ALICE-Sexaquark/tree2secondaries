#pragma once

#include "common/POD_Track.hpp"

#include "Seeder/BaseSeeder.hxx"

namespace R2DS::Seeder::HelixVertex {

struct Cache {
    // filled @ `FastPCA_XY` //
    double x0{}, y0{}, z0{};
    double px0{}, py0{}, pz0{};
    double pt2{};
    double dx{}, dy{}, dz{};
    double a{};
    double bq{};
    double abq{};
    // filled @ `CorrectPCA_Z` //
    double bbq{};
    double cbq{};
    double sz{};
    // filled @ `ComputeDerivatives_XY` //
    double den{};
    // status flag //
    int pca_dz_worked{0};
};

// Main Methods //

Seed FastPCA_XY(const POD::Track& q, const std::array<double, 3>& v, double bz, Cache* cache = nullptr);
Seed CorrectPCA_Z(const Seed& s_xy, Cache& c);

Deriv ComputeDerivatives_XY(Cache& c);
Deriv UpdateDerivatives_Z(const Seed& s_xy, const Deriv& d_xy, const Cache& c);

// Inline Methods //

inline std::tuple<Seed, Cache> FastCorrectPCA(const POD::Track& q, const std::array<double, 3>& v, double bz) {
    Cache cache;
    auto seed_xy = FastPCA_XY(q, v, bz, &cache);
    auto seed = CorrectPCA_Z(seed_xy, cache);
    return {seed, cache};
}

inline Deriv ComputeDerivatives(const Seed& seed_xy, Cache& cache) {
    auto deriv_xy = ComputeDerivatives_XY(cache);
    return UpdateDerivatives_Z(seed_xy, deriv_xy, cache);
}

}  // namespace R2DS::Seeder::HelixVertex
