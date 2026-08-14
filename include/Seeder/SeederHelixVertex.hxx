#pragma once

#include <array>
#include <tuple>

#include "common/POD_Track.hpp"

#include "Seeder/BaseSeeder.hxx"

namespace T2DS::Seeder::HelixVertex {

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

Seed FastPCA_XY(double x0, double y0, double z0, double px, double py, double pz, int charge,  //
                const std::array<double, 3>& v, double bz, Cache* cache = nullptr);
Seed CorrectPCA_Z(const Seed& s_xy, Cache& c);

Deriv ComputeDerivatives_XY(Cache& c);
Deriv UpdateDerivatives_Z(const Seed& s_xy, const Deriv& d_xy, const Cache& c);

// Inline Methods //

inline Seed FastPCA_XY(const POD::Track& q, const std::array<double, 3>& v, double bz, Cache* cache = nullptr) {
    return FastPCA_XY(static_cast<double>(q.X), static_cast<double>(q.Y), static_cast<double>(q.Z), static_cast<double>(q.Px),
                      static_cast<double>(q.Py), static_cast<double>(q.Pz), q.Charge, v, bz, cache);
}

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

inline Result FullPCA(double x0, double y0, double z0, double px, double py, double pz, int charge, const std::array<double, 3>& xyz, double bz) {
    Cache cache;
    Seed seed_xy = FastPCA_XY(x0, y0, z0, px, py, pz, charge, xyz, bz, &cache);
    Seed seed = CorrectPCA_Z(seed_xy, cache);
    Deriv deriv_xy = ComputeDerivatives_XY(cache);
    return {seed, UpdateDerivatives_Z(seed_xy, deriv_xy, cache)};
}

}  // namespace T2DS::Seeder::HelixVertex
