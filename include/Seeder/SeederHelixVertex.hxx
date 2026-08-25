#pragma once

#include <array>

#include "Seeder/SeederTypes.hxx"

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
    bool pca_dz_worked{false};
};

// Main Methods //

Seed FastPCA_XY(const State& s0, const std::array<double, 3>& vtx, double bz, Cache& c);
Seed CorrectPCA_Z(const Seed& s_xy, Cache& c);

Deriv ComputeDerivatives_XY(Cache& c);
Deriv UpdateDerivatives_Z(const Seed& s_xy, const Deriv& d_xy, const Cache& c);

// Inline Methods //

inline Seed FastPCA(const State& s0, const std::array<double, 3>& vtx, double bz, Cache& c, Seed& s_xy) {
    s_xy = FastPCA_XY(s0, vtx, bz, c);
    return CorrectPCA_Z(s_xy, c);
}
inline Seed FastPCA(const POD::Track& q, const std::array<double, 3>& vtx, double bz, Cache& c, Seed& s_xy) {  //
    return FastPCA(State::FromTrack(q), vtx, bz, c, s_xy);
}

inline Deriv ComputeDerivatives(const Seed& s_xy, Cache& c) {
    auto d_xy = ComputeDerivatives_XY(c);
    return UpdateDerivatives_Z(s_xy, d_xy, c);
}

inline Result FullPCA(const State& s0, const std::array<double, 3>& vtx, double bz) {
    Cache c;
    Seed s_xy;
    Seed seed = FastPCA(s0, vtx, bz, c, s_xy);
    return {seed, ComputeDerivatives(s_xy, c)};
}
inline Result FullPCA(const POD::Track& q, const std::array<double, 3>& vtx, double bz) {  //
    return FullPCA(State::FromTrack(q), vtx, bz);
}

}  // namespace T2DS::Seeder::HelixVertex
