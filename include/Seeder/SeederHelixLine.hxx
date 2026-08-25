#pragma once

#include <utility>

#include "Seeder/SeederTypes.hxx"

namespace T2DS::Seeder::HelixLine {

struct Cache {
    // filled @ `FastPCAs_XY` //
    double bq1{};
    double x01{}, y01{}, z01{};
    double px01{}, py01{}, pz01{};
    double pt12{};
    double x02{}, y02{}, z02{};
    double px02{}, py02{}, pz02{};
    double pt22{};
    double dx0{}, dy0{};
    double drp2{};
    double dxyp2{};
    double p1p2{};
    double dp1p2{};
    double k11{}, k21{}, k12{}, k22{};
    double kp{}, kd{};
    double c1{}, c2{};
    double d1{};
    double dx{}, dy{}, dz{};
    // filled @ `CorrectPCAs_Z` //
    double px1{}, py1{};
    double p12{}, p22{}, lp1p2{}, detp{};
    double ldrp1{}, ldrp2{};
    double a1{}, a2{};
    // filled @ `ComputeDerivatives_XY` //
    double aa1{}, bb1{}, cc1{}, dd1{};
    double na{}, nb{}, nb2{};
    // filled @ `FastPCAs_XY`, but here for memory alignment reasons //
    int w_sign{};
    // status flag //
    bool pca_dz_worked{false};
};

// Main Methods //

SeedsPair FastPCAs_XY(const State& s01, const State& s02, double bz, Cache& c);
SeedsPair CorrectPCAs_Z(const SeedsPair& seeds_xy, Cache& c);

std::pair<Deriv, Deriv> ComputeDerivatives_XY(Cache& c);
std::pair<Deriv, Deriv> UpdateDerivatives_Z(const SeedsPair& seeds_xy, const Deriv& d1_xy, const Deriv& d2_xy, const Cache& c);

// Inline Methods //

inline SeedsPair FastPCAs(const State& s01, const State& s02, double bz, Cache& c, SeedsPair& seeds_xy) {
    seeds_xy = FastPCAs_XY(s01, s02, bz, c);
    return CorrectPCAs_Z(seeds_xy, c);
}
inline SeedsPair FastPCAs(const POD::Track& q1, const POD::V0& n2, double bz, Cache& c, SeedsPair& seeds_xy) {
    return FastPCAs(State::FromTrack(q1), State::FromV0(n2), bz, c, seeds_xy);
}

inline std::pair<Deriv, Deriv> ComputeDerivatives(const SeedsPair& seeds_xy, Cache& c) {
    auto [d1_xy, d2_xy] = ComputeDerivatives_XY(c);
    return UpdateDerivatives_Z(seeds_xy, d1_xy, d2_xy, c);
}

inline std::pair<Result, Result> FullPCAs(const State& s01, const State& s02, double bz) {
    Cache c;
    SeedsPair seeds_xy;
    auto [s1, s2] = FastPCAs(s01, s02, bz, c, seeds_xy);
    auto [d1, d2] = ComputeDerivatives(seeds_xy, c);
    return {{s1, d1}, {s2, d2}};
}
inline std::pair<Result, Result> FullPCAs(const POD::Track& q1, const POD::V0& n2, double bz) {
    return FullPCAs(State::FromTrack(q1), State::FromV0(n2), bz);
}

}  // namespace T2DS::Seeder::HelixLine
