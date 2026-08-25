#pragma once

#include <utility>

#include "Seeder/SeederTypes.hxx"

namespace T2DS::Seeder::LineLine {

struct Cache {
    // filled @ `FastPCAs` //
    double px01{}, py01{}, pz01{};
    double px02{}, py02{}, pz02{};
    double dx{}, dy{}, dz{};
    double p12{};
    double p22{};
    double p1p2{};
    double drp1{};
    double drp2{};
    double detp{};
};

// Main Methods //

SeedsPair FastPCAs(const State& s01, const State& s02, Cache& c);
std::pair<Deriv, Deriv> ComputeDerivatives(const Cache& c);

// Inline Methods //

inline SeedsPair FastPCAs(const POD::V0& v01, const POD::V0& v02, Cache& c) {  //
    return FastPCAs(State::FromV0(v01), State::FromV0(v02), c);
}
inline SeedsPair FastPCAs(const POD::Extended::PreFoundLambda& l1, const POD::Extended::PreFoundLambda& l2, Cache& c) {
    return FastPCAs(State::FromPreFoundLambda(l1), State::FromPreFoundLambda(l2), c);
}

inline std::pair<Result, Result> FullPCAs(const State& s01, const State& s02) {
    Cache c;
    auto [seed1, seed2] = FastPCAs(s01, s02, c);
    auto [d1, d2] = ComputeDerivatives(c);
    return {{seed1, d1}, {seed2, d2}};
}
inline std::pair<Result, Result> FullPCAs(const POD::V0& v01, const POD::V0& v02) {  //
    return FullPCAs(State::FromV0(v01), State::FromV0(v02));
}
inline std::pair<Result, Result> FullPCAs(const POD::Extended::PreFoundLambda& l1, const POD::Extended::PreFoundLambda& l2) {
    return FullPCAs(State::FromPreFoundLambda(l1), State::FromPreFoundLambda(l2));
}

}  // namespace T2DS::Seeder::LineLine
