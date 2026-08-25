#pragma once

#include <array>

#include "Seeder/SeederTypes.hxx"

namespace T2DS::Seeder::LineVertex {

struct Cache {
    // filled @ `FastPCA` //
    double px0{}, py0{}, pz0{};
    double dx{}, dy{}, dz{};
    double p2{};
    double a{};
};

// Main Methods //

Seed FastPCA(const State& s0, const std::array<double, 3>& vtx, Cache& c);
Deriv ComputeDerivatives(const Cache& c);

// Inline Methods //

inline Seed FastPCA(const POD::V0& v0, const std::array<double, 3>& vtx, Cache& c) {  //
    return FastPCA(State::FromV0(v0), vtx, c);
}
inline Seed FastPCA(const POD::Extended::PreFoundLambda& l, const std::array<double, 3>& vtx, Cache& c) {
    return FastPCA(State::FromPreFoundLambda(l), vtx, c);
}

inline Result FullPCA(const State& s0, const std::array<double, 3>& vtx) {
    Cache c;
    Seed seed = FastPCA(s0, vtx, c);
    return {seed, ComputeDerivatives(c)};
}
inline Result FullPCA(const POD::V0& v0, const std::array<double, 3>& vtx) {  //
    return FullPCA(State::FromV0(v0), vtx);
}
inline Result FullPCA(const POD::Extended::PreFoundLambda& l, const std::array<double, 3>& vtx) {  //
    return FullPCA(State::FromPreFoundLambda(l), vtx);
}

}  // namespace T2DS::Seeder::LineVertex
