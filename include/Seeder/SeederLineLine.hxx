#pragma once

#include "Math/Constants.hxx"
#include "Seeder/BaseSeeder.hxx"
#include "View/ViewVectorV0s.hxx"

namespace Tree2Secondaries::Seeder::LineLine {

struct alignas(T2S_SIMD_ALIGN) Cache {
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

std::pair<Seed, Seed> FastPCAs(const View::VecV0s& n1, const View::VecV0s& n2, Cache* cache = nullptr);
std::pair<Deriv, Deriv> ComputeDerivatives(const Cache& c);

// Inline Method //

inline std::tuple<Result, Result> FullPCAs(const View::VecV0s& n1, const View::VecV0s& n2) {
    Cache cache;
    auto [seed1, seed2] = FastPCAs(n1, n2, &cache);
    auto [deriv1, deriv2] = ComputeDerivatives(cache);
    return {{seed1, deriv1}, {seed2, deriv2}};
}

}  // namespace Tree2Secondaries::Seeder::LineLine
