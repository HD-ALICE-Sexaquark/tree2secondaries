#pragma once

#include "common/VC_OnTheFlyLambdaView.hpp"
#include "common/VC_V0View.hpp"

#include "Seeder/BaseSeeder.hxx"

namespace R2DS::Seeder::LineLine {

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

std::pair<Seed, Seed> FastPCAs(const Vector::V0View& n1, const Vector::V0View& n2, Cache* cache = nullptr);
std::pair<Deriv, Deriv> ComputeDerivatives(const Cache& c);

// Inline Method //

inline std::pair<Result, Result> FullPCAs(const Vector::V0View& n1, const Vector::V0View& n2) {
    Cache cache;
    auto [seed1, seed2] = FastPCAs(n1, n2, &cache);
    auto [deriv1, deriv2] = ComputeDerivatives(cache);
    return {{seed1, deriv1}, {seed2, deriv2}};
}

// Alternative Method //

std::pair<Seed, Seed> FullPCAs(const Vector::OnTheFlyLambdaView& l1, const Vector::OnTheFlyLambdaView& l2);

}  // namespace R2DS::Seeder::LineLine
