#pragma once

#include <array>

#include "common/VC_OnTheFlyLambdaView.hpp"
#include "common/VC_V0View.hpp"

#include "Seeder/BaseSeeder.hxx"

namespace R2DS::Seeder::LineVertex {

struct Cache {
    // filled @ `FastPCA` //
    double px0{}, py0{}, pz0{};
    double dx{}, dy{}, dz{};
    double p2{};
    double a{};
};

// Main Method //

Seed FastPCA(const Vector::V0View& n, const std::array<double, 3>& v, Cache* cache = nullptr);
Deriv ComputeDerivatives(const Cache& c);

// Inline Method //

inline Result FullPCA(const Vector::V0View& n, const std::array<double, 3>& v) {
    Cache cache;
    Seed seed = FastPCA(n, v, &cache);
    Deriv deriv = ComputeDerivatives(cache);
    return {seed, deriv};
}

// Alternative Method //

Seed FullPCA(const Vector::OnTheFlyLambdaView& l, const std::array<double, 3>& v);

}  // namespace R2DS::Seeder::LineVertex
