#pragma once

#include "Math/Constants.hxx"
#include "Seeder/BaseSeeder.hxx"
#include "View/Reconstructed/ViewRecV0.hxx"

namespace Tree2Secondaries::Seeder::LineVertex {

struct alignas(T2S_SIMD_ALIGN) Cache {
    // filled @ `FastPCA` //
    double px0{}, py0{}, pz0{};
    double dx{}, dy{}, dz{};
    double p2{};
    double a{};
};

// Main Method //

Seed FastPCA(const View::Rec::V0& n, const std::array<double, 3>& v, Cache* cache = nullptr);
Deriv ComputeDerivatives(const Cache& c);

// Inline Method //

inline Result FullPCA(const View::Rec::V0& n, const std::array<double, 3>& v) {
    Cache cache;
    Seed seed = FastPCA(n, v, &cache);
    Deriv deriv = ComputeDerivatives(cache);
    return {seed, deriv};
}

}  // namespace Tree2Secondaries::Seeder::LineVertex
