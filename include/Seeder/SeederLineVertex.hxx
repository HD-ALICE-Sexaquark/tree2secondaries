#pragma once

#include <array>

#include "common/POD_PreFoundLambda.hpp"
#include "common/POD_V0.hpp"

#include "Seeder/BaseSeeder.hxx"

namespace T2DS::Seeder::LineVertex {

struct Cache {
    // filled @ `FastPCA` //
    double px0{}, py0{}, pz0{};
    double dx{}, dy{}, dz{};
    double p2{};
    double a{};
};

// Main Methods //

Seed FastPCA(double x0, double y0, double z0, double px, double py, double pz,  //
             const std::array<double, 3>& vtx, Cache* cache = nullptr);

Deriv ComputeDerivatives(const Cache& c);

// Inline Methods //

inline Seed FastPCA(const POD::V0& v0, const std::array<double, 3>& vtx, Cache* cache = nullptr) {
    return FastPCA(static_cast<double>(v0.Decay_X), static_cast<double>(v0.Decay_Y), static_cast<double>(v0.Decay_Z), static_cast<double>(v0.Px),
                   static_cast<double>(v0.Py), static_cast<double>(v0.Pz), vtx, cache);
}

inline Seed FastPCA(const POD::Extended::PreFoundLambda& lambda, const std::array<double, 3>& vtx, Cache* cache = nullptr) {
    return FastPCA(static_cast<double>(lambda.Decay_X), static_cast<double>(lambda.Decay_Y), static_cast<double>(lambda.Decay_Z),
                   static_cast<double>(lambda.Px), static_cast<double>(lambda.Py), static_cast<double>(lambda.Pz), vtx, cache);
}

inline Result FullPCA(const POD::V0& v0, const std::array<double, 3>& vtx) {
    Cache cache;
    Seed seed = FastPCA(v0, vtx, &cache);
    return {seed, ComputeDerivatives(cache)};
}

inline Result FullPCA(const POD::Extended::PreFoundLambda& lambda, const std::array<double, 3>& vtx) {
    Cache cache;
    Seed seed = FastPCA(lambda, vtx, &cache);
    return {seed, ComputeDerivatives(cache)};
}

inline Result FullPCA(double x0, double y0, double z0, double px, double py, double pz, const std::array<double, 3>& vtx) {
    Cache cache;
    Seed seed = FastPCA(x0, y0, z0, px, py, pz, vtx, &cache);
    return {seed, ComputeDerivatives(cache)};
}

}  // namespace T2DS::Seeder::LineVertex
