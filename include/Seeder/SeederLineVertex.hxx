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
    return FastPCA(v0.Decay_X, v0.Decay_Y, v0.Decay_Z, v0.Px, v0.Py, v0.Pz, vtx, cache);
}

inline Seed FastPCA(const POD::Extended::PreFoundLambda& lambda, const std::array<double, 3>& vtx, Cache* cache = nullptr) {
    return FastPCA(lambda.Decay_X, lambda.Decay_Y, lambda.Decay_Z, lambda.Px, lambda.Py, lambda.Pz, vtx, cache);
}

inline Result FullPCA(const POD::V0& v0, const std::array<double, 3>& vtx) {
    Cache cache;
    Seed seed = FastPCA(v0, vtx, &cache);
    Deriv deriv = ComputeDerivatives(cache);
    return {seed, deriv};
}

inline Result FullPCA(const POD::Extended::PreFoundLambda& lambda, const std::array<double, 3>& vtx) {
    Cache cache;
    Seed seed = FastPCA(lambda, vtx, &cache);
    Deriv deriv = ComputeDerivatives(cache);
    return {seed, deriv};
}

}  // namespace T2DS::Seeder::LineVertex
