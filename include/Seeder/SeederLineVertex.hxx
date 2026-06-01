#pragma once

#include <array>

#include "common/POD_OnTheFlyLambda.hpp"
#include "common/POD_V0.hpp"

#include "Seeder/BaseSeeder.hxx"

namespace R2DS::Seeder::LineVertex {

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

inline Result FullPCA(const POD::V0& v0, const std::array<double, 3>& vtx) {
    Cache cache;
    Seed seed = FastPCA(v0, vtx, &cache);
    Deriv deriv = ComputeDerivatives(cache);
    return {seed, deriv};
}

inline Seed FullPCA(const POD::OnTheFlyLambda& l, const std::array<double, 3>& vtx) {
    return FastPCA(l.Decay_X, l.Decay_Y, l.Decay_Z, l.Px, l.Py, l.Pz, vtx);
}

}  // namespace R2DS::Seeder::LineVertex
