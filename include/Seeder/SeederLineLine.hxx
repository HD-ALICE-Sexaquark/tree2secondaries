#pragma once

#include "common/POD_OnTheFlyLambda.hpp"
#include "common/POD_V0.hpp"

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

std::pair<Seed, Seed> FastPCAs(double x01, double y01, double z01, double px01, double py01, double pz01,  //
                               double x02, double y02, double z02, double px02, double py02, double pz02,  //
                               Cache* cache = nullptr);

std::pair<Deriv, Deriv> ComputeDerivatives(const Cache& c);

// Inline Methods //

inline std::pair<Seed, Seed> FastPCAs(const POD::V0& v01, const POD::V0& v02, Cache* cache = nullptr) {
    return FastPCAs(v01.Decay_X, v01.Decay_Y, v01.Decay_Z, v01.Px, v01.Py, v01.Pz,  //
                    v02.Decay_X, v02.Decay_Y, v02.Decay_Z, v02.Px, v02.Py, v02.Pz,  //
                    cache);
}

inline std::pair<Result, Result> FullPCAs(const POD::V0& v01, const POD::V0& v02) {
    Cache cache;
    auto [seed1, seed2] = FastPCAs(v01, v02, &cache);
    auto [deriv1, deriv2] = ComputeDerivatives(cache);
    return {{seed1, deriv1}, {seed2, deriv2}};
}

inline std::pair<Seed, Seed> FullPCAs(const POD::OnTheFlyLambda& l1, const POD::OnTheFlyLambda& l2) {
    return FastPCAs(l1.Decay_X, l1.Decay_Y, l1.Decay_Z, l1.Px, l1.Py, l1.Pz,  //
                    l2.Decay_X, l2.Decay_Y, l2.Decay_Z, l2.Px, l2.Py, l2.Pz,  //
                    nullptr);
}

}  // namespace R2DS::Seeder::LineLine
