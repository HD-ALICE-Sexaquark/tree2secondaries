#pragma once

#include <tuple>

#include "common/POD_Track.hpp"

#include "Seeder/BaseSeeder.hxx"

namespace T2DS::Seeder::HelixHelix {

struct Cache {
    // filled @ `FastPCAs_XY` //
    double bq1{}, bq2{};
    double x01{}, y01{}, z01{};
    double px01{}, py01{}, pz01{};
    double pt12{};
    double x02{}, y02{}, z02{};
    double px02{}, py02{}, pz02{};
    double pt22{};
    double dx0{}, dy0{};
    double dr02{};
    double drp1{};
    double dxyp1{};
    double drp2{};
    double dxyp2{};
    double p1p2{};
    double dp1p2{};
    double k11{}, k21{}, k12{}, k22{};
    double kp{}, kd{};
    double c1{}, c2{};
    double d1{};
    double dx{}, dy{}, dz{};
    // filled @ `CorrectPCAs_Z` //
    double px1{}, py1{};
    double px2{}, py2{};
    double p12{}, p22{};
    double lp1p2{};
    double detp{};
    double ldrp1{}, ldrp2{};
    double a1{}, a2{};
    // filled @ `ComputeDerivatives_XY` //
    double aa1{}, bb1{}, cc1{}, dd1{};
    double aa2{}, bb2{}, cc2{}, dd2{};
    // filled @ `FastPCAs_XY`, but here for memory alignment reasons //
    int w_sign{};
    // status flag //
    int pca_dz_worked{0};
};

// Main Methods //

std::pair<Seed, Seed> FastPCAs_XY(const POD::Track& q1, const POD::Track& q2, double bz, Cache* cache = nullptr);
std::pair<Seed, Seed> CorrectPCAs_Z(const Seed& s1_xy, const Seed& s2_xy, Cache& c);

std::pair<Deriv, Deriv> ComputeDerivatives_XY(Cache& c);
std::pair<Deriv, Deriv> UpdateDerivatives_Z(const Seed& s1_xy, const Seed& s2_xy, const Deriv& d1_xy, const Deriv& d2_xy, const Cache& c);

// Inline Methods //

inline std::tuple<Seed, Seed, Cache> FastCorrectPCAs(const POD::Track& q1, const POD::Track& q2, double bz) {
    Cache cache;
    auto [seed1_xy, seed2_xy] = FastPCAs_XY(q1, q2, bz, &cache);
    auto [seed1, seed2] = CorrectPCAs_Z(seed1_xy, seed2_xy, cache);
    return {seed1, seed2, cache};
}

inline std::tuple<Deriv, Deriv> ComputeDerivatives(const Seed& seed1_xy, const Seed& seed2_xy, Cache& cache) {
    auto [deriv1_xy, deriv2_xy] = ComputeDerivatives_XY(cache);
    auto [deriv1, deriv2] = UpdateDerivatives_Z(seed1_xy, seed2_xy, deriv1_xy, deriv2_xy, cache);
    return {deriv1, deriv2};
}

}  // namespace T2DS::Seeder::HelixHelix
