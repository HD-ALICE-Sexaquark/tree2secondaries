#include <array>
#include <cmath>
#include <utility>

#include "common/Constants.hpp"

#if T2DS_DEBUG
#include "Utils/Logger.hxx"
#endif
#include "Seeder/SeederLineVertex.hxx"
#include "Seeder/SeederTransport.hxx"
#include "Seeder/SeederTypes.hxx"

#include "Seeder/SeederLineLine.hxx"

namespace T2DS::Seeder::LineLine {

// Find points of closest approach (PCAs) between two V0s.
// Assume they transport as straight lines.
// Arguments:
// - `s01` -- [input] neutral particle 1
// - `s02` -- [input] neutral particle 2
// - `c`   -- [output] cache, consumed by `ComputeDerivatives`
// Return: (packed as a pair of `Seed` structs)
// - `pca.xyz`, `pca.mom`              -- position and momentum at their PCAs
// - `ds`                              -- transport parameters needed to reach their PCAs
// - `theta`, `sin`, `cos`, `sB`, `cB` -- cached quantities
// NOTE: when V0s are parallel (very unlikely), return the PCA to the origin {0, 0, 0} (for now)
SeedsPair FastPCAs(const State& s01, const State& s02, Cache& c) {

    // cache //

    c.px01 = s01.px;
    c.py01 = s01.py;
    c.pz01 = s01.pz;

    c.px02 = s02.px;
    c.py02 = s02.py;
    c.pz02 = s02.pz;

    c.dx = s02.x - s01.x;
    c.dy = s02.y - s01.y;
    c.dz = s02.z - s01.z;

    c.p12 = s01.px * s01.px + s01.py * s01.py + s01.pz * s01.pz;
    c.p22 = s02.px * s02.px + s02.py * s02.py + s02.pz * s02.pz;
    c.p1p2 = s01.px * s02.px + s01.py * s02.py + s01.pz * s02.pz;

    c.drp1 = s01.px * c.dx + s01.py * c.dy + s01.pz * c.dz;
    c.drp2 = s02.px * c.dx + s02.py * c.dy + s02.pz * c.dz;

    c.detp = c.p1p2 * c.p1p2 - c.p12 * c.p22;

    // protection : fully parallel v0s //

    if (std::abs(c.detp) < Common::AbsAlmostZero) [[unlikely]] {
        LineVertex::Cache to_origin_1;
        LineVertex::Cache to_origin_2;
        return {LineVertex::FastPCA(s01, {0., 0., 0.}, to_origin_1),  //
                LineVertex::FastPCA(s02, {0., 0., 0.}, to_origin_2)};
    }

    // seeds //

    Seed seed1 = TransportLine(s01, (c.drp2 * c.p1p2 - c.drp1 * c.p22) / c.detp);
    Seed seed2 = TransportLine(s02, (c.drp2 * c.p12 - c.drp1 * c.p1p2) / c.detp);

#if T2DS_DEBUG
    Logger::Debug(__FUNCTION__, "seed1.ds = {:13.6e}", seed1.ds);
    Logger::Debug(__FUNCTION__, "seed1.(x,y,z) = {}", seed1.pca.xyz);
    Logger::Debug(__FUNCTION__, "seed1.(px,py,pz) = {}", seed1.pca.mom);
    Logger::Debug(__FUNCTION__, "seed2.ds = {:13.6e}", seed2.ds);
    Logger::Debug(__FUNCTION__, "seed2.(x,y,z) = {}", seed2.pca.xyz);
    Logger::Debug(__FUNCTION__, "seed2.(px,py,pz) = {}", seed2.pca.mom);
#endif

    return {seed1, seed2};
}

// Compute derivatives of the previous PCAs calculation.
// Argument:
// - `c` -- [input] cache filled in `FastPCAs`
// Return: (packed as a pair of `Deriv` structs)
// - `ds_dr`  -- partial derivatives of current particle's ds w.r.t. current particle's state parameters = d(ds1)/dr1
// - `ds_dr1` -- partial derivatives of current particle's ds w.r.t. other particle's state parameters = d(ds1)/dr2, d(ds2)/dr1
std::pair<Deriv, Deriv> ComputeDerivatives(const Cache& c) {

    Deriv deriv1;
    Deriv deriv2;

    if (std::abs(c.detp) < Common::AbsAlmostZero) [[unlikely]] {
        return {deriv1, deriv2};
    }

    std::array<double, 6> drp1_dr1{-c.px01, -c.py01, -c.pz01, c.dx, c.dy, c.dz};
    std::array<double, 6> drp1_dr2{c.px01, c.py01, c.pz01, 0., 0., 0.};
    std::array<double, 6> drp2_dr1{-c.px02, -c.py02, -c.pz02, 0., 0., 0.};
    std::array<double, 6> drp2_dr2{c.px02, c.py02, c.pz02, c.dx, c.dy, c.dz};
    std::array<double, 6> dp1p2_dr1{0., 0., 0., c.px02, c.py02, c.pz02};
    std::array<double, 6> dp1p2_dr2{0., 0., 0., c.px01, c.py01, c.pz01};
    std::array<double, 6> dp12_dr1{0., 0., 0., 2. * c.px01, 2. * c.py01, 2. * c.pz01};
    // std::array<double, 6> dp12_dr2{0., 0., 0., 0., 0., 0};  // NOTE: kept commented for historical reasons
    // std::array<double, 6> dp22_dr1{0., 0., 0., 0., 0., 0};  // NOTE: kept commented for historical reasons
    std::array<double, 6> dp22_dr2{0., 0., 0., 2. * c.px02, 2. * c.py02, 2. * c.pz02};
    std::array<double, 6> ddetp_dr1{0.,
                                    0.,
                                    0.,
                                    -2. * c.p22 * c.px01 + 2. * c.p1p2 * c.px02,
                                    -2. * c.p22 * c.py01 + 2. * c.p1p2 * c.py02,
                                    -2. * c.p22 * c.pz01 + 2. * c.p1p2 * c.pz02};
    std::array<double, 6> ddetp_dr2{0.,
                                    0.,
                                    0.,
                                    2. * c.p1p2 * c.px01 - 2. * c.p12 * c.px02,
                                    2. * c.p1p2 * c.py01 - 2. * c.p12 * c.py02,
                                    2. * c.p1p2 * c.pz01 - 2. * c.p12 * c.pz02};

    double a1 = c.drp2 * c.p1p2 - c.drp1 * c.p22;
    double a2 = c.drp2 * c.p12 - c.drp1 * c.p1p2;
    double detp2 = c.detp * c.detp;

    for (std::size_t i = 0; i < 6; ++i) {
        double da1_dr1 = drp2_dr1[i] * c.p1p2 + c.drp2 * dp1p2_dr1[i] - drp1_dr1[i] * c.p22;  // - drp1 * dp22_dr1[i] = 0
        double da1_dr2 = drp2_dr2[i] * c.p1p2 + c.drp2 * dp1p2_dr2[i] - drp1_dr2[i] * c.p22 - c.drp1 * dp22_dr2[i];
        double da2_dr1 = drp2_dr1[i] * c.p12 + c.drp2 * dp12_dr1[i] - drp1_dr1[i] * c.p1p2 - c.drp1 * dp1p2_dr1[i];
        double da2_dr2 = drp2_dr2[i] * c.p12 - drp1_dr2[i] * c.p1p2 - c.drp1 * dp1p2_dr2[i];  // + drp2 * dp12_dr2[i] = 0

        deriv1.ds_dr[i] = da1_dr1 / c.detp - a1 * ddetp_dr1[i] / detp2;
        deriv1.ds_dr1[i] = da1_dr2 / c.detp - a1 * ddetp_dr2[i] / detp2;
        deriv2.ds_dr1[i] = da2_dr1 / c.detp - a2 * ddetp_dr1[i] / detp2;
        deriv2.ds_dr[i] = da2_dr2 / c.detp - a2 * ddetp_dr2[i] / detp2;
    }

#if T2DS_DEBUG
    Logger::Debug(__FUNCTION__, "deriv1.ds_dr = {}", deriv1.ds_dr);
    Logger::Debug(__FUNCTION__, "deriv1.ds_dr1 = {}", deriv1.ds_dr1);
    Logger::Debug(__FUNCTION__, "deriv2.ds_dr = {}", deriv2.ds_dr);
    Logger::Debug(__FUNCTION__, "deriv2.ds_dr1 = {}", deriv2.ds_dr1);
#endif

    return {deriv1, deriv2};
}

}  // namespace T2DS::Seeder::LineLine
