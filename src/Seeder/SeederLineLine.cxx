#include "Seeder/SeederLineLine.hxx"

#include <array>
#include <cmath>
#include <utility>

#include "Math/Constants.hxx"
#include "Seeder/BaseSeeder.hxx"
#include "Seeder/SeederLineVertex.hxx"
#include "View/ViewVectorV0s.hxx"
#if T2S_DEBUG
#include "App/Logger.hxx"
#endif

namespace Tree2Secondaries::Seeder::LineLine {

// First phase. Find points of closest approach (PCAs) between two neutral particles.
// Assume they transport as straight lines.
// Arguments:
// - `n1`    -- [input] neutral particle 1
// - `n2`    -- [input] neutral particle 2
// - `cache` -- [output,optional]
// Return: (packed as a pair of `Seed` structs)
// - `pca.xyz`, `pca.mom`              -- position and momentum at their PCAs
// - `ds`                              -- transport parameters needed to reach their PCAs
// - `theta`, `sin`, `cos`, `sB`, `cB` -- cached quantities
// NOTE: when tracks are parallel, return anything. For now, the PCA to the origin {0, 0, 0}.
std::pair<Seed, Seed> FastPCAs(const View::VecV0s& n1, const View::VecV0s& n2, Cache* cache) {

    double x01 = n1.X();
    double y01 = n1.Y();
    double z01 = n1.Z();

    double x02 = n2.X();
    double y02 = n2.Y();
    double z02 = n2.Z();

    // cache //

    Cache local;
    Cache& c = cache != nullptr ? *cache : local;

    c.px01 = n1.Px();
    c.py01 = n1.Py();
    c.pz01 = n1.Pz();

    c.px02 = n2.Px();
    c.py02 = n2.Py();
    c.pz02 = n2.Pz();

    c.dx = x02 - x01;
    c.dy = y02 - y01;
    c.dz = z02 - z01;

    c.p12 = c.px01 * c.px01 + c.py01 * c.py01 + c.pz01 * c.pz01;
    c.p22 = c.px02 * c.px02 + c.py02 * c.py02 + c.pz02 * c.pz02;
    c.p1p2 = c.px01 * c.px02 + c.py01 * c.py02 + c.pz01 * c.pz02;

    c.drp1 = c.px01 * c.dx + c.py01 * c.dy + c.pz01 * c.dz;
    c.drp2 = c.px02 * c.dx + c.py02 * c.dy + c.pz02 * c.dz;

    c.detp = c.p1p2 * c.p1p2 - c.p12 * c.p22;

    // protection : fully parallel v0s //

    if (std::abs(c.detp) < Const::AbsAlmostZero) [[unlikely]] {
        return {LineVertex::FastPCA(n1, {0., 0., 0.}), LineVertex::FastPCA(n2, {0., 0., 0.})};
    }

    // seed 1 //

    Seed seed1;

    seed1.ds = (c.drp2 * c.p1p2 - c.drp1 * c.p22) / c.detp;

    seed1.theta = 0.;
    seed1.sin = 0.;
    seed1.cos = 1.;
    seed1.sB = seed1.ds;
    seed1.cB = 0.;

    seed1.pca.xyz = {x01 + c.px01 * seed1.ds, y01 + c.py01 * seed1.ds, z01 + c.pz01 * seed1.ds};
    seed1.pca.mom = {c.px01, c.py01, c.pz01};

    // seed 2 //

    Seed seed2;

    seed2.ds = (c.drp2 * c.p12 - c.drp1 * c.p1p2) / c.detp;

    seed2.theta = 0.;
    seed2.sin = 0.;
    seed2.cos = 1.;
    seed2.sB = seed2.ds;
    seed2.cB = 0.;

    seed2.pca.xyz = {x02 + c.px02 * seed2.ds, y02 + c.py02 * seed2.ds, z02 + c.pz02 * seed2.ds};
    seed2.pca.mom = {c.px02, c.py02, c.pz02};

#if T2S_DEBUG
    Logger::Debug(__FUNCTION__, "seed1.ds = {:13.6e}", seed1.ds);
    Logger::Debug(__FUNCTION__, "seed1.(x,y,z) = {}", seed1.pca.xyz);
    Logger::Debug(__FUNCTION__, "seed1.(px,py,pz) = {}", seed1.pca.mom);
    Logger::Debug(__FUNCTION__, "seed2.ds = {:13.6e}", seed2.ds);
    Logger::Debug(__FUNCTION__, "seed2.(x,y,z) = {}", seed2.pca.xyz);
    Logger::Debug(__FUNCTION__, "seed2.(px,py,pz) = {}", seed2.pca.mom);
#endif

    return {seed1, seed2};
}

// Second phase. Compute derivatives of previous PCAs calculation.
// Argument:
// - `c`      -- [input] cache filled in `FastPCAs`
// Return: (packed as a pair of `Deriv` structs)
// - `ds_dr`  -- partial derivatives of current particle's ds w.r.t. current particle's state parameters = d(ds1)/dr1
// - `ds_dr1` -- partial derivatives of current particle's ds w.r.t. other particle's state parameters = d(ds1)/dr2, d(ds2)/dr1
std::pair<Deriv, Deriv> ComputeDerivatives(const Cache& c) {

    Deriv deriv1;
    Deriv deriv2;

    if (std::abs(c.detp) < Const::AbsAlmostZero) [[unlikely]] {
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

    for (size_t i = 0; i < 6; ++i) {
        double da1_dr1 = drp2_dr1[i] * c.p1p2 + c.drp2 * dp1p2_dr1[i] - drp1_dr1[i] * c.p22;  // - drp1 * dp22_dr1[i] = 0
        double da1_dr2 = drp2_dr2[i] * c.p1p2 + c.drp2 * dp1p2_dr2[i] - drp1_dr2[i] * c.p22 - c.drp1 * dp22_dr2[i];
        double da2_dr1 = drp2_dr1[i] * c.p12 + c.drp2 * dp12_dr1[i] - drp1_dr1[i] * c.p1p2 - c.drp1 * dp1p2_dr1[i];
        double da2_dr2 = drp2_dr2[i] * c.p12 - drp1_dr2[i] * c.p1p2 - c.drp1 * dp1p2_dr2[i];  // + drp2 * dp12_dr2[i] = 0

        deriv1.ds_dr[i] = da1_dr1 / c.detp - a1 * ddetp_dr1[i] / detp2;
        deriv1.ds_dr1[i] = da1_dr2 / c.detp - a1 * ddetp_dr2[i] / detp2;
        deriv2.ds_dr1[i] = da2_dr1 / c.detp - a2 * ddetp_dr1[i] / detp2;
        deriv2.ds_dr[i] = da2_dr2 / c.detp - a2 * ddetp_dr2[i] / detp2;
    }

#if T2S_DEBUG
    Logger::Debug(__FUNCTION__, "deriv1.ds_dr = {}", deriv1.ds_dr);
    Logger::Debug(__FUNCTION__, "deriv1.ds_dr1 = {}", deriv1.ds_dr1);
    Logger::Debug(__FUNCTION__, "deriv2.ds_dr = {}", deriv2.ds_dr);
    Logger::Debug(__FUNCTION__, "deriv2.ds_dr1 = {}", deriv2.ds_dr1);
#endif

    return {deriv1, deriv2};
}

}  // namespace Tree2Secondaries::Seeder::LineLine
