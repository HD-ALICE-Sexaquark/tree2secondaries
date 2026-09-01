#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <utility>

#include "common/Constants.hpp"

#if T2DS_DEBUG
#include "Utils/Logger.hxx"
#endif
#include "Seeder/SeederTransport.hxx"
#include "Seeder/SeederTypes.hxx"

#include "Seeder/SeederHelixLine.hxx"

namespace T2DS::Seeder::HelixLine {

// First phase. Find the points of closest approach (PCA) between two particles in the XY plane.
// One particle is charged under a constant magnetic field and the other one is neutral.
// Transport the former as a helix and the latter as a line.
// Arguments:
// - `s01` -- [input] charged particle, transports as helix
// - `s02` -- [input] neutral particle, transports as straight line
// - `bz`  -- [input] z-component of homogeneous magnetic field
// - `c`   -- [output] cache
// Return: (packed as a pair of `Seed` structs)
// - `ds`  -- transport parameters
// - `pca` -- points of closest approach (position and momentum)
// - `theta`, `sin`, `cos`, `sB`, `cB` -- cached ds computation variables
SeedsPair FastPCAs_XY(const State& s01, const State& s02, double bz, Cache& c) {

    // cache //

    c.bq1 = bz * static_cast<double>(s01.charge) * Common::Kappa;

    c.x01 = s01.x;
    c.y01 = s01.y;
    c.z01 = s01.z;
    c.px01 = s01.px;
    c.py01 = s01.py;
    c.pz01 = s01.pz;
    c.pt12 = c.px01 * c.px01 + c.py01 * c.py01;

    c.x02 = s02.x;
    c.y02 = s02.y;
    c.z02 = s02.z;
    c.px02 = s02.px;
    c.py02 = s02.py;
    c.pz02 = s02.pz;
    c.pt22 = c.px02 * c.px02 + c.py02 * c.py02;

    c.dx0 = c.x01 - c.x02;
    c.dy0 = c.y01 - c.y02;
    c.drp2 = c.dx0 * c.px02 + c.dy0 * c.py02;
    c.dxyp2 = c.dx0 * c.py02 - c.dy0 * c.px02;
    c.p1p2 = c.px01 * c.px02 + c.py01 * c.py02;
    c.dp1p2 = c.px01 * c.py02 - c.px02 * c.py01;

    c.k11 = -c.dp1p2;
    c.k21 = -c.bq1 * c.p1p2;
    c.k12 = c.bq1 * c.drp2 - c.dp1p2;
    c.k22 = -c.bq1 * c.pt22;

    c.kp = -c.dxyp2 * c.bq1 - c.p1p2;
    c.kd = c.kp;
    c.c1 = -c.bq1 * c.kd;
    c.c2 = c.pt22 * c.bq1;

    c.d1 = std::sqrt(std::max(c.pt12 * c.pt22 - c.kd * c.kd, 0.));

    // prepare seeds //
    // -- by testing minimum via calculating 3D distance

    Seed seed1_xy;
    Seed seed2_xy;

    double dca_sq = std::numeric_limits<double>::max();

    for (auto sign : {+1, -1}) {

        // charged particle //

        Seed tmp1 = TransportHelixFromTheta(
            s01, std::atan2(c.bq1 * (c.k11 * c.c1 + sign * c.k21 * c.d1), sign * c.bq1 * c.k11 * c.d1 * c.bq1 - c.k21 * c.c1), c.bq1);

        // neutral particle //

        double denom = -c.k22 * c.c2;
        if (std::abs(denom) < Common::AbsAlmostZero) continue;  // skip this sign

        Seed tmp2 = TransportLine(s02, (c.k12 * c.c2 + sign * c.k22 * c.d1) / denom);

        // distance of closest approach (DCA) //
        double tmp_dx = tmp2.pca.xyz[0] - tmp1.pca.xyz[0];
        double tmp_dy = tmp2.pca.xyz[1] - tmp1.pca.xyz[1];
        double tmp_dz = tmp2.pca.xyz[2] - tmp1.pca.xyz[2];
        double tmp_dca_sq = tmp_dx * tmp_dx + tmp_dy * tmp_dy + tmp_dz * tmp_dz;

        if (tmp_dca_sq < dca_sq) {
            seed1_xy = tmp1;
            seed2_xy = tmp2;

            c.dx = tmp_dx;
            c.dy = tmp_dy;
            c.dz = tmp_dz;

            c.w_sign = sign;
            dca_sq = tmp_dca_sq;
        }
    }

#if T2DS_DEBUG
    Logger::Debug(__FUNCTION__, "seed1_xy.ds = {:13.6e}", seed1_xy.ds);
    Logger::Debug(__FUNCTION__, "seed1_xy.(x,y,z) = {}", seed1_xy.pca.xyz);
    Logger::Debug(__FUNCTION__, "seed1_xy.(px,py,pz) = {}", seed1_xy.pca.mom);
    Logger::Debug(__FUNCTION__, "seed2_xy.ds = {:13.6e}", seed2_xy.ds);
    Logger::Debug(__FUNCTION__, "seed2_xy.(x,y,z) = {}", seed2_xy.pca.xyz);
    Logger::Debug(__FUNCTION__, "seed2_xy.(px,py,pz) = {}", seed2_xy.pca.mom);
#endif

    return {seed1_xy, seed2_xy};
}

// Second phase. Add z-component as small correction.
// Arguments:
// - `seeds_xy` -- [input] seeds of particle 1 and 2 calculated in `FastPCAs_XY`
// - `c`        -- [input/output] cache filled in `FastPCAs_XY`
// Return: (packed as a pair of `Seed` structs)
// - `ds`  -- transport parameters
// - `pca` -- points of closest approach (position and momentum)
// - `theta`, `sin`, `cos`, `sB`, `cB` -- cached ds computation variables
SeedsPair CorrectPCAs_Z(const SeedsPair& seeds_xy, Cache& c) {

    // cache (1) //

    c.px1 = seeds_xy.first.pca.mom[0];
    c.py1 = seeds_xy.first.pca.mom[1];

    c.p12 = c.px1 * c.px1 + c.py1 * c.py1 + c.pz01 * c.pz01;
    c.p22 = c.px02 * c.px02 + c.py02 * c.py02 + c.pz02 * c.pz02;
    c.lp1p2 = c.px1 * c.px02 + c.py1 * c.py02 + c.pz01 * c.pz02;

    c.detp = c.lp1p2 * c.lp1p2 - c.p12 * c.p22;

    if (std::abs(c.detp) < Common::AbsAlmostZero) return seeds_xy;  // protection

    // cache (2) //

    c.ldrp2 = c.px02 * c.dx + c.py02 * c.dy + c.pz02 * c.dz;
    c.ldrp1 = c.px1 * c.dx + c.py1 * c.dy + c.pz01 * c.dz;
    c.a2 = c.ldrp2 * c.p12 - c.ldrp1 * c.lp1p2;
    c.a1 = c.ldrp2 * c.lp1p2 - c.ldrp1 * c.p22;

    // -- charged particle
    Seed s1 = TransportHelixFromDs(c.x01, c.y01, c.z01, c.px01, c.py01, c.pz01, seeds_xy.first.ds + c.a1 / c.detp, c.bq1);

    // -- neutral particle
    Seed s2 = TransportLine(c.x02, c.y02, c.z02, c.px02, c.py02, c.pz02, seeds_xy.second.ds + c.a2 / c.detp);

    // update status flag //

    c.pca_dz_worked = true;

#if T2DS_DEBUG
    Logger::Debug(__FUNCTION__, "seed1.ds = {:13.6e}", s1.ds);
    Logger::Debug(__FUNCTION__, "seed1.(x,y,z) = {}", s1.pca.xyz);
    Logger::Debug(__FUNCTION__, "seed1.(px,py,pz) = {}", s1.pca.mom);
    Logger::Debug(__FUNCTION__, "seed2.ds = {:13.6e}", s2.ds);
    Logger::Debug(__FUNCTION__, "seed2.(x,y,z) = {}", s2.pca.xyz);
    Logger::Debug(__FUNCTION__, "seed2.(px,py,pz) = {}", s2.pca.mom);
#endif

    return {s1, s2};
}

// Third phase. Compute derivatives of `FastPCAs_XY`.
// Argument:
// - `c` -- [input/output] cache filled in `FastPCAs_XY`
// Return: (packed as a pair of `Deriv` structs)
// - `ds_dr`  -- partial derivatives of current particle's ds w.r.t. current particle's state parameters = d(ds1)/dr1
// - `ds_dr1` -- partial derivatives of current particle's ds w.r.t. other particle's state parameters = d(ds1)/dr2, d(ds2)/dr1
std::pair<Deriv, Deriv> ComputeDerivatives_XY(Cache& c) {

    Deriv deriv1;
    Deriv deriv2;

    std::array<double, 6> dk11dr1{0., 0., 0., -c.py02, c.px02, 0.};
    std::array<double, 6> dk11dr2{0., 0., 0., c.py01, -c.px01, 0.};
    std::array<double, 6> dk12dr1{c.bq1 * c.px02, c.bq1 * c.py02, 0., -c.py02, c.px02, 0.};
    std::array<double, 6> dk12dr2{-c.bq1 * c.px02, -c.bq1 * c.py02, 0., c.bq1 * c.dx0 + c.py01, c.bq1 * c.dy0 - c.px01, 0.};
    std::array<double, 6> dk21dr1{0., 0., 0., -c.bq1 * c.px02, -c.bq1 * c.py02, 0.};
    std::array<double, 6> dk21dr2{0., 0., 0., -c.bq1 * c.px01, -c.bq1 * c.py01, 0.};
    // std::array<double, 6> dk22dr1{0., 0., 0., 0., 0., 0.};  // NOTE: kept commented for historical reasons
    std::array<double, 6> dk22dr2{0., 0., 0., -2. * c.bq1 * c.px02, -2. * c.bq1 * c.py02, 0.};

    std::array<double, 6> dkddr1{-c.bq1 * c.py02, c.bq1 * c.px02, 0., -c.px02, -c.py02, 0.};
    std::array<double, 6> dkddr2{c.bq1 * c.py02, -c.bq1 * c.px02, 0., c.bq1 * c.dy0 - c.px01, -c.bq1 * c.dx0 - c.py01, 0.};

    std::array<double, 6> dc1dr1{c.bq1 * c.bq1 * c.py02, -c.bq1 * c.bq1 * c.px02, 0., c.bq1 * c.px02, c.bq1 * c.py02, 0.};
    std::array<double, 6> dc1dr2{-c.bq1 * c.bq1 * c.py02,           c.bq1 * c.bq1 * c.px02,           0.,
                                 -c.bq1 * (c.bq1 * c.dy0 - c.px01), c.bq1 * (c.bq1 * c.dx0 + c.py01), 0.};

    // std::array<double, 6> dc2dr1{0., 0., 0., 0., 0., 0.};  // NOTE: kept commented for historical reasons
    std::array<double, 6> dc2dr2{0., 0., 0., 2. * c.bq1 * c.px02, 2. * c.bq1 * c.py02, 0.};

    std::array<double, 6> dd1dr1{};
    std::array<double, 6> dd1dr2{};
    if (c.d1 > 0.) {
        double kd_sq_d1{-c.kd / c.d1};
        for (std::size_t i = 0; i < 6; ++i) {
            dd1dr1[i] = kd_sq_d1 * dkddr1[i];
            dd1dr2[i] = kd_sq_d1 * dkddr2[i];
        }
        dd1dr1[3] += c.px01 * c.pt22 / c.d1;
        dd1dr1[4] += c.py01 * c.pt22 / c.d1;
        dd1dr2[3] += c.px02 * c.pt12 / c.d1;
        dd1dr2[4] += c.py02 * c.pt12 / c.d1;
    }

    // charged particle //

    c.aa1 = c.bq1 * (c.k11 * c.c1 + c.w_sign * c.k21 * c.d1);
    c.bb1 = c.w_sign * c.bq1 * c.k11 * c.d1 * c.bq1 - c.k21 * c.c1;
    c.cc1 = c.aa1 * c.aa1 + c.bb1 * c.bb1;
    c.dd1 = c.cc1 > 0. ? (1. / c.bq1 * 1. / c.cc1) : 0.;

    for (std::size_t i = 0; i < 6; ++i) {
        double daa1_dr1 = c.bq1 * (dk11dr1[i] * c.c1 + c.k11 * dc1dr1[i] + c.w_sign * dk21dr1[i] * c.d1 + c.w_sign * c.k21 * dd1dr1[i]);
        double daa1_dr2 = c.bq1 * (dk11dr2[i] * c.c1 + c.k11 * dc1dr2[i] + c.w_sign * dk21dr2[i] * c.d1 + c.w_sign * c.k21 * dd1dr2[i]);
        double dbb1_dr1 = c.w_sign * c.bq1 * c.bq1 * (dk11dr1[i] * c.d1 + c.k11 * dd1dr1[i]) - (dk21dr1[i] * c.c1 + c.k21 * dc1dr1[i]);
        double dbb1_dr2 = c.w_sign * c.bq1 * c.bq1 * (dk11dr2[i] * c.d1 + c.k11 * dd1dr2[i]) - (dk21dr2[i] * c.c1 + c.k21 * dc1dr2[i]);

        deriv1.ds_dr[i] = c.dd1 * (daa1_dr1 * c.bb1 - dbb1_dr1 * c.aa1);
        deriv1.ds_dr1[i] = c.dd1 * (daa1_dr2 * c.bb1 - dbb1_dr2 * c.aa1);
    }

    // neutral particle //

    c.na = c.k12 * c.c2 + c.w_sign * c.k22 * c.d1;
    c.nb = -c.k22 * c.c2;
    c.nb2 = c.nb * c.nb;

    if (std::abs(c.nb) < Common::AbsAlmostZero) return {deriv1, deriv2};  // protection

    for (std::size_t i = 0; i < 6; ++i) {
        double dna_dr1 = dk12dr1[i] * c.c2 + c.w_sign * c.k22 * dd1dr1[i];  // + c.k12 * dc2dr1[i] + c.w_sign * dk22dr1[i] * c.d1 = 0
        double dna_dr2 = dk12dr2[i] * c.c2 + c.k12 * dc2dr2[i] + c.w_sign * (c.k22 * dd1dr2[i] + dk22dr2[i] * c.d1);
        // double dnb2_dr1 = -dk22dr1[i] * c.c2 - c.k22 * dc2dr1[i]; = 0
        double dnb2_dr2 = -dk22dr2[i] * c.c2 - c.k22 * dc2dr2[i];

        deriv2.ds_dr1[i] = dna_dr1 / c.nb;  // - dnb2_dr1 * c.na / c.nb2 = 0
        deriv2.ds_dr[i] = dna_dr2 / c.nb - dnb2_dr2 * c.na / c.nb2;
    }

#if T2DS_DEBUG
    Logger::Debug(__FUNCTION__, "deriv1.ds_dr = {}", deriv1.ds_dr);
    Logger::Debug(__FUNCTION__, "deriv1.ds_dr1 = {}", deriv1.ds_dr1);
    Logger::Debug(__FUNCTION__, "deriv2.ds_dr = {}", deriv2.ds_dr);
    Logger::Debug(__FUNCTION__, "deriv2.ds_dr1 = {}", deriv2.ds_dr1);
#endif

    return {deriv1, deriv2};
}

// Final phase. Update derivatives considering the correction made in `CorrectPCAs_Z`.
// Arguments:
// - `seeds_xy` -- [input] seeds of particle 1 and 2 calculated in `FastPCAs_XY`
// - `d1_xy`    -- [input] derivatives of particle 1 calculated in `ComputeDerivatives_XY`
// - `d2_xy`    -- [input] derivatives of particle 2 calculated in `ComputeDerivatives_XY`
// - `c`        -- [input] cache filled in `FastPCAs_XY`, `CorrectPCAs_Z` and `ComputeDerivatives_XY`
// Return: (packed as a pair of `Deriv` structs)
// - `ds_dr`  -- partial derivatives of current particle's ds w.r.t. current particle's state parameters = d(ds1)/dr1
// - `ds_dr1` -- partial derivatives of current particle's ds w.r.t. other particle's state parameters = d(ds1)/dr2, d(ds2)/dr1
std::pair<Deriv, Deriv> UpdateDerivatives_Z(const SeedsPair& seeds_xy, const Deriv& d1_xy, const Deriv& d2_xy, const Cache& c) {

    // if `CorrectPCAs_Z` didn't work, skip this method //

    if (!c.pca_dz_worked) return {d1_xy, d2_xy};

    // update derivatives //

    double lp1p2_ds0 = c.bq1 * (c.px02 * c.py1 - c.py02 * c.px1);
    // double lp1p2_ds1 = 0.; // commented out for historical reasons
    double ldrp1_ds0 = -c.p12 + c.bq1 * (c.py1 * c.dx - c.px1 * c.dy);
    double ldrp1_ds1 = c.lp1p2;
    double ldrp2_ds0 = -c.lp1p2;
    double ldrp2_ds1 = c.p22;
    double detp_ds0 = 2. * c.lp1p2 * lp1p2_ds0;
    // double detp_ds1 = 2. * c.lp1p2 * lp1p2_ds1; // commented out for historical reasons
    double a1_ds0 = ldrp2_ds0 * c.lp1p2 + c.ldrp2 * lp1p2_ds0 - ldrp1_ds0 * c.p22;
    double a1_ds1 = ldrp2_ds1 * c.lp1p2 - ldrp1_ds1 * c.p22;  // + c.ldrp2 * lp1p2_ds1 = 0
    double a2_ds0 = ldrp2_ds0 * c.p12 - ldrp1_ds0 * c.lp1p2 - c.ldrp1 * lp1p2_ds0;
    double a2_ds1 = ldrp2_ds1 * c.p12 - ldrp1_ds1 * c.lp1p2;  // - c.ldrp1 * lp1p2_ds1 = 0

    double detp2 = c.detp * c.detp;
    double dsl1ds0 = a1_ds0 / c.detp - c.a1 * detp_ds0 / detp2;
    double dsl1ds1 = a1_ds1 / c.detp;  // - c.a1 * detp_ds1 / detp2 = 0
    double dsl2ds0 = a2_ds0 / c.detp - c.a2 * detp_ds0 / detp2;
    double dsl2ds1 = a2_ds1 / c.detp;  // - c.a2 * detp_ds1 / detp2 = 0

    std::array<double, 6> dsldr0{};
    std::array<double, 6> dsldr1{};
    std::array<double, 6> dsldr2{};
    std::array<double, 6> dsldr3{};

    for (std::size_t i = 0; i < 6; ++i) {
        dsldr0[i] = dsl1ds0 * d1_xy.ds_dr[i] + dsl1ds1 * d2_xy.ds_dr1[i];
        dsldr1[i] = dsl1ds0 * d1_xy.ds_dr1[i] + dsl1ds1 * d2_xy.ds_dr[i];
        dsldr2[i] = dsl2ds0 * d1_xy.ds_dr[i] + dsl2ds1 * d2_xy.ds_dr1[i];
        dsldr3[i] = dsl2ds0 * d1_xy.ds_dr1[i] + dsl2ds1 * d2_xy.ds_dr[i];
    }

    Deriv d1 = d1_xy;
    Deriv d2 = d2_xy;

    for (std::size_t i = 0; i < 6; ++i) {
        d1.ds_dr[i] += dsldr0[i];
        d1.ds_dr1[i] += dsldr1[i];
        d2.ds_dr1[i] += dsldr2[i];
        d2.ds_dr[i] += dsldr3[i];
    }

    const auto& [s1_xy, s2_xy] = seeds_xy;

    std::array<double, 6> lp1p2_dr0{0., 0., 0., s1_xy.cos * c.px02 - c.py02 * s1_xy.sin, s1_xy.cos * c.py02 + c.px02 * s1_xy.sin, c.pz02};
    std::array<double, 6> lp1p2_dr1{0., 0., 0., c.px1, c.py1, c.pz01};
    std::array<double, 6> ldrp1_dr0{-c.px1,
                                    -c.py1,
                                    -c.pz01,
                                    s1_xy.cB * c.py1 - c.px1 * s1_xy.sB + s1_xy.cos * c.dx - s1_xy.sin * c.dy,
                                    -s1_xy.cB * c.px1 - c.py1 * s1_xy.sB + s1_xy.sin * c.dx + s1_xy.cos * c.dy,
                                    -s1_xy.ds * c.pz01 + c.dz};
    std::array<double, 6> ldrp1_dr1{c.px1, c.py1, c.pz01, c.px1 * s2_xy.ds, c.py1 * s2_xy.ds, s2_xy.ds * c.pz01};
    std::array<double, 6> ldrp2_dr0{
        -c.px02, -c.py02, -c.pz02, s1_xy.cB * c.py02 - c.px02 * s1_xy.sB, -s1_xy.cB * c.px02 - c.py02 * s1_xy.sB, -s1_xy.ds * c.pz02};
    std::array<double, 6> ldrp2_dr1{c.px02, c.py02, c.pz02, c.px02 * s2_xy.ds + c.dx, c.py02 * s2_xy.ds + c.dy, c.dz + s2_xy.ds * c.pz02};
    std::array<double, 6> p12_dr0{0., 0., 0., 2. * c.px01, 2. * c.py01, 2. * c.pz01};
    std::array<double, 6> p22_dr1{0., 0., 0., 2. * c.px02, 2. * c.py02, 2. * c.pz02};

    for (std::size_t i = 0; i < 6; ++i) {
        double a1_dr0 = ldrp2_dr0[i] * c.lp1p2 + c.ldrp2 * lp1p2_dr0[i] - ldrp1_dr0[i] * c.p22;
        double a1_dr1 = ldrp2_dr1[i] * c.lp1p2 + c.ldrp2 * lp1p2_dr1[i] - ldrp1_dr1[i] * c.p22 - c.ldrp1 * p22_dr1[i];
        double a2_dr0 = ldrp2_dr0[i] * c.p12 + c.ldrp2 * p12_dr0[i] - ldrp1_dr0[i] * c.lp1p2 - c.ldrp1 * lp1p2_dr0[i];
        double a2_dr1 = ldrp2_dr1[i] * c.p12 - ldrp1_dr1[i] * c.lp1p2 - c.ldrp1 * lp1p2_dr1[i];
        double detp_dr0 = 2. * c.lp1p2 * lp1p2_dr0[i] - p12_dr0[i] * c.p22;
        double detp_dr1 = 2. * c.lp1p2 * lp1p2_dr1[i] - c.p12 * p22_dr1[i];

        d1.ds_dr[i] += a1_dr0 / c.detp - c.a1 * detp_dr0 / detp2;
        d1.ds_dr1[i] += a1_dr1 / c.detp - c.a1 * detp_dr1 / detp2;
        d2.ds_dr1[i] += a2_dr0 / c.detp - c.a2 * detp_dr0 / detp2;
        d2.ds_dr[i] += a2_dr1 / c.detp - c.a2 * detp_dr1 / detp2;
    }

#if T2DS_DEBUG
    Logger::Debug(__FUNCTION__, "deriv1.ds_dr  = {}", d1.ds_dr);
    Logger::Debug(__FUNCTION__, "deriv1.ds_dr1 = {}", d1.ds_dr1);
    Logger::Debug(__FUNCTION__, "deriv2.ds_dr  = {}", d2.ds_dr);
    Logger::Debug(__FUNCTION__, "deriv2.ds_dr1 = {}", d2.ds_dr1);
#endif

    return {d1, d2};
}

}  // namespace T2DS::Seeder::HelixLine
