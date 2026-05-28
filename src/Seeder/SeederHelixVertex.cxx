#include "Seeder/SeederHelixVertex.hxx"

#include <array>
#include <cmath>

#include "common/Constants.hpp"
#include "common/Math.hpp"
#include "common/VC_TrackView.hpp"

#include "Seeder/BaseSeeder.hxx"
#if R2DS_DEBUG
#include "App/Logger.hxx"
#endif

namespace R2DS::Seeder::HelixVertex {

// First phase. Find point of closest approach (PCA) of this particle w.r.t. an arbitrary vertex in the XY plane.
// Arguments:
// - `q`     -- [input] charged particle
// - `v`     -- [input] arbitrary vertex
// - `bz`    -- [input] z-component of homogeneous magnetic field
// - `cache` -- [output,optional]
// Return: (packed in a single `Seed` struct)
// - `pca.xyz`, `pca.mom`              -- position and momentum at their PCAs
// - `ds`                              -- transport parameters needed to reach their PCAs
// - `theta`, `sin`, `cos`, `sB`, `cB` -- cache related-quantities
Seed FastPCA_XY(const Vector::TrackView& q, const std::array<double, 3>& v, double bz, Cache* cache) {

    // cache //

    Cache local;
    Cache& c = cache != nullptr ? *cache : local;

    c.x0 = q.X();
    c.y0 = q.Y();
    c.z0 = q.Z();

    c.px0 = q.Px();
    c.py0 = q.Py();
    c.pz0 = q.Pz();

    c.pt2 = c.px0 * c.px0 + c.py0 * c.py0;

    c.dx = v[0] - c.x0;
    c.dy = v[1] - c.y0;
    c.dz = v[2] - c.z0;

    c.a = c.dx * c.px0 + c.dy * c.py0;

    c.bq = bz * static_cast<double>(q.Charge()) * Common::Kappa;
    c.abq = c.a * c.bq;

    // prepare seed //

    Seed seed;

    seed.theta = std::atan2(c.abq, c.pt2 + c.bq * (c.dy * c.px0 - c.dx * c.py0));
    std::tie(seed.sin, seed.cos) = Common::Math::sincos(seed.theta);
    seed.sB = seed.sin / c.bq;
    seed.cB = (1. - seed.cos) / c.bq;

    seed.ds = seed.theta / c.bq;

    seed.pca.xyz = {c.x0 + seed.sB * c.px0 + seed.cB * c.py0, c.y0 - seed.cB * c.px0 + seed.sB * c.py0, c.z0 + seed.ds * c.pz0};
    seed.pca.mom = {seed.cos * c.px0 + seed.sin * c.py0, -seed.sin * c.px0 + seed.cos * c.py0, c.pz0};

#if R2DS_DEBUG
    Logger::Debug(__FUNCTION__, "seed.ds = {:13.6e}", seed.ds);
    Logger::Debug(__FUNCTION__, "seed.(x,y,z) = {}", seed.pca.xyz);
    Logger::Debug(__FUNCTION__, "seed.(px,py,pz) = {}", seed.pca.mom);
#endif

    return seed;
}

// Second phase. Correct previous seed by adding the z-component as a correction.
// Arguments:
// - `s_xy` -- [input] seed calculated in `FastPCA_XY`
// - `c`    -- [input/output] cache filled in `FastPCA_XY`
Seed CorrectPCA_Z(const Seed& s_xy, Cache& c) {

    // cache (1) //

    c.bbq = c.bq * (c.dx * c.py0 - c.dy * c.px0) - c.pt2;
    c.cbq = c.bbq * s_xy.cos - c.abq * s_xy.sin - c.pz0 * c.pz0;

    // protection //

    if (std::abs(c.cbq) < Common::AbsAlmostZero) return s_xy;

    // cache (2) //

    c.sz = (s_xy.ds * c.pz0 - c.dz) * c.pz0 / c.cbq;

    // update seed //

    Seed out;

    out.ds = s_xy.ds + c.sz;

    out.theta = c.bq * out.ds;
    std::tie(out.sin, out.cos) = Common::Math::sincos(out.theta);
    out.sB = out.sin / c.bq;
    out.cB = (1. - out.cos) / c.bq;

    out.pca.xyz = {c.x0 + out.sB * c.px0 + out.cB * c.py0, c.y0 - out.cB * c.px0 + out.sB * c.py0, c.z0 + out.ds * c.pz0};
    out.pca.mom = {out.cos * c.px0 + out.sin * c.py0, -out.sin * c.px0 + out.cos * c.py0, c.pz0};

    // update status flag //

    c.pca_dz_worked = 1;

#if R2DS_DEBUG
    Logger::Debug(__FUNCTION__, "seed.ds = {:13.6e}", out.ds);
    Logger::Debug(__FUNCTION__, "seed.(x,y,z) = {}", out.pca.xyz);
    Logger::Debug(__FUNCTION__, "seed.(px,py,pz) = {}", out.pca.mom);
#endif

    return out;
}

// Third phase. Compute derivatives from `FastPCA_XY`.
// Argument:
// - `c` -- [input/output] cache filled in `FastPCA_XY` and `CorrectPCA_Z`
// Return: (packed in a single `Deriv` struct)
// - `ds_dr`  -- partial derivatives of current particle's ds w.r.t. current particle's state parameters = d(ds1)/dr1, d(ds2)/dr2
// - `ds_dr1` -- partial derivatives of current particle's ds w.r.t. other particle's state parameters = d(ds2)/dr1, d(ds1)/dr2
Deriv ComputeDerivatives_XY(Cache& c) {

    Deriv out;

    // cache //

    c.den = c.abq * c.abq + c.bbq * c.bbq;

    // protection //

    if (c.den < Common::AbsAlmostZero) return out;

    // compute deriv //

    out.ds_dr = {(c.px0 * c.bbq - c.py0 * c.abq) / c.den,
                 (c.px0 * c.abq + c.py0 * c.bbq) / c.den,
                 0.,
                 -(c.dx * c.bbq + c.dy * c.abq + 2. * c.px0 * c.a) / c.den,
                 (c.dx * c.abq - c.dy * c.bbq - 2. * c.py0 * c.a) / c.den,
                 0.};

    out.ds_dr1 = {-out.ds_dr[0], -out.ds_dr[1], -out.ds_dr[2], 0., 0., 0.};

#if R2DS_DEBUG
    Logger::Debug(__FUNCTION__, "(end) deriv.ds_dr = {}", out.ds_dr);
    Logger::Debug(__FUNCTION__, "(end) deriv.ds_dr1 = {}", out.ds_dr1);
#endif

    return out;
}

// Final phase. Update derivatives considering correction made in `CorrectPCA_Z`.
// Arguments:
// - `s_xy` -- [input] seed calculated in `FastPCA_XY`
// - `d_xy` -- [input] derivatives calculated in `ComputeDerivatives_XY`
// - `c`    -- [input] cache filled in `FastPCA_XY`, `CorrectPCA_Z` and `ComputeDerivatives_XY`
Deriv UpdateDerivatives_Z(const Seed& s_xy, const Deriv& d_xy, const Cache& c) {

    // if `CorrectPCAs_Z` didn't work, skip this method //

    if (c.pca_dz_worked == 0) return d_xy;

    std::array<double, 6> dc_dr = {
        -c.bq * c.py0 * s_xy.cos - c.bbq * s_xy.sin * c.bq * d_xy.ds_dr[0] + c.px0 * c.bq * s_xy.sin - c.abq * s_xy.cos * c.bq * d_xy.ds_dr[0],
        c.bq * c.px0 * s_xy.cos - c.bbq * s_xy.sin * c.bq * d_xy.ds_dr[1] + c.py0 * c.bq * s_xy.sin - c.abq * s_xy.cos * c.bq * d_xy.ds_dr[1],
        0.,
        (-c.bq * c.dy - 2. * c.px0) * s_xy.cos - c.bbq * s_xy.sin * c.bq * d_xy.ds_dr[3] - c.dx * c.bq * s_xy.sin -
            c.abq * s_xy.cos * c.bq * d_xy.ds_dr[3],
        (c.bq * c.dx - 2. * c.py0) * s_xy.cos - c.bbq * s_xy.sin * c.bq * d_xy.ds_dr[4] - c.dy * c.bq * s_xy.sin -
            c.abq * s_xy.cos * c.bq * d_xy.ds_dr[4],
        -2. * c.pz0};

    Deriv out = d_xy;

    for (size_t i = 0; i < 6; ++i) {
        out.ds_dr[i] += c.pz0 * c.pz0 * d_xy.ds_dr[i] / c.cbq - c.sz * dc_dr[i] / c.cbq;
    }
    out.ds_dr[2] += c.pz0 / c.cbq;
    out.ds_dr[5] += (2. * c.pz0 * s_xy.ds - c.dz) / c.cbq;

    out.ds_dr1 = {-out.ds_dr[0], -out.ds_dr[1], -out.ds_dr[2], 0., 0., 0.};

#if R2DS_DEBUG
    Logger::Debug(__FUNCTION__, "(end) deriv.ds_dr = {}", out.ds_dr);
    Logger::Debug(__FUNCTION__, "(end) deriv.ds_dr1 = {}", out.ds_dr1);
#endif

    return out;
}

}  // namespace R2DS::Seeder::HelixVertex
