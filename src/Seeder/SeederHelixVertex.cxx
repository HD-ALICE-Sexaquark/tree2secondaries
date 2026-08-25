#include <array>
#include <cmath>

#include "common/Constants.hpp"

#if T2DS_DEBUG
#include "App/Logger.hxx"
#endif
#include "Seeder/SeederTransport.hxx"
#include "Seeder/SeederTypes.hxx"

#include "Seeder/SeederHelixVertex.hxx"

namespace T2DS::Seeder::HelixVertex {

// First phase. Find point of closest approach (PCA) of this particle w.r.t. an arbitrary vertex in the XY plane.
// Arguments:
// - `s0`  -- [input] charged particle
// - `vtx` -- [input] arbitrary vertex
// - `bz`  -- [input] z-component of homogeneous magnetic field
// - `c`   -- [output] cache
// Return: (packed in a single `Seed` struct)
// - `pca.xyz`, `pca.mom`              -- position and momentum at the PCA
// - `ds`                              -- transport parameters needed to reach PCA
// - `theta`, `sin`, `cos`, `sB`, `cB` -- cache related-quantities
Seed FastPCA_XY(const State& s0, const std::array<double, 3>& vtx, double bz, Cache& c) {

    // cache //

    c.x0 = s0.x;
    c.y0 = s0.y;
    c.z0 = s0.z;

    c.px0 = s0.px;
    c.py0 = s0.py;
    c.pz0 = s0.pz;

    c.pt2 = c.px0 * c.px0 + c.py0 * c.py0;

    c.dx = vtx[0] - c.x0;
    c.dy = vtx[1] - c.y0;
    c.dz = vtx[2] - c.z0;

    c.a = c.dx * c.px0 + c.dy * c.py0;

    c.bq = bz * static_cast<double>(s0.charge) * Common::Kappa;
    c.abq = c.a * c.bq;

    // prepare seed //

    Seed out = TransportHelixFromTheta(s0, std::atan2(c.abq, c.pt2 + c.bq * (c.dy * c.px0 - c.dx * c.py0)), c.bq);

#if T2DS_DEBUG
    Logger::Debug(__FUNCTION__, "out.ds = {:13.6e}", out.ds);
    Logger::Debug(__FUNCTION__, "out.(x,y,z) = {}", out.pca.xyz);
    Logger::Debug(__FUNCTION__, "out.(px,py,pz) = {}", out.pca.mom);
#endif

    return out;
}

// Second phase. Correct previous seed by adding the z-component as a correction.
// Arguments:
// - `s_xy` -- [input] seed calculated in `FastPCA_XY`
// - `c`    -- [input/output] cache filled in `FastPCA_XY`
Seed CorrectPCA_Z(const Seed& s_xy, Cache& c) {

    // cache (1) //

    c.bbq = c.bq * (c.dx * c.py0 - c.dy * c.px0) - c.pt2;
    c.cbq = c.bbq * s_xy.cos - c.abq * s_xy.sin - c.pz0 * c.pz0;

    if (std::abs(c.cbq) < Common::AbsAlmostZero) return s_xy;  // protection

    // cache (2) //

    c.sz = (s_xy.ds * c.pz0 - c.dz) * c.pz0 / c.cbq;

    // update seed //

    Seed out = TransportHelixFromDs(c.x0, c.y0, c.z0, c.px0, c.py0, c.pz0, s_xy.ds + c.sz, c.bq);

    // update status flag //

    c.pca_dz_worked = true;

#if T2DS_DEBUG
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
// - `ds_dr`  -- partial derivatives of current particle's ds w.r.t. current particle's state parameters
// - `ds_dr1` -- partial derivatives of current particle's ds w.r.t. the vertex; i.e., negated position derivatives
Deriv ComputeDerivatives_XY(Cache& c) {

    Deriv out;

    // cache //

    c.den = c.abq * c.abq + c.bbq * c.bbq;

    if (c.den < Common::AbsAlmostZero) return out;  // protection

    // compute deriv //

    out.ds_dr = {(c.px0 * c.bbq - c.py0 * c.abq) / c.den,
                 (c.px0 * c.abq + c.py0 * c.bbq) / c.den,
                 0.,
                 -(c.dx * c.bbq + c.dy * c.abq + 2. * c.px0 * c.a) / c.den,
                 (c.dx * c.abq - c.dy * c.bbq - 2. * c.py0 * c.a) / c.den,
                 0.};

    out.ds_dr1 = {-out.ds_dr[0], -out.ds_dr[1], -out.ds_dr[2], 0., 0., 0.};

#if T2DS_DEBUG
    Logger::Debug(__FUNCTION__, "(end) deriv.ds_dr = {}", out.ds_dr);
    Logger::Debug(__FUNCTION__, "(end) deriv.ds_dr1 = {}", out.ds_dr1);
#endif

    return out;
}

// Final phase. Update derivatives considering the correction made in `CorrectPCA_Z`.
// Arguments:
// - `s_xy` -- [input] seed calculated in `FastPCA_XY`
// - `d_xy` -- [input] derivatives calculated in `ComputeDerivatives_XY`
// - `c`    -- [input] cache filled in `FastPCA_XY`, `CorrectPCA_Z` and `ComputeDerivatives_XY`
Deriv UpdateDerivatives_Z(const Seed& s_xy, const Deriv& d_xy, const Cache& c) {

    // if `CorrectPCA_Z` didn't work, skip this method //

    if (!c.pca_dz_worked) return d_xy;

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

    for (std::size_t i = 0; i < 6; ++i) {
        out.ds_dr[i] += c.pz0 * c.pz0 * d_xy.ds_dr[i] / c.cbq - c.sz * dc_dr[i] / c.cbq;
    }
    out.ds_dr[2] += c.pz0 / c.cbq;
    out.ds_dr[5] += (2. * c.pz0 * s_xy.ds - c.dz) / c.cbq;

    out.ds_dr1 = {-out.ds_dr[0], -out.ds_dr[1], -out.ds_dr[2], 0., 0., 0.};

#if T2DS_DEBUG
    Logger::Debug(__FUNCTION__, "(end) deriv.ds_dr = {}", out.ds_dr);
    Logger::Debug(__FUNCTION__, "(end) deriv.ds_dr1 = {}", out.ds_dr1);
#endif

    return out;
}

}  // namespace T2DS::Seeder::HelixVertex
