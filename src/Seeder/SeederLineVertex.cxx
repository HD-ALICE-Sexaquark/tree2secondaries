#include <array>
#include <cmath>

#include "common/Constants.hpp"

#include "Seeder/BaseSeeder.hxx"
#if T2DS_DEBUG
#include "App/Logger.hxx"
#endif

#include "Seeder/SeederLineVertex.hxx"

namespace T2DS::Seeder::LineVertex {

// First phase. Find the point of closest approach (PCA) of a neutral particle w.r.t. a vertex, assuming it transports as a neutral particle.
// Arguments:
// - `n`     -- [input] neutral particle
// - `v`     -- [input] arbitrary vertex
// - `cache` -- [output,optional]
// Return: (packed in a single `Seed` struct)
// - `pca.xyz`, `pca.mom`              -- position and momentum at the PCA
// - `ds`                              -- transport parameter needed to reach PCA
// - `theta`, `sin`, `cos`, `sB`, `cB` -- cached quantities
Seed FastPCA(double x0, double y0, double z0, double px, double py, double pz,  //
             const std::array<double, 3>& vtx, Cache* cache) {

    // cache //

    Cache local;
    Cache& c = cache != nullptr ? *cache : local;

    c.px0 = px;
    c.py0 = py;
    c.pz0 = pz;

    c.dx = vtx[0] - x0;
    c.dy = vtx[1] - y0;
    c.dz = vtx[2] - z0;

    c.p2 = c.px0 * c.px0 + c.py0 * c.py0 + c.pz0 * c.pz0;
    c.a = c.px0 * c.dx + c.py0 * c.dy + c.pz0 * c.dz;

    // particle's seed //

    Seed seed;

    if (c.p2 < Common::AbsAlmostZero) return seed;  // protection

    seed.ds = c.a / c.p2;

    seed.theta = 0.;
    seed.sin = 0.;
    seed.cos = 1.;
    seed.sB = seed.ds;
    seed.cB = 0.;

    seed.pca.xyz = {x0 + c.px0 * seed.ds, y0 + c.py0 * seed.ds, z0 + c.pz0 * seed.ds};
    seed.pca.mom = {c.px0, c.py0, c.pz0};

#if T2DS_DEBUG
    Logger::Debug(__FUNCTION__, "seed.ds = {:13.6e}", seed.ds);
    Logger::Debug(__FUNCTION__, "seed.(x,y,z) = {}", seed.pca.xyz);
    Logger::Debug(__FUNCTION__, "seed.(px,py,pz) = {}", seed.pca.mom);
#endif

    return seed;
}

// Second phase. Compute derivatives of `FastPCA`.
// Argument:
// - `c` -- [input] cache filled in `FastPCA`
// Return: (packed in a single `Deriv` struct)
// - `ds_dr`  -- partial derivatives of current particle's ds w.r.t. current particle's state parameters = d(ds1)/dr1
// - `ds_dr1` -- partial derivatives of current particle's ds w.r.t. other particle's state parameters = d(ds1)/dr2, d(ds2)/dr1
Deriv ComputeDerivatives(const Cache& c) {

    Deriv out;

    if (c.p2 < Common::AbsAlmostZero) return out;  // protection

    double p4 = c.p2 * c.p2;

    out.ds_dr = {-c.px0 / c.p2,
                 -c.py0 / c.p2,
                 -c.pz0 / c.p2,
                 (c.dx * c.p2 - 2. * c.px0 * c.a) / p4,
                 (c.dy * c.p2 - 2. * c.py0 * c.a) / p4,
                 (c.dz * c.p2 - 2. * c.pz0 * c.a) / p4};
    out.ds_dr1 = {-out.ds_dr[0], -out.ds_dr[1], -out.ds_dr[2], 0., 0., 0.};

#if T2DS_DEBUG
    Logger::Debug(__FUNCTION__, "deriv.ds_dr = {}", out.ds_dr);
    Logger::Debug(__FUNCTION__, "deriv.ds_dr1 = {}", out.ds_dr1);
#endif

    return out;
}

}  // namespace T2DS::Seeder::LineVertex
