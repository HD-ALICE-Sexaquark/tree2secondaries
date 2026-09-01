#include <array>

#include "common/Constants.hpp"

#if T2DS_DEBUG
#include "Utils/Logger.hxx"
#endif
#include "Seeder/SeederTransport.hxx"
#include "Seeder/SeederTypes.hxx"

#include "Seeder/SeederLineVertex.hxx"

namespace T2DS::Seeder::LineVertex {

// Find the point of closest approach (PCA) of a neutral particle w.r.t. a vertex, assuming it transports as a straight line.
// Arguments:
// - `s0`  -- [input] neutral particle
// - `vtx` -- [input] arbitrary vertex
// - `c`   -- [output] cache, consumed by `ComputeDerivatives`
// Return: (packed in a single `Seed` struct)
// - `pca.xyz`, `pca.mom`              -- position and momentum at the PCA
// - `ds`                              -- transport parameter needed to reach PCA
// - `theta`, `sin`, `cos`, `sB`, `cB` -- cached quantities
Seed FastPCA(const State& s0, const std::array<double, 3>& vtx, Cache& c) {

    // cache //

    c.px0 = s0.px;
    c.py0 = s0.py;
    c.pz0 = s0.pz;

    c.dx = vtx[0] - s0.x;
    c.dy = vtx[1] - s0.y;
    c.dz = vtx[2] - s0.z;

    c.p2 = s0.px * s0.px + s0.py * s0.py + s0.pz * s0.pz;
    c.a = s0.px * c.dx + s0.py * c.dy + s0.pz * c.dz;

    // particle's seed //

    if (c.p2 < Common::AbsAlmostZero) [[unlikely]] {
        return TransportLine(s0, 0.);  // protection; degenerate seed with cos=1.
    }

    Seed out = TransportLine(s0, c.a / c.p2);

#if T2DS_DEBUG
    Logger::Debug(__FUNCTION__, "out.ds = {:13.6e}", out.ds);
    Logger::Debug(__FUNCTION__, "out.(x,y,z) = {}", out.pca.xyz);
    Logger::Debug(__FUNCTION__, "out.(px,py,pz) = {}", out.pca.mom);
#endif

    return out;
}

// Compute derivatives of `FastPCA`.
// Argument:
// - `c` -- [input] cache filled in `FastPCA`
// Return: (packed in a single `Deriv` struct)
// - `ds_dr`  -- partial derivatives of current particle's ds w.r.t. current particle's state parameters = d(ds1)/dr1
// - `ds_dr1` -- partial derivatives of current particle's ds w.r.t. the vertex; i.e., negated position derivatives
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
