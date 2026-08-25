#pragma once

#include <cassert>
#include <cmath>
#include <tuple>

#include "common/Constants.hpp"  // needed for `Common::AbsAlmostZero` below
#include "common/Math.hpp"

#include "Seeder/SeederTypes.hxx"

namespace T2DS::Seeder {

// Transport `s0` along a straight line by the path parameter `ds`.
// The trigonometric cache degenerates to theta = 0, sin = 0, cos = 1, sB = ds, cB = 0, which is exactly what
// `KF::Detail::Daughter::PrepareJacobAndCorr` needs to build the q = 0 Jacobian.
[[nodiscard]] inline Seed TransportLine(double x0, double y0, double z0, double px0, double py0, double pz0, double ds) {
    Seed out;

    out.theta = 0.;
    out.sin = 0.;
    out.cos = 1.;
    out.sB = ds;
    out.cB = 0.;
    out.ds = ds;
    out.pca.xyz = {x0 + px0 * ds, y0 + py0 * ds, z0 + pz0 * ds};
    out.pca.mom = {px0, py0, pz0};

    return out;
}
[[nodiscard]] inline Seed TransportLine(const State& s0, double ds) {  //
    return TransportLine(s0.x, s0.y, s0.z, s0.px, s0.py, s0.pz, ds);
}

// Transport `s0` along its helix through the turn angle `theta`, using this parametrization:
//     x' = x + sB * px + cB * py
//     y' = y - cB * px + sB * py
//     z' = z + ds * pz
//     px' =  cos * px + sin * py
//     py' = -sin * px + cos * py
//     pz' = pz
// NOTE: both `theta` and `ds` are taken as arguments, despite theta == bq * ds;
//       this is done to prevent deriving each other again through division;
//       they can be input independently, though, in the functions below
[[nodiscard]] inline Seed TransportHelix(double x0, double y0, double z0, double px0, double py0, double pz0, double theta, double ds, double bq) {
    assert(std::abs(bq) > Common::AbsAlmostZero);
    Seed out;

    out.theta = theta;
    std::tie(out.sin, out.cos) = Common::Math::sincos(theta);
    out.sB = out.sin / bq;
    out.cB = (1. - out.cos) / bq;
    out.ds = ds;
    out.pca.xyz = {x0 + out.sB * px0 + out.cB * py0, y0 - out.cB * px0 + out.sB * py0, z0 + out.ds * pz0};
    out.pca.mom = {out.cos * px0 + out.sin * py0, -out.sin * px0 + out.cos * py0, pz0};

    return out;
}
[[nodiscard]] inline Seed TransportHelix(const State& s0, double theta, double ds, double bq) {
    return TransportHelix(s0.x, s0.y, s0.z, s0.px, s0.py, s0.pz, theta, ds, bq);
}

[[nodiscard]] inline Seed TransportHelixFromTheta(const State& s0, double theta, double bq) {  //
    assert(std::abs(bq) > Common::AbsAlmostZero);
    return TransportHelix(s0, theta, theta / bq, bq);
}
[[nodiscard]] inline Seed TransportHelixFromTheta(double x0, double y0, double z0, double px0, double py0, double pz0, double theta, double bq) {
    assert(std::abs(bq) > Common::AbsAlmostZero);
    return TransportHelix(x0, y0, z0, px0, py0, pz0, theta, theta / bq, bq);
}

[[nodiscard]] inline Seed TransportHelixFromDs(const State& s0, double ds, double bq) {  //
    return TransportHelix(s0, bq * ds, ds, bq);
}
[[nodiscard]] inline Seed TransportHelixFromDs(double x0, double y0, double z0, double px0, double py0, double pz0, double ds, double bq) {
    return TransportHelix(x0, y0, z0, px0, py0, pz0, bq * ds, ds, bq);
}

}  // namespace T2DS::Seeder
