#pragma once

#include <array>

#include "App/Logger.hxx"
#include "Legacy/LegacyParticle.hxx"
#include "Seeder/BaseSeeder.hxx"
#include "View/Reconstructed/ViewRecV0.hxx"

namespace Tree2Secondaries::Legacy::LinePoint {

// Formatted function with minimal changes.
double GetDStoPointLine(const Particle& n, const double xyz[3], double dsdr[6]);

// Interface with `View::Reconstructed::V0` + `Seeder::Seed` + `Seeder::Deriv`
inline std::tuple<Seeder::Seed, Seeder::Deriv> FullPCAs(const View::Rec::V0& n, const double xyz[3]) {

    Particle p = Particle::FromV0(n);

    double dsdr[6]{};

    Seeder::Seed seed;
    seed.ds = GetDStoPointLine(p, xyz, dsdr);

#if T2S_DEBUG
    Logger::Debug(__FUNCTION__, "seed.ds = {:13.6e}", seed.ds);
#endif

    Seeder::Deriv deriv;

    deriv.ds_dr = std::to_array(dsdr);

#if T2S_DEBUG
    Logger::Debug(__FUNCTION__, "deriv.ds_dr = {}", deriv.ds_dr);
#endif

    return {seed, deriv};
}

}  // namespace Tree2Secondaries::Legacy::LinePoint
