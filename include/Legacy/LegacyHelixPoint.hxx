#pragma once

#include <array>

#include "App/Logger.hxx"
#include "Legacy/LegacyParticle.hxx"
#include "Seeder/BaseSeeder.hxx"
#include "View/Reconstructed/ViewRecTrack.hxx"

namespace R2DS::Legacy::HelixPoint {

// Formatted function with minimal changes.
double GetDStoPointBz(double B, const Particle& part, const double xyz[3], double dsdr[6]);

// Interface with `View::Reconstructed::Track` + `Seeder::Seed` + `Seeder::Deriv`
inline std::tuple<Seeder::Seed, Seeder::Deriv> FullPCAs(const View::Rec::Track& q, double mass, const double xyz[3], double bz) {

    Particle p = Particle::FromTrack(q, mass);

    double dsdr[6]{};

    Seeder::Seed seed;
    seed.ds = GetDStoPointBz(bz, p, xyz, dsdr);

#if R2DS_DEBUG
    Logger::Debug(__FUNCTION__, "seed.ds = {:13.6e}", seed.ds);
#endif

    Seeder::Deriv deriv;

    deriv.ds_dr = std::to_array(dsdr);

#if R2DS_DEBUG
    Logger::Debug(__FUNCTION__, "deriv.ds_dr = {}", deriv.ds_dr);
#endif

    return {seed, deriv};
}

}  // namespace R2DS::Legacy::HelixPoint
