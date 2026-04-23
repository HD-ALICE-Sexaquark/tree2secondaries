#pragma once

#include <array>

#include "App/Logger.hxx"
#include "Legacy/LegacyParticle.hxx"
#include "Seeder/BaseSeeder.hxx"
#include "View/Reconstructed/ViewRecV0.hxx"

namespace Tree2Secondaries::Legacy::LineLine {

// Formatted function with minimal changes.
void GetDStoParticleLine(const Particle& n1, const Particle& n2, double dS[2], double dsdr[4][6]);

// Interface with `View::Reconstructed::V0` + `Seeder::Seed` + `Seeder::Deriv`
inline std::tuple<Seeder::Seed, Seeder::Seed, Seeder::Deriv, Seeder::Deriv> FullPCAs(const View::Rec::V0& n1, const View::Rec::V0& n2) {

    Particle p1 = Particle::FromV0(n1);
    Particle p2 = Particle::FromV0(n2);

    double dS[2]{};
    double dsdr[4][6]{};
    GetDStoParticleLine(p1, p2, dS, dsdr);

    Seeder::Seed s1;
    Seeder::Seed s2;

    s1.ds = dS[0];
    s2.ds = dS[1];

#if T2S_DEBUG
    Logger::Debug(__FUNCTION__, "seed1.ds = {:13.6e}", s1.ds);
    Logger::Debug(__FUNCTION__, "seed2.ds = {:13.6e}", s2.ds);
#endif

    Seeder::Deriv d1;
    Seeder::Deriv d2;

    d1.ds_dr = std::to_array(dsdr[0]);
    d1.ds_dr1 = std::to_array(dsdr[1]);
    d2.ds_dr1 = std::to_array(dsdr[2]);
    d2.ds_dr = std::to_array(dsdr[3]);

#if T2S_DEBUG
    Logger::Debug(__FUNCTION__, "deriv1.ds_dr = {}", d1.ds_dr);
    Logger::Debug(__FUNCTION__, "deriv1.ds_dr1 = {}", d1.ds_dr1);
    Logger::Debug(__FUNCTION__, "deriv2.ds_dr = {}", d2.ds_dr);
    Logger::Debug(__FUNCTION__, "deriv2.ds_dr1 = {}", d2.ds_dr1);
#endif

    return {s1, s2, d1, d2};
}

}  // namespace Tree2Secondaries::Legacy::LineLine
