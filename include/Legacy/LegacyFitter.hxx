#pragma once

#include <cmath>

#include <Math/Point3D.h>
#include <Math/Vector3D.h>

#include "Legacy/LegacyParticle.hxx"
#include "View/Reconstructed/ViewRecTrack.hxx"

namespace Tree2Secondaries::Legacy {

// Main Methods //

void AddDaughterWithEnergyFit(Particle& part, const Particle& Daughter, double bz);
void AddDaughter(Particle& part, const Particle& Daughter, double bz);

bool GetMeasurement(Particle& p1, const Particle& p2, double m[], double V[], double D[3][3], double bz);

void TransportLine(const Particle& part, double dS, const double* dsdr, double P[], double C[], double* dsdr1, double* F, double* F1);
void TransportBz(const Particle& part, double Bz, double dS, const double* dsdr, double P[], double C[], double* dsdr1, double* F, double* F1);

void MultQSQt(const double Q[], const double S[], double SOut[], const int kN);
void InvertCholetsky3(double a[6]);

// Inline Methods //

inline void Transport(const Particle& part, double Bz, double dS, const double* dsdr, double P[], double C[], double* dsdr1, double* F, double* F1) {
    if (part.Charge() != 0) {
        TransportBz(part, Bz, dS, dsdr, P, C, dsdr1, F, F1);
    } else {
        TransportLine(part, dS, dsdr, P, C, dsdr1, F, F1);
    }
}

// -- Interface //

inline Particle Fit(const View::Rec::Track& track_1, const View::Rec::Track& track_2, double mass_1, double mass_2, double bz) {
    const auto kf_1 = Particle::FromTrack(track_1, mass_1);
    const auto kf_2 = Particle::FromTrack(track_2, mass_2);
    Particle out;
    AddDaughter(out, kf_1, bz);
    AddDaughter(out, kf_2, bz);
    return out;
}

}  // namespace Tree2Secondaries::Legacy
