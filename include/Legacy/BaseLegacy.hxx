#pragma once

#include "Legacy/LegacyHelixHelix.hxx"
#include "Legacy/LegacyHelixPoint.hxx"
#include "Legacy/LegacyLineLine.hxx"
#include "Legacy/LegacyLinePoint.hxx"
#include "Legacy/LegacyParticle.hxx"

namespace Tree2Secondaries::Legacy {

inline void GetDStoParticle(const Particle& p1, const Particle& p2, double ds[2], double dsdr[4][6], double bz = 0.) {
    if (p1.Charge() == 0 && p2.Charge() == 0) {
        LineLine::GetDStoParticleLine(p1, p2, ds, dsdr);
    } else {
        HelixHelix::GetDStoParticleBz(bz, p1, p2, ds, dsdr);
    }
}

inline double GetDStoPoint(const Particle& part, const double xyz[3], double dsdr[6], double bz = 0.) {
    if (part.Charge() == 0) {
        return LinePoint::GetDStoPointLine(part, xyz, dsdr);
    }  // early return
    return HelixPoint::GetDStoPointBz(bz, part, xyz, dsdr);
}

}  // namespace Tree2Secondaries::Legacy
