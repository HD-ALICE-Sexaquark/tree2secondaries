// Compare with: `original/GetDStoPointLine.txt`

#include <cmath>

#include "Legacy/LegacyLinePoint.hxx"
#include "Legacy/LegacyParticle.hxx"

namespace Tree2Secondaries::Legacy::LinePoint {
// clang-format off

// Returns dS = l/p parameter, where
// 1) l - signed distance to the DCA point with the input xyz point
// 2) p - momentum of the particle;
// assuming the straigth line trajectory. Is used for particles with charge 0 or in case of zero magnetic field.
// Also calculates partial derivatives dsdr of the parameter dS over the state vector of the current particle.
// \param[in] xyz[3] - point where particle should be transported
// \param[out] dsdr[6] = ds/dr partial derivatives of the parameter dS over the state vector of the current particle
double GetDStoPointLine(const Particle& n, const double xyz[3], double dsdr[6])
{

  double p2 = n.Px()*n.Px() + n.Py()*n.Py() + n.Pz()*n.Pz();
  if( p2<1.e-4 ) p2 = 1;

  double a = n.Px()*(xyz[0]-n.X()) + n.Py()*(xyz[1]-n.Y()) + n.Pz()*(xyz[2]-n.Z());
  dsdr[0] = -n.Px()/p2;
  dsdr[1] = -n.Py()/p2;
  dsdr[2] = -n.Pz()/p2;
  dsdr[3] = ((xyz[0]-n.X())*p2 - 2.* n.Px()*a)/(p2*p2);
  dsdr[4] = ((xyz[1]-n.Y())*p2 - 2.* n.Py()*a)/(p2*p2);
  dsdr[5] = ((xyz[2]-n.Z())*p2 - 2.* n.Pz()*a)/(p2*p2);

  return a/p2;
}

// clang-format on
}  // namespace Tree2Secondaries::Legacy::LinePoint
