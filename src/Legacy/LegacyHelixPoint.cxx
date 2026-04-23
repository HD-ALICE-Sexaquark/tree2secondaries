// Compare with: `original/GetDStoPointBz.txt`

#include <cmath>

#include "Legacy/LegacyHelixPoint.hxx"
#include "Legacy/LegacyParticle.hxx"

namespace Tree2Secondaries::Legacy::HelixPoint {
// clang-format off

// Returns dS = l/p parameter, where
// 1) l - signed distance to the DCA point with the input xyz point;
// 2) p - momentum of the particle;
// under the assumption of the constant homogeneous field Bz.
// Also calculates partial derivatives dsdr of the parameter dS over the state vector of the current particle.
// \param[in] B - magnetic field Bz
// \param[in] xyz[3] - point, to which particle should be transported
// \param[out] dsdr[6] = ds/dr partial derivatives of the parameter dS over the state vector of the current particle
double GetDStoPointBz(double B, const Particle& part, const double xyz[3], double dsdr[6])
{

  // if(!param)
    // param = fP;

  double x  = part.X();
  double y  = part.Y();
  double z  = part.Z();
  double px = part.Px();
  double py = part.Py();
  double pz = part.Pz();

  double kCLight = 0.000299792458;
  double bq = B*part.Charge()*kCLight;
  double pt2 = px*px + py*py;
  double p2 = pt2 + pz*pz;

  double dx = xyz[0] - x;
  double dy = xyz[1] - y;
  double dz = xyz[2] - z;
  double a = dx*px+dy*py;
  double dS(0.);

  double abq = bq*a;

  double LocalSmall = 1.e-8;
  bool mask = ( std::abs(bq)<LocalSmall );
  if(mask && p2>1.e-4)
  {
    dS = (a + dz*pz)/p2;

    dsdr[0] = -px/p2;
    dsdr[1] = -py/p2;
    dsdr[2] = -pz/p2;
    dsdr[3] = (dx*p2 - 2.* px *(a + dz *pz))/(p2*p2);
    dsdr[4] = (dy*p2 - 2.* py *(a + dz *pz))/(p2*p2);
    dsdr[5] = (dz*p2 - 2.* pz *(a + dz *pz))/(p2*p2);
  }
  if(mask)
  {
    return dS;
  }

  dS = std::atan2( abq, pt2 + bq*(dy*px -dx*py) )/bq;

  double bs= bq*dS;

  double s = std::sin(bs), c = std::cos(bs);

  if(std::abs(bq) < LocalSmall)
    bq = LocalSmall;
  double bbq = bq*(dx*py - dy*px) - pt2;

  double den = (abq*abq + bbq*bbq);
  den = den < LocalSmall ? LocalSmall : den;

  dsdr[0] = (px*bbq - py*abq)/den;
  dsdr[1] = (px*abq + py*bbq)/den;
  dsdr[2] = 0;
  dsdr[3] = -(dx*bbq + dy*abq + 2.*px*a)/den;
  dsdr[4] = (dx*abq - dy*bbq - 2.*py*a)/den;
  dsdr[5] = 0;

  double sz(0.);
  double cCoeff =  (bbq*c - abq*s) - pz*pz ;
  if(std::abs(cCoeff) > LocalSmall)
    sz = (dS*pz - dz)*pz / cCoeff;

  double dcdr[6] = {0.};
  dcdr[0] = -bq*py*c - bbq*s*bq*dsdr[0] + px*bq*s - abq*c*bq*dsdr[0];
  dcdr[1] =  bq*px*c - bbq*s*bq*dsdr[1] + py*bq*s - abq*c*bq*dsdr[1];
  dcdr[3] = (-bq*dy-2*px)*c - bbq*s*bq*dsdr[3] - dx*bq*s - abq*c*bq*dsdr[3];
  dcdr[4] = ( bq*dx-2*py)*c - bbq*s*bq*dsdr[4] - dy*bq*s - abq*c*bq*dsdr[4];
  dcdr[5] = -2*pz;

  for(int iP=0; iP<6; iP++)
    dsdr[iP] += pz*pz/cCoeff*dsdr[iP] - sz/cCoeff*dcdr[iP];
  dsdr[2] += pz/cCoeff;
  dsdr[5] += (2.*pz*dS - dz)/cCoeff;

  dS += sz;

  bs= bq*dS;
  s = std::sin(bs), c = std::cos(bs);

  double sB, cB;
  //  double kOvSqr6 = 1./std::sqrt(double(6.));

  // if(LocalSmall < std::abs(bs))
  // {
    sB = s/bq;
    cB = (1.-c)/bq;
  // }
  // else
  // {
    // sB = (1.f-bs*kOvSqr6)*(1.f+bs*kOvSqr6)*dS;
    // cB = .5f*sB*bs;
  // }

  double p[5];
  p[0] = x + sB*px + cB*py;
  p[1] = y - cB*px + sB*py;
  p[2] = z +  dS*pz;
  p[3] =          c*px + s*py;
  p[4] =         -s*px + c*py;

  dx = xyz[0] - p[0];
  dy = xyz[1] - p[1];
  dz = xyz[2] - p[2];
  a = dx*p[3]+dy*p[4] + dz*pz;
  abq = bq*a;

  dS += std::atan2( abq, p2 + bq*(dy*p[3] -dx*p[4]) )/bq;

  return dS;
}

// clang-format on
}  // namespace Tree2Secondaries::Legacy::HelixPoint
