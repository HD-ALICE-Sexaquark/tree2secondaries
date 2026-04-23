// Compare with: `original/GetDStoParticleLine.txt`

#include <cmath>

#include "Legacy/LegacyLineLine.hxx"
#include "Legacy/LegacyParticle.hxx"

namespace Tree2Secondaries::Legacy::LineLine {
// clang-format off

// Calculates dS = l/p parameters for two particles, where
// 1) l - signed distance to the DCA point with the other particle;
// 2) p - momentum of the particle;
// under the assumption of the straight line trajectory. Is used for particles with charge 0 or in case of zero magnetic field.
// dS[0] is the transport parameter for the current particle, dS[1] - for the particle "p".
// Also calculates partial derivatives dsdr of the parameters dS[0] and dS[1] over the state vectors of the particles:
// 1) dsdr[0][6] = d(dS[0])/d(param1);
// 2) dsdr[1][6] = d(dS[0])/d(param2);
// 3) dsdr[2][6] = d(dS[1])/d(param1);
// 4) dsdr[3][6] = d(dS[1])/d(param2);
// where param1 are parameters of the current particle fP and
// param2 are parameters of the second particle p.fP.
// \param[in] p - second particle
// \param[out] dS[2] - transport parameters dS for the current particle (dS[0]) and the second particle "p" (dS[1])
// \param[out] dsdr[4][6] - partial derivatives of the parameters dS[0] and dS[1] over the state vectors of the both particles
// void GetDStoParticleLine( const KFParticleBase &p, double dS[2], double dsdr[4][6] )
void GetDStoParticleLine(const Particle& n1, const Particle& n2, double dS[2], double dsdr[4][6])
{

  double p12 = n1.Px()*n1.Px() + n1.Py()*n1.Py() + n1.Pz()*n1.Pz();
  double p22 = n2.Px()*n2.Px() + n2.Py()*n2.Py() + n2.Pz()*n2.Pz();
  double p1p2 = n1.Px()*n2.Px() + n1.Py()*n2.Py() + n1.Pz()*n2.Pz();

  double drp1 = n1.Px()*(n2.X()-n1.X()) + n1.Py()*(n2.Y()-n1.Y()) + n1.Pz()*(n2.Z()-n1.Z());
  double drp2 = n2.Px()*(n2.X()-n1.X()) + n2.Py()*(n2.Y()-n1.Y()) + n2.Pz()*(n2.Z()-n1.Z());

  double detp =  p1p2*p1p2 - p12*p22;
  if( std::abs(detp)<1.e-4 ) detp = 1; //TODO correct!!!

  dS[0]  = (drp2*p1p2 - drp1*p22) /detp;
  dS[1] = (drp2*p12  - drp1*p1p2)/detp;

  double x01 = n1.X();
  double y01 = n1.Y();
  double z01 = n1.Z();
  double px1 = n1.Px();
  double py1 = n1.Py();
  double pz1 = n1.Pz();

  double x02 = n2.X();
  double y02 = n2.Y();
  double z02 = n2.Z();
  double px2 = n2.Px();
  double py2 = n2.Py();
  double pz2 = n2.Pz();

  double drp1_dr1[6] = {-px1, -py1, -pz1, -x01 + x02, -y01 + y02, -z01 + z02};
  double drp1_dr2[6] = {px1, py1, pz1, 0, 0, 0};
  double drp2_dr1[6] = {-px2, -py2, -pz2, 0, 0, 0};
  double drp2_dr2[6] = {px2, py2, pz2, -x01 + x02, -y01 + y02, -z01 + z02};
  double dp1p2_dr1[6] = {0, 0, 0, px2, py2, pz2};
  double dp1p2_dr2[6] = {0, 0, 0, px1, py1, pz1};
  double dp12_dr1[6] = {0, 0, 0, 2*px1, 2*py1, 2*pz1};
  double dp12_dr2[6] = {0, 0, 0, 0, 0, 0};
  double dp22_dr1[6] = {0, 0, 0, 0, 0, 0};
  double dp22_dr2[6] = {0, 0, 0, 2*px2, 2*py2, 2*pz2};
  double ddetp_dr1[6] = {0, 0, 0, -2*p22*px1 + 2*p1p2*px2, -2*p22*py1 + 2*p1p2*py2, -2*p22*pz1 + 2*p1p2*pz2};
  double ddetp_dr2[6] = {0, 0, 0, 2*p1p2*px1 - 2*p12*px2,   2*p1p2*py1 - 2*p12*py2, 2*p1p2*pz1 - 2*p12*pz2};


  double da1_dr1[6], da1_dr2[6], da2_dr1[6], da2_dr2[6];

  double a1 = drp2*p1p2 - drp1*p22;
  double a2 = drp2*p12  - drp1*p1p2;
  for(int i=0; i<6; i++)
  {
    da1_dr1[i] = drp2_dr1[i]*p1p2 + drp2*dp1p2_dr1[i] - drp1_dr1[i]*p22 - drp1*dp22_dr1[i];
    da1_dr2[i] = drp2_dr2[i]*p1p2 + drp2*dp1p2_dr2[i] - drp1_dr2[i]*p22 - drp1*dp22_dr2[i];

    da2_dr1[i] = drp2_dr1[i]*p12 + drp2*dp12_dr1[i] - drp1_dr1[i]*p1p2 - drp1*dp1p2_dr1[i];
    da2_dr2[i] = drp2_dr2[i]*p12 + drp2*dp12_dr2[i] - drp1_dr2[i]*p1p2 - drp1*dp1p2_dr2[i];

    dsdr[0][i] = da1_dr1[i]/detp - a1/(detp*detp)*ddetp_dr1[i];
    dsdr[1][i] = da1_dr2[i]/detp - a1/(detp*detp)*ddetp_dr2[i];

    dsdr[2][i] = da2_dr1[i]/detp - a2/(detp*detp)*ddetp_dr1[i];
    dsdr[3][i] = da2_dr2[i]/detp - a2/(detp*detp)*ddetp_dr2[i];
  }
}

// clang-format off
}
