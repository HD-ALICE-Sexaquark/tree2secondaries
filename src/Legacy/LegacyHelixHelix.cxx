// Compare with: `original/GetDStoParticleBz.txt`

#include <cmath>

#include "Legacy/LegacyHelixHelix.hxx"
#include "Legacy/LegacyParticle.hxx"

namespace Tree2Secondaries::Legacy::HelixHelix {
// clang-format off

// Calculates dS = l/p parameters for two particles, where
// 1) l - signed distance to the DCA point with the other particle;
// 2) p - momentum of the particle;
// under the assumption of the constant homogeneous field Bz. dS[0] is the transport parameter for the current particle,
// dS[1] - for the particle "p".
// Also calculates partial derivatives dsdr of the parameters dS[0] and dS[1] over the state vectors of the particles:
// 1) dsdr[0][6] = d(dS[0])/d(param1);
// 2) dsdr[1][6] = d(dS[0])/d(param2);
// 3) dsdr[2][6] = d(dS[1])/d(param1);
// 4) dsdr[3][6] = d(dS[1])/d(param2);
// where param1 are parameters of the current particle (if the pointer is not provided it is initialised with fP) and
// param2 are parameters of the second particle "p" (if the pointer is not provided it is initialised with p.fP). Parameters
// param1 and param2 should be either provided both or both set to null pointers.
// \param[in] Bz - magnetic field Bz
// \param[in] p - second particle
// \param[out] dS[2] - transport parameters dS for the current particle (dS[0]) and the second particle "p" (dS[1])
// \param[out] dsdr[4][6] - partial derivatives of the parameters dS[0] and dS[1] over the state vectors of the both particles
void GetDStoParticleBz(double Bz, const Particle& p1, const Particle& p2, double dS[2], double dsdr[4][6]) {

//   if(!param1)
//   {
    // param1 = fP;
    // param2 = p.fP;
//   }

  // Get dS to another particle for Bz field
//   double kOvSqr6 = 1./std::sqrt(6.);
  double kCLight = 0.000299792458;

  //in XY plane
  //first root
  double bq1 = Bz*p1.Charge()*kCLight;
  double bq2 = Bz*p2.Charge()*kCLight;

  bool isStraight1 = std::abs(bq1) < 1.e-8;
  bool isStraight2 = std::abs(bq2) < 1.e-8;

//   if( isStraight1 && isStraight2 )
//   {
    // GetDStoParticleLine(p, dS, dsdr);
    // return;
//   }

  double px1 = p1.Px();
  double py1 = p1.Py();
  double pz1 = p1.Pz();

  double px2 = p2.Px();
  double py2 = p2.Py();
  double pz2 = p2.Pz();

  double pt12 = px1*px1 + py1*py1;
  double pt22 = px2*px2 + py2*py2;

  double x01 = p1.X();
  double y01 = p1.Y();
  double z01 = p1.Z();

  double x02 = p2.X();
  double y02 = p2.Y();
  double z02 = p2.Z();

  double dS1[2] = {0.}, dS2[2]={0.};

  double dx0 = (x01 - x02);
  double dy0 = (y01 - y02);
  double dr02 = dx0*dx0 + dy0*dy0;
  double drp1  = dx0*px1 + dy0*py1;
  double dxyp1 = dx0*py1 - dy0*px1;
  double drp2  = dx0*px2 + dy0*py2;
  double dxyp2 = dx0*py2 - dy0*px2;
  double p1p2 = px1*px2 + py1*py2;
  double dp1p2 = px1*py2 - px2*py1;

  double k11 = (bq2*drp1 - dp1p2);
  double k21 = (bq1*(bq2*dxyp1 - p1p2) + bq2*pt12);
  double k12 = ((bq1*drp2 - dp1p2));
  double k22 = (bq2*(bq1*dxyp2 + p1p2) - bq1*pt22);

  double kp = (dxyp1*bq2 - dxyp2*bq1 - p1p2);
  double kd = dr02/2.*bq1*bq2 + kp;
  double c1 = -(bq1*kd + pt12*bq2);
  double c2 = bq2*kd + pt22*bq1;

  double d1 = pt12*pt22 - kd*kd;
  if(d1<0.)
    d1 = 0.;
  d1 = std::sqrt( d1 );
  double d2 = pt12*pt22 - kd*kd;
  if(d2<0.)
    d2 = 0.;
  d2 = std::sqrt( d2 );

  // find two points of closest approach in XY plane

  double dS1dR1[2][6];
  double dS2dR2[2][6];

  double dS1dR2[2][6];
  double dS2dR1[2][6];

  double dk11dr1[6] = {bq2*px1, bq2*py1, 0, bq2*dx0 - py2, bq2*dy0 + px2, 0};
  double dk11dr2[6] = {-bq2*px1, -bq2*py1, 0, py1, -px1, 0};
  double dk12dr1[6] = {bq1*px2, bq1*py2, 0, -py2, px2, 0};
  double dk12dr2[6] = {-bq1*px2, -bq1*py2, 0, bq1*dx0 + py1, bq1*dy0 - px1, 0};
  double dk21dr1[6] = {bq1*bq2*py1, -bq1*bq2*px1, 0, 2*bq2*px1 + bq1*(-(bq2*dy0) - px2), 2*bq2*py1 + bq1*(bq2*dx0 - py2), 0};
  double dk21dr2[6] = {-(bq1*bq2*py1), bq1*bq2*px1, 0, -(bq1*px1), -(bq1*py1), 0};
  double dk22dr1[6] = {bq1*bq2*py2, -(bq1*bq2*px2), 0, bq2*px2, bq2*py2, 0};
  double dk22dr2[6] = {-(bq1*bq2*py2), bq1*bq2*px2, 0, bq2*(-(bq1*dy0) + px1) - 2*bq1*px2, bq2*(bq1*dx0 + py1) - 2*bq1*py2, 0};

  double dkddr1[6] = {bq1*bq2*dx0 + bq2*py1 - bq1*py2, bq1*bq2*dy0 - bq2*px1 + bq1*px2, 0, -bq2*dy0 - px2, bq2*dx0 - py2, 0};
  double dkddr2[6] = {-bq1*bq2*dx0 - bq2*py1 + bq1*py2, -bq1*bq2*dy0 + bq2*px1 - bq1*px2, 0, bq1*dy0 - px1, -bq1*dx0 - py1, 0};

  double dc1dr1[6] = {-(bq1*(bq1*bq2*dx0 + bq2*py1 - bq1*py2)), -(bq1*(bq1*bq2*dy0 - bq2*px1 + bq1*px2)), 0, -2*bq2*px1 - bq1*(-(bq2*dy0) - px2), -2*bq2*py1 - bq1*(bq2*dx0 - py2), 0};
  double dc1dr2[6] = {-(bq1*(-(bq1*bq2*dx0) - bq2*py1 + bq1*py2)), -(bq1*(-(bq1*bq2*dy0) + bq2*px1 - bq1*px2)), 0, -(bq1*(bq1*dy0 - px1)), -(bq1*(-(bq1*dx0) - py1)), 0};

  double dc2dr1[6] = {bq2*(bq1*bq2*dx0 + bq2*py1 - bq1*py2), bq2*(bq1*bq2*dy0 - bq2*px1 + bq1*px2), 0, bq2*(-(bq2*dy0) - px2), bq2*(bq2*dx0 - py2), 0};
  double dc2dr2[6] = {bq2*(-(bq1*bq2*dx0) - bq2*py1 + bq1*py2), bq2*(-(bq1*bq2*dy0) + bq2*px1 - bq1*px2), 0, bq2*(bq1*dy0 - px1) + 2*bq1*px2, bq2*(-(bq1*dx0) - py1) + 2*bq1*py2, 0};

  double dd1dr1[6] = {0,0,0,0,0,0};
  double dd1dr2[6] = {0,0,0,0,0,0};
  if(d1>0)
  {
    for(int i=0; i<6; i++)
    {
      dd1dr1[i] = -kd/d1*dkddr1[i];
      dd1dr2[i] = -kd/d1*dkddr2[i];
    }
    dd1dr1[3] += px1/d1*pt22; dd1dr1[4] += py1/d1*pt22;
    dd1dr2[3] += px2/d1*pt12; dd1dr2[4] += py2/d1*pt12;
  }

  if(!isStraight1)
  {
    dS1[0] = std::atan2( bq1*(k11*c1 + k21*d1), (bq1*k11*d1*bq1 - k21*c1) )/bq1;
    dS1[1] = std::atan2( bq1*(k11*c1 - k21*d1), (-bq1*k11*d1*bq1 - k21*c1) )/bq1;

    double a = bq1*(k11*c1 + k21*d1);
    double b = bq1*k11*d1*bq1 - k21*c1;
    for(int iP=0; iP<6; iP++)
    {
      if(( b*b + a*a ) > 0)
      {
        double dadr1 = bq1*( dk11dr1[iP]*c1 + k11*dc1dr1[iP] + dk21dr1[iP]*d1 + k21*dd1dr1[iP] );
        double dadr2 = bq1*( dk11dr2[iP]*c1 + k11*dc1dr2[iP] + dk21dr2[iP]*d1 + k21*dd1dr2[iP] );
        double dbdr1 = bq1*bq1*( dk11dr1[iP]*d1 + k11*dd1dr1[iP] ) - ( dk21dr1[iP]*c1 + k21*dc1dr1[iP] );
        double dbdr2 = bq1*bq1*( dk11dr2[iP]*d1 + k11*dd1dr2[iP] ) - ( dk21dr2[iP]*c1 + k21*dc1dr2[iP] );

        dS1dR1[0][iP] = 1/bq1 * 1/( b*b + a*a ) * ( dadr1*b - dbdr1*a );
        dS1dR2[0][iP] = 1/bq1 * 1/( b*b + a*a ) * ( dadr2*b - dbdr2*a );
      }
      else
      {
        dS1dR1[0][iP] = 0;
        dS1dR2[0][iP] = 0;
      }
    }

    a = bq1*(k11*c1 - k21*d1);
    b = -bq1*k11*d1*bq1 - k21*c1;
    for(int iP=0; iP<6; iP++)
    {
      if(( b*b + a*a ) > 0)
      {
        double dadr1 = bq1*( dk11dr1[iP]*c1 + k11*dc1dr1[iP] - (dk21dr1[iP]*d1 + k21*dd1dr1[iP]) );
        double dadr2 = bq1*( dk11dr2[iP]*c1 + k11*dc1dr2[iP] - (dk21dr2[iP]*d1 + k21*dd1dr2[iP]) );
        double dbdr1 = -bq1*bq1*( dk11dr1[iP]*d1 + k11*dd1dr1[iP] ) - ( dk21dr1[iP]*c1 + k21*dc1dr1[iP] );
        double dbdr2 = -bq1*bq1*( dk11dr2[iP]*d1 + k11*dd1dr2[iP] ) - ( dk21dr2[iP]*c1 + k21*dc1dr2[iP] );

        dS1dR1[1][iP] = 1/bq1 * 1/( b*b + a*a ) * ( dadr1*b - dbdr1*a );
        dS1dR2[1][iP] = 1/bq1 * 1/( b*b + a*a ) * ( dadr2*b - dbdr2*a );
      }
      else
      {
        dS1dR1[1][iP] = 0;
        dS1dR2[1][iP] = 0;
      }
    }
  }
  if(!isStraight2)
  {
    dS2[0] = std::atan2( (bq2*k12*c2 + k22*d2*bq2), (bq2*k12*d2*bq2 - k22*c2) )/bq2;
    dS2[1] = std::atan2( (bq2*k12*c2 - k22*d2*bq2), (-bq2*k12*d2*bq2 - k22*c2) )/bq2;

    double a = bq2*(k12*c2 + k22*d2);
    double b = bq2*k12*d2*bq2 - k22*c2;
    for(int iP=0; iP<6; iP++)
    {
      if(( b*b + a*a ) > 0)
      {
        double dadr1 = bq2*( dk12dr1[iP]*c2 + k12*dc2dr1[iP] + dk22dr1[iP]*d1 + k22*dd1dr1[iP] );
        double dadr2 = bq2*( dk12dr2[iP]*c2 + k12*dc2dr2[iP] + dk22dr2[iP]*d1 + k22*dd1dr2[iP] );
        double dbdr1 = bq2*bq2*( dk12dr1[iP]*d1 + k12*dd1dr1[iP] ) - (dk22dr1[iP]*c2 + k22*dc2dr1[iP]);
        double dbdr2 = bq2*bq2*( dk12dr2[iP]*d1 + k12*dd1dr2[iP] ) - (dk22dr2[iP]*c2 + k22*dc2dr2[iP]);

        dS2dR1[0][iP] = 1/bq2 * 1/( b*b + a*a ) * ( dadr1*b - dbdr1*a );
        dS2dR2[0][iP] = 1/bq2 * 1/( b*b + a*a ) * ( dadr2*b - dbdr2*a );
      }
      else
      {
        dS2dR1[0][iP] = 0;
        dS2dR2[0][iP] = 0;
      }
    }

    a = bq2*(k12*c2 - k22*d2);
    b = -bq2*k12*d2*bq2 - k22*c2;
    for(int iP=0; iP<6; iP++)
    {
      if(( b*b + a*a ) > 0)
      {
        double dadr1 = bq2*( dk12dr1[iP]*c2 + k12*dc2dr1[iP] - (dk22dr1[iP]*d1 + k22*dd1dr1[iP]) );
        double dadr2 = bq2*( dk12dr2[iP]*c2 + k12*dc2dr2[iP] - (dk22dr2[iP]*d1 + k22*dd1dr2[iP]) );
        double dbdr1 = -bq2*bq2*( dk12dr1[iP]*d1 + k12*dd1dr1[iP] ) - (dk22dr1[iP]*c2 + k22*dc2dr1[iP]);
        double dbdr2 = -bq2*bq2*( dk12dr2[iP]*d1 + k12*dd1dr2[iP] ) - (dk22dr2[iP]*c2 + k22*dc2dr2[iP]);

        dS2dR1[1][iP] = 1/bq2 * 1/( b*b + a*a ) * ( dadr1*b - dbdr1*a );
        dS2dR2[1][iP] = 1/bq2 * 1/( b*b + a*a ) * ( dadr2*b - dbdr2*a );
      }
      else
      {
        dS2dR1[1][iP] = 0;
        dS2dR2[1][iP] = 0;
      }
    }
  }
  if(isStraight1 && (pt12>0.) )
  {
    dS1[0] = (k11*c1 + k21*d1)/(- k21*c1);
    dS1[1] = (k11*c1 - k21*d1)/(- k21*c1);

    double a = k11*c1 + k21*d1;
    double b = -k21*c1;

    for(int iP=0; iP<6; iP++)
    {
      if(b*b > 0)
      {
        double dadr1 = ( dk11dr1[iP]*c1 + k11*dc1dr1[iP] + dk21dr1[iP]*d1 + k21*dd1dr1[iP] );
        double dadr2 = ( dk11dr2[iP]*c1 + k11*dc1dr2[iP] + dk21dr2[iP]*d1 + k21*dd1dr2[iP] );
        double dbdr1 = -( dk21dr1[iP]*c1 + k21*dc1dr1[iP] );
        double dbdr2 = -( dk21dr2[iP]*c1 + k21*dc1dr2[iP] );

        dS1dR1[0][iP] = dadr1/b - dbdr1*a/(b*b) ;
        dS1dR2[0][iP] = dadr2/b - dbdr2*a/(b*b) ;
      }
      else
      {
        dS1dR1[0][iP] = 0;
        dS1dR2[0][iP] = 0;
      }
    }

    a = k11*c1 - k21*d1;
    for(int iP=0; iP<6; iP++)
    {
      if(b*b > 0)
      {
        double dadr1 = ( dk11dr1[iP]*c1 + k11*dc1dr1[iP] - dk21dr1[iP]*d1 - k21*dd1dr1[iP] );
        double dadr2 = ( dk11dr2[iP]*c1 + k11*dc1dr2[iP] - dk21dr2[iP]*d1 - k21*dd1dr2[iP] );
        double dbdr1 = -( dk21dr1[iP]*c1 + k21*dc1dr1[iP] );
        double dbdr2 = -( dk21dr2[iP]*c1 + k21*dc1dr2[iP] );

        dS1dR1[1][iP] = dadr1/b - dbdr1*a/(b*b) ;
        dS1dR2[1][iP] = dadr2/b - dbdr2*a/(b*b) ;
      }
      else
      {
        dS1dR1[1][iP] = 0;
        dS1dR2[1][iP] = 0;
      }
    }
  }
  if(isStraight2 && (pt22>0.) )
  {
    dS2[0] = (k12*c2 + k22*d2)/(- k22*c2);
    dS2[1] = (k12*c2 - k22*d2)/(- k22*c2);

    double a = k12*c2 + k22*d1;
    double b = -k22*c2;

    for(int iP=0; iP<6; iP++)
    {
      if(b*b > 0)
      {
        double dadr1 = ( dk12dr1[iP]*c2 + k12*dc2dr1[iP] + dk22dr1[iP]*d1 + k22*dd1dr1[iP] );
        double dadr2 = ( dk12dr2[iP]*c2 + k12*dc2dr2[iP] + dk22dr2[iP]*d1 + k22*dd1dr2[iP] );
        double dbdr1 = -( dk22dr1[iP]*c2 + k22*dc2dr1[iP] );
        double dbdr2 = -( dk22dr2[iP]*c2 + k22*dc2dr2[iP] );

        dS2dR1[0][iP] = dadr1/b - dbdr1*a/(b*b) ;
        dS2dR2[0][iP] = dadr2/b - dbdr2*a/(b*b) ;
      }
      else
      {
        dS2dR1[0][iP] = 0;
        dS2dR2[0][iP] = 0;
      }
    }

    a = k12*c2 - k22*d1;
    for(int iP=0; iP<6; iP++)
    {
      if(b*b > 0)
      {
        double dadr1 = ( dk12dr1[iP]*c2 + k12*dc2dr1[iP] - dk22dr1[iP]*d1 - k22*dd1dr1[iP] );
        double dadr2 = ( dk12dr2[iP]*c2 + k12*dc2dr2[iP] - dk22dr2[iP]*d1 - k22*dd1dr2[iP] );
        double dbdr1 = -( dk22dr1[iP]*c2 + k22*dc2dr1[iP] );
        double dbdr2 = -( dk22dr2[iP]*c2 + k22*dc2dr2[iP] );

        dS2dR1[1][iP] = dadr1/b - dbdr1*a/(b*b) ;
        dS2dR2[1][iP] = dadr2/b - dbdr2*a/(b*b) ;
      }
      else
      {
        dS2dR1[1][iP] = 0;
        dS2dR2[1][iP] = 0;
      }
    }
  }

  //select a point which is close to the primary vertex (with the smallest r)

  double dr2[2];
  for(int iP = 0; iP<2; iP++)
  {
    double bs1 = bq1*dS1[iP];
    double bs2 = bq2*dS2[iP];
    double sss = std::sin(bs1), ccc = std::cos(bs1);

    // const bool& bs1Big = std::abs(bs1) > 1.e-8;
    // const bool& bs2Big = std::abs(bs2) > 1.e-8;

    double sB(0.), cB(0.);
    // if(bs1Big)
    // {
      sB = sss/bq1;
      cB = (1.-ccc)/bq1;
    // }
    // else
    // {
    //   sB = ((1.-bs1*kOvSqr6)*(1.+bs1*kOvSqr6)*dS1[iP]);
    //   cB = .5*sB*bs1;
    // }

    double x1 = x01 + sB*px1 + cB*py1;
    double y1 = y01 - cB*px1 + sB*py1;
    double z1 = z01 + dS1[iP]*pz1;

    sss = std::sin(bs2), ccc = std::cos(bs2);

    // if(bs2Big)
    // {
      sB = sss/bq2;
      cB = (1.-ccc)/bq2;
    // }
    // else
    // {
    //   sB = ((1.-bs2*kOvSqr6)*(1.+bs2*kOvSqr6)*dS2[iP]);
    //   cB = .5*sB*bs2;
    // }

    double x2 = x02 + sB*px2 + cB*py2;
    double y2 = y02 - cB*px2 + sB*py2;
    double z2 = z02 + dS2[iP]*pz2;

    double dx = (x1-x2);
    double dy = (y1-y2);
    double dz = (z1-z2);

    dr2[iP] = dx*dx + dy*dy + dz*dz;
  }

  const bool isFirstRoot = dr2[0] < dr2[1];
  if(isFirstRoot)
  {
    dS[0]  = dS1[0];
    dS[1] = dS2[0];

    for(int iP=0; iP<6; iP++)
    {
      dsdr[0][iP] = dS1dR1[0][iP];
      dsdr[1][iP] = dS1dR2[0][iP];
      dsdr[2][iP] = dS2dR1[0][iP];
      dsdr[3][iP] = dS2dR2[0][iP];
    }
  }
  else
  {
    dS[0]  = dS1[1];
    dS[1] = dS2[1];

    for(int iP=0; iP<6; iP++)
    {
      dsdr[0][iP] = dS1dR1[1][iP];
      dsdr[1][iP] = dS1dR2[1][iP];
      dsdr[2][iP] = dS2dR1[1][iP];
      dsdr[3][iP] = dS2dR2[1][iP];
    }
  }

  //Line correction
//   {
    double bs1 = bq1*dS[0];
    double bs2 = bq2*dS[1];
    double sss = std::sin(bs1), ccc = std::cos(bs1);

    // const bool& bs1Big = std::abs(bs1) > 1.e-8;
    // const bool& bs2Big = std::abs(bs2) > 1.e-8;

    double sB(0.), cB(0.);
    // if(bs1Big)
    // {
      sB = sss/bq1;
      cB = (1.-ccc)/bq1;
    // }
    // else
    // {
    //   sB = ((1.-bs1*kOvSqr6)*(1.+bs1*kOvSqr6)*dS[0]);
    //   cB = .5*sB*bs1;
    // }

    double x1 = x01 + sB*px1 + cB*py1;
    double y1 = y01 - cB*px1 + sB*py1;
    double z1 = z01 + dS[0]*pz1;
    double ppx1 =  ccc*px1 + sss*py1;
    double ppy1 = -sss*px1 + ccc*py1;
    double ppz1 = pz1;

    double sss1 = sin(bs2), ccc1 = cos(bs2);

    double sB1(0.), cB1(0.);
    // if(bs2Big)
    // {
      sB1 = sss1/bq2;
      cB1 = (1.-ccc1)/bq2;
    // }
    // else
    // {
    //   sB1 = ((1.-bs2*kOvSqr6)*(1.+bs2*kOvSqr6)*dS[1]);
    //   cB1 = .5*sB1*bs2;
    // }

    double x2 = x02 + sB1*px2 + cB1*py2;
    double y2 = y02 - cB1*px2 + sB1*py2;
    double z2 = z02 + dS[1]*pz2;
    double ppx2 =  ccc1*px2 + sss1*py2;
    double ppy2 = -sss1*px2 + ccc1*py2;
    double ppz2 = pz2;

    double p12  = ppx1*ppx1 + ppy1*ppy1 + ppz1*ppz1;
    double p22  = ppx2*ppx2 + ppy2*ppy2 + ppz2*ppz2;
    double lp1p2 = ppx1*ppx2 + ppy1*ppy2 + ppz1*ppz2;

    double dx = (x2 - x1);
    double dy = (y2 - y1);
    double dz = (z2 - z1);

    double ldrp1 = ppx1*dx + ppy1*dy + ppz1*dz;
    double ldrp2 = ppx2*dx + ppy2*dy + ppz2*dz;

    double detp =  lp1p2*lp1p2 - p12*p22;
    if( std::abs(detp)<1.e-4 ) detp = 1;

    //dsdr calculation
    double a1 = ldrp2*lp1p2 - ldrp1*p22;
    double a2 = ldrp2*p12 - ldrp1*lp1p2;
    double lp1p2_ds0 = bq1*( ppx2*ppy1 - ppy2*ppx1);
    double lp1p2_ds1 = bq2*( ppx1*ppy2 - ppy1*ppx2);
    double ldrp1_ds0 = -p12 + bq1*(ppy1*dx - ppx1*dy);
    double ldrp1_ds1 =  lp1p2;
    double ldrp2_ds0 = -lp1p2;
    double ldrp2_ds1 =  p22 + bq2*(ppy2*dx - ppx2*dy);
    double detp_ds0 = 2*lp1p2*lp1p2_ds0;
    double detp_ds1 = 2*lp1p2*lp1p2_ds1;
    double a1_ds0 = ldrp2_ds0*lp1p2 + ldrp2*lp1p2_ds0 - ldrp1_ds0*p22;
    double a1_ds1 = ldrp2_ds1*lp1p2 + ldrp2*lp1p2_ds1 - ldrp1_ds1*p22;
    double a2_ds0 = ldrp2_ds0*p12 - ldrp1_ds0*lp1p2 - ldrp1*lp1p2_ds0;
    double a2_ds1 = ldrp2_ds1*p12 - ldrp1_ds1*lp1p2 - ldrp1*lp1p2_ds1;

    double dsl1ds0 = a1_ds0/detp - a1*detp_ds0/(detp*detp);
    double dsl1ds1 = a1_ds1/detp - a1*detp_ds1/(detp*detp);
    double dsl2ds0 = a2_ds0/detp - a2*detp_ds0/(detp*detp);
    double dsl2ds1 = a2_ds1/detp - a2*detp_ds1/(detp*detp);

    double dsldr[4][6];
    for(int iP=0; iP<6; iP++)
    {
      dsldr[0][iP] = dsl1ds0*dsdr[0][iP] + dsl1ds1*dsdr[2][iP];
      dsldr[1][iP] = dsl1ds0*dsdr[1][iP] + dsl1ds1*dsdr[3][iP];
      dsldr[2][iP] = dsl2ds0*dsdr[0][iP] + dsl2ds1*dsdr[2][iP];
      dsldr[3][iP] = dsl2ds0*dsdr[1][iP] + dsl2ds1*dsdr[3][iP];
    }

    for(int iDS=0; iDS<4; iDS++)
      for(int iP=0; iP<6; iP++)
        dsdr[iDS][iP] += dsldr[iDS][iP];

    double lp1p2_dr0[6] = {0, 0, 0, ccc*ppx2 - ppy2*sss, ccc*ppy2 + ppx2*sss, pz2};
    double lp1p2_dr1[6] = {0, 0, 0, ccc1*ppx1 - ppy1*sss1, ccc1*ppy1 + ppx1*sss1, pz1};
    double ldrp1_dr0[6] = {-ppx1, -ppy1, -pz1,  cB*ppy1 - ppx1*sB + ccc*dx - sss*dy, -cB*ppx1-ppy1*sB + sss*dx + ccc*dy, -dS[0]*pz1 + dz};
    double ldrp1_dr1[6] = { ppx1,  ppy1,  pz1, -cB1*ppy1 + ppx1*sB1, cB1*ppx1 + ppy1*sB1, dS[1]*pz1};
    double ldrp2_dr0[6] = {-ppx2, -ppy2, -pz2, cB*ppy2 - ppx2*sB, -cB*ppx2-ppy2*sB, -dS[0]*pz2};
    double ldrp2_dr1[6] = {ppx2, ppy2, pz2, -cB1*ppy2 + ppx2*sB1 + ccc1*dx- sss1*dy, cB1*ppx2 + ppy2*sB1 + sss1*dx + ccc1*dy, dz + dS[1]*pz2};
    double p12_dr0[6] = {0, 0, 0, 2*px1, 2*py1, 2*pz1};
    double p22_dr1[6] = {0, 0, 0, 2*px2, 2*py2, 2*pz2};
    double a1_dr0[6], a1_dr1[6], a2_dr0[6], a2_dr1[6], detp_dr0[6], detp_dr1[6];
    for(int iP=0; iP<6; iP++)
    {
      a1_dr0[iP] = ldrp2_dr0[iP]*lp1p2 + ldrp2*lp1p2_dr0[iP] - ldrp1_dr0[iP]*p22;
      a1_dr1[iP] = ldrp2_dr1[iP]*lp1p2 + ldrp2*lp1p2_dr1[iP] - ldrp1_dr1[iP]*p22 - ldrp1*p22_dr1[iP];
      a2_dr0[iP] = ldrp2_dr0[iP]*p12 + ldrp2*p12_dr0[iP] - ldrp1_dr0[iP]*lp1p2 - ldrp1*lp1p2_dr0[iP];
      a2_dr1[iP] = ldrp2_dr1[iP]*p12 - ldrp1_dr1[iP]*lp1p2 - ldrp1*lp1p2_dr1[iP];
      detp_dr0[iP] = 2*lp1p2*lp1p2_dr0[iP] - p12_dr0[iP]*p22;
      detp_dr1[iP] = 2*lp1p2*lp1p2_dr1[iP] - p12*p22_dr1[iP];

      dsdr[0][iP] += a1_dr0[iP]/detp - a1*detp_dr0[iP]/(detp*detp);
      dsdr[1][iP] += a1_dr1[iP]/detp - a1*detp_dr1[iP]/(detp*detp);
      dsdr[2][iP] += a2_dr0[iP]/detp - a2*detp_dr0[iP]/(detp*detp);
      dsdr[3][iP] += a2_dr1[iP]/detp - a2*detp_dr1[iP]/(detp*detp);
    }

    dS[0] += (ldrp2*lp1p2 - ldrp1*p22) /detp;
    dS[1] += (ldrp2*p12 - ldrp1*lp1p2)/detp;
//   }
}

// clang-format on
}  // namespace Tree2Secondaries::Legacy::HelixHelix
