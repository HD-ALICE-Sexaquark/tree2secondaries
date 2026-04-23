// Compare this file with `original/TransportBz.txt`.

#include <cmath>

#include "Legacy/LegacyFitter.hxx"

namespace Tree2Secondaries::Legacy {
// clang-format off

// Transports the parameters and their covariance matrix of the current particle assuming constant homogeneous
// magnetic field Bz on the length defined by the transport parameter dS = l/p, where l is the signed distance and p is
// the momentum of the current particle. The obtained parameters and covariance matrix are stored to the arrays P and
// C respectively. P and C can be set to the parameters fP and covariance matrix fC of the current particle. In this
// case the particle parameters will be modified. Dependence of the transport parameter dS on the state vector of the
// current particle is taken into account in the covariance matrix using partial derivatives dsdr = d(dS)/d(fP). If
// a pointer to F is initialised the transport jacobian F = d(fP new)/d(fP old) is stored.
// Since dS can depend on the state vector r1 of other particle or vertex, the corelation matrix
// F1 = d(fP new)/d(r1) can be optionally calculated if a pointer F1 is provided.
// Parameters F and F1 should be either both initialised or both set to null pointer.
// \param[in] Bz - z-component of the constant homogeneous magnetic field Bz
// \param[in] dS - transport parameter which defines the distance to which particle should be transported
// \param[in] dsdr[6] = ds/dr - partial derivatives of the parameter dS over the state vector of the current particle
// \param[out] P[8] - array, where transported parameters should be stored
// \param[out] C[36] - array, where transported covariance matrix (8x8) should be stored in the lower triangular form
// \param[in] dsdr1[6] = ds/dr - partial derivatives of the parameter dS over the state vector of another particle
// or vertex
// \param[out] F[36] - optional parameter, transport jacobian, 6x6 matrix F = d(fP new)/d(fP old)
// \param[out] F1[36] - optional parameter, corelation 6x6 matrix betweeen the current particle and particle or vertex
// with the state vector r1, to which the current particle is being transported, F1 = d(fP new)/d(r1)
void TransportBz( const Particle& part, double Bz, double dS, const double* dsdr, double P[], double C[], double* dsdr1, double* F, double* F1 )
{

  const double kCLight = 0.000299792458;
  Bz = Bz*part.fQ*kCLight;
  double bs= Bz*dS;
  double s = std::sin(bs), c = std::cos(bs);
  double sB, cB;
//   if( fabs(bs)>1.e-10){
    sB= s/Bz;
    cB= (1-c)/Bz;
//   }else{
    // const double kOvSqr6 = 1./sqrt(6.);
    // sB = (1.-bs*kOvSqr6)*(1.+bs*kOvSqr6)*dS;
    // cB = .5*sB*bs;
//   }

  double px = part.fP[3];
  double py = part.fP[4];
  double pz = part.fP[5];

  P[0] = part.fP[0] + sB*px + cB*py;
  P[1] = part.fP[1] - cB*px + sB*py;
  P[2] = part.fP[2] +  dS*pz;
  P[3] =          c*px + s*py;
  P[4] =         -s*px + c*py;
  P[5] = part.fP[5];
  P[6] = part.fP[6];
  P[7] = part.fP[7];

  double mJ[8][8];
  for( int i=0; i<8; i++ ) for( int j=0; j<8; j++) mJ[i][j]=0;

  for(int i=0; i<8; i++) mJ[i][i]=1;
  mJ[0][3] =  sB; mJ[0][4] = cB;
  mJ[1][3] = -cB; mJ[1][4] = sB;
  mJ[2][5] = dS;
  mJ[3][3] =  c; mJ[3][4] = s;
  mJ[4][3] = -s; mJ[4][4] = c;


  double mJds[6][6];
  for( int i=0; i<6; i++ ) for( int j=0; j<6; j++) mJds[i][j]=0;
  mJds[0][3] =  c; mJds[0][4] = s;
  mJds[1][3] = -s; mJds[1][4] = c;
  mJds[2][5] = 1;
  mJds[3][3] = -Bz*s; mJds[3][4] =  Bz*c;
  mJds[4][3] = -Bz*c; mJds[4][4] = -Bz*s;

  for(int i1=0; i1<6; i1++)
    for(int i2=0; i2<6; i2++)
      mJ[i1][i2] += mJds[i1][3]*px*dsdr[i2] + mJds[i1][4]*py*dsdr[i2] + mJds[i1][5]*pz*dsdr[i2];

  MultQSQt( mJ[0], part.fC, C, 8);

  if(F)
  {
    for(int i=0; i<6; i++)
      for(int j=0; j<6; j++)
        F[i*6+j] = mJ[i][j];

    for(int i1=0; i1<6; i1++)
      for(int i2=0; i2<6; i2++)
        F1[i1*6 + i2] = mJds[i1][3]*px*dsdr1[i2] + mJds[i1][4]*py*dsdr1[i2] + mJds[i1][5]*pz*dsdr1[i2];
  }
}

// clang-format on
}  // namespace Tree2Secondaries::Legacy
