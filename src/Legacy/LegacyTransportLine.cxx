// Compare this file with `original/TransportLine.txt`.

#include "Legacy/LegacyFitter.hxx"
#include "Legacy/LegacyParticle.hxx"

namespace Tree2Secondaries::Legacy {
// clang-format off

//Transports the parameters and their covariance matrix of the current particle assuming the straight line trajectory
// on the length defined by the transport parameter dS = l/p, where l is the signed distance and p is
// the momentum of the current particle. The obtained parameters and covariance matrix are stored to the arrays P and
// C respectively. P and C can be set to the parameters fP and covariance matrix fC of the current particle. In this
// case the particle parameters will be modified. Dependence of the transport parameter dS on the state vector of the
// current particle is taken into account in the covariance matrix using partial derivatives dsdr = d(dS)/d(fP). If
// a pointer to F is initialised the transport jacobian F = d(fP new)/d(fP old) is stored.
// Since dS can depend on the state vector r1 of other particle or vertex, the corelation matrix
// F1 = d(fP new)/d(r1) can be optionally calculated if a pointer F1 is provided.
// Parameters F and F1 should be either both initialised or both set to null pointer.
// \param[in] dS - transport parameter which defines the distance to which particle should be transported
// \param[in] dsdr[6] = ds/dr - partial derivatives of the parameter dS over the state vector of the current particle
// \param[out] P[8] - array, where transported parameters should be stored
// \param[out] C[36] - array, where transported covariance matrix (8x8) should be stored in the lower triangular form
// \param[in] dsdr1[6] = ds/dr - partial derivatives of the parameter dS over the state vector of another particle
// or vertex
// \param[out] F[36] - optional parameter, transport jacobian, 6x6 matrix F = d(fP new)/d(fP old)
// \param[out] F1[36] - optional parameter, corelation 6x6 matrix betweeen the current particle and particle or vertex
// with the state vector r1, to which the current particle is being transported, F1 = d(fP new)/d(r1)
void TransportLine( const Particle& part, double dS, const double* dsdr, double P[], double C[], double* dsdr1, double* F, double* F1 )
{

  double mJ[8][8];
  for( int i=0; i<8; i++ ) for( int j=0; j<8; j++) mJ[i][j]=0;

  mJ[0][0]=1; mJ[0][1]=0; mJ[0][2]=0; mJ[0][3]=dS;  mJ[0][4]=0;  mJ[0][5]=0;
  mJ[1][0]=0; mJ[1][1]=1; mJ[1][2]=0; mJ[1][3]=0;     mJ[1][4]=dS;  mJ[1][5]=0;
  mJ[2][0]=0; mJ[2][1]=0; mJ[2][2]=1; mJ[2][3]=0; mJ[2][4]=0; mJ[2][5]=dS;

  mJ[3][0]=0; mJ[3][1]=0; mJ[3][2]=0; mJ[3][3]=1;   mJ[3][4]=0;  mJ[3][5]=0;
  mJ[4][0]=0; mJ[4][1]=0; mJ[4][2]=0; mJ[4][3]=0;     mJ[4][4]=1;   mJ[4][5]=0;
  mJ[5][0]=0; mJ[5][1]=0; mJ[5][2]=0; mJ[5][3]=0; mJ[5][4]=0; mJ[5][5]=1;
  mJ[6][6] = mJ[7][7] = 1;

  double px = part.fP[3], py = part.fP[4], pz = part.fP[5];

  P[0] = part.fP[0] + dS*part.fP[3];
  P[1] = part.fP[1] + dS*part.fP[4];
  P[2] = part.fP[2] + dS*part.fP[5];
  P[3] = part.fP[3];
  P[4] = part.fP[4];
  P[5] = part.fP[5];
  P[6] = part.fP[6];
  P[7] = part.fP[7];

  double mJds[6][6];
  for( int i=0; i<6; i++ ) for( int j=0; j<6; j++) mJds[i][j]=0;

  mJds[0][3]= 1;
  mJds[1][4]= 1;
  mJds[2][5]= 1;

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
