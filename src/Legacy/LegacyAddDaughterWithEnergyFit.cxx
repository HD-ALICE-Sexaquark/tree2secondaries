// Compare with: `original/AddDaughterWithEnergyFit.txt`

#include "Legacy/LegacyFitter.hxx"
#include "Legacy/LegacyParticle.hxx"

namespace Tree2Secondaries::Legacy {
// clang-format off

// Adds daughter to the current particle. Uses simplifyed fast mathematics which consideres momentum
// and energy as independent variables and thus ignores constraint on the fixed mass.
// In this case the mass of the daughter particle can be corrupted when the constructed vertex
// is added as the measurement and the mass of the output short-lived particle can become
// unphysical - smaller then the threshold.
// \param[in] Daughter - the daughter particle
void AddDaughterWithEnergyFit( Particle &part, const Particle& Daughter, double bz )
{

  // int maxIter = 1;

  // for( int iter=0; iter<maxIter; iter++ ){

    double m[8], mV[36];

    double D[3][3];
    if(! GetMeasurement(part, Daughter, m, mV, D, bz) )
      return;

    double mS[6]= { part.fC[0]+mV[0],
                   part.fC[1]+mV[1], part.fC[2]+mV[2],
                   part.fC[3]+mV[3], part.fC[4]+mV[4], part.fC[5]+mV[5] };

    InvertCholetsky3(mS);

    //* Residual (measured - estimated)

    double zeta[3] = { m[0]-part.fP[0], m[1]-part.fP[1], m[2]-part.fP[2] };

    double dChi2 = (mS[0]*zeta[0] + mS[1]*zeta[1] + mS[3]*zeta[2])*zeta[0]
           +      (mS[1]*zeta[0] + mS[2]*zeta[1] + mS[4]*zeta[2])*zeta[1]
           +      (mS[3]*zeta[0] + mS[4]*zeta[1] + mS[5]*zeta[2])*zeta[2];
    if (dChi2 > 1e9) return;
//     if(fNDF > 100 && dChi2 > 9) return;

    double K[3][3];
    for(int i=0; i<3; i++)
      for(int j=0; j<3; j++)
      {
        K[i][j] = 0;
        for(int k=0; k<3; k++)
          K[i][j] += part.fC[IJ(i,k)] * mS[IJ(k,j)];
      }

    //* CHt = CH' - D'
    double mCHt0[7], mCHt1[7], mCHt2[7];

    mCHt0[0]=part.fC[ 0] ;       mCHt1[0]=part.fC[ 1] ;       mCHt2[0]=part.fC[ 3] ;
    mCHt0[1]=part.fC[ 1] ;       mCHt1[1]=part.fC[ 2] ;       mCHt2[1]=part.fC[ 4] ;
    mCHt0[2]=part.fC[ 3] ;       mCHt1[2]=part.fC[ 4] ;       mCHt2[2]=part.fC[ 5] ;
    mCHt0[3]=part.fC[ 6]-mV[ 6]; mCHt1[3]=part.fC[ 7]-mV[ 7]; mCHt2[3]=part.fC[ 8]-mV[ 8];
    mCHt0[4]=part.fC[10]-mV[10]; mCHt1[4]=part.fC[11]-mV[11]; mCHt2[4]=part.fC[12]-mV[12];
    mCHt0[5]=part.fC[15]-mV[15]; mCHt1[5]=part.fC[16]-mV[16]; mCHt2[5]=part.fC[17]-mV[17];
    mCHt0[6]=part.fC[21]-mV[21]; mCHt1[6]=part.fC[22]-mV[22]; mCHt2[6]=part.fC[23]-mV[23];

    //* Kalman gain K = mCH'*S

    double k0[7], k1[7], k2[7];

    for(int i=0;i<7;++i){
      k0[i] = mCHt0[i]*mS[0] + mCHt1[i]*mS[1] + mCHt2[i]*mS[3];
      k1[i] = mCHt0[i]*mS[1] + mCHt1[i]*mS[2] + mCHt2[i]*mS[4];
      k2[i] = mCHt0[i]*mS[3] + mCHt1[i]*mS[4] + mCHt2[i]*mS[5];
    }

    //* Add the daughter momentum to the particle momentum

    part.fP[ 3] += m[ 3];
    part.fP[ 4] += m[ 4];
    part.fP[ 5] += m[ 5];
    part.fP[ 6] += m[ 6];

    part.fC[ 9] += mV[ 9];
    part.fC[13] += mV[13];
    part.fC[14] += mV[14];
    part.fC[18] += mV[18];
    part.fC[19] += mV[19];
    part.fC[20] += mV[20];
    part.fC[24] += mV[24];
    part.fC[25] += mV[25];
    part.fC[26] += mV[26];
    part.fC[27] += mV[27];


   //* New estimation of the vertex position r += K*zeta

    for(int i=0;i<7;++i)
      part.fP[i] = part.fP[i] + k0[i]*zeta[0] + k1[i]*zeta[1] + k2[i]*zeta[2];

    //* New covariance matrix C -= K*(mCH')'

    for(int i=0, k=0;i<7;++i){
      for(int j=0;j<=i;++j,++k){
	part.fC[k] = part.fC[k] - (k0[i]*mCHt0[j] + k1[i]*mCHt1[j] + k2[i]*mCHt2[j] );
      }
    }

    double K2[3][3];
    for(int i=0; i<3; i++)
    {
      for(int j=0; j<3; j++)
        K2[i][j] = -K[j][i];
      K2[i][i] += 1;
    }

    double A[3][3];
    for(int i=0; i<3; i++)
      for(int j=0; j<3; j++)
      {
        A[i][j] = 0;
        for(int k=0; k<3; k++)
        {
          A[i][j] += D[i][k] * K2[k][j];
        }
      }

    double M[3][3];
    for(int i=0; i<3; i++)
      for(int j=0; j<3; j++)
      {
        M[i][j] = 0;
        for(int k=0; k<3; k++)
        {
          M[i][j] += K[i][k] * A[k][j];
        }
      }

    part.fC[0] += 2*M[0][0];
    part.fC[1] += M[0][1] + M[1][0];
    part.fC[2] += 2*M[1][1];
    part.fC[3] += M[0][2] + M[2][0];
    part.fC[4] += M[1][2] + M[2][1];
    part.fC[5] += 2*M[2][2];

    //* Calculate Chi^2

    part.fNDF  += 2;
    part.fQ    +=  Daughter.fQ;
    // part.fSFromDecay = 0;
    part.fChi2 += dChi2;

  // }
}

// clang-format on
}  // namespace Tree2Secondaries::Legacy
