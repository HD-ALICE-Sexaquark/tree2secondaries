// Compare with: `original/GetMeasurement.txt`

#include "Legacy/BaseLegacy.hxx"
#include "Legacy/LegacyFitter.hxx"

namespace Tree2Secondaries::Legacy {
// clang-format off

// Obtains the measurements from the current particle and the daughter to be added for the Kalman filter
// mathematics. If these are two first daughters they are transported to the point of the closest approach,
// if the third or higher daughter is added it is transported to the DCA point of the already constructed
// vertex. The correlations are taken into account in the covariance matrices of both measurements,
// the correlation matrix of two measurements is also calculated. Parameters of the current particle are
// modified by this function, the daughter is not changed, its parameters are stored to the output arrays
// after modifications.
// \param[in] daughter - the daughter particle to be added, stays unchanged
// \param[out] m[8] - the output parameters of the daughter particle at the DCA point
// \param[out] V[36] - the output covariance matrix of the daughter parameters, takes into account the correlation
// \param[out] D[3][3] - the correlation matrix between the current and daughter particles
bool GetMeasurement( Particle& p1, const Particle& daughter, double m[], double V[], double D[3][3], double bz )
{

  if(p1.fNDF == -1)
  {
    double ds[2] = {0.,0.};
    double dsdr[4][6];
    double F1[36], F2[36], F3[36], F4[36];
    for(int i1=0; i1<36; i1++)
    {
      F1[i1] = 0;
      F2[i1] = 0;
      F3[i1] = 0;
      F4[i1] = 0;
    }
    GetDStoParticle( p1, daughter, ds, dsdr, bz );

    if( std::abs(ds[0]*p1.fP[5]) > 1000. || std::abs(ds[1]*daughter.fP[5]) > 1000.)
      return false;

    double V0Tmp[36] = {0.};
    double V1Tmp[36] = {0.};

    double C[36];
    for(int iC=0; iC<36; iC++)
      C[iC] = p1.fC[iC];

             Transport(p1, bz, ds[0], dsdr[0], p1.fP, p1.fC, dsdr[1], F1, F2);
    Transport(daughter, bz, ds[1], dsdr[3],  m,  V, dsdr[2], F4, F3);

    MultQSQt(F2, daughter.fC, V0Tmp, 6);
    MultQSQt(F3, C, V1Tmp, 6);

    for(int iC=0; iC<21; iC++)
    {
      p1.fC[iC] += V0Tmp[iC];
      V[iC]  += V1Tmp[iC];
    }

    double C1F1T[6][6];
    for(int i=0; i<6; i++)
      for(int j=0; j<6; j++)
      {
        C1F1T[i][j] = 0;
        for(int k=0; k<6; k++)
        {
          C1F1T[i][j] +=  C[IJ(i,k)] * F1[j*6+k];
        }
      }
    double F3C1F1T[6][6];
    for(int i=0; i<6; i++)
      for(int j=0; j<6; j++)
      {
        F3C1F1T[i][j] = 0;
        for(int k=0; k<6; k++)
        {
          F3C1F1T[i][j] += F3[i*6+k] * C1F1T[k][j];
        }
      }
    double C2F2T[6][6];
    for(int i=0; i<6; i++)
      for(int j=0; j<6; j++)
      {
        C2F2T[i][j] = 0;
        for(int k=0; k<6; k++)
        {
          C2F2T[i][j] +=  daughter.fC[IJ(i,k)] * F2[j*6+k];
        }
      }
    for(int i=0; i<3; i++)
      for(int j=0; j<3; j++)
      {
        D[i][j] = F3C1F1T[i][j];
        for(int k=0; k<6; k++)
        {
          D[i][j] += F4[i*6+k] * C2F2T[k][j];
        }
      }
  }
  else
  {
    double dsdr[6];
    double dS = GetDStoPoint(daughter, p1.fP, dsdr);

    double dsdp[6] = {-dsdr[0], -dsdr[1], -dsdr[2], 0, 0, 0};

    double F[36], F1[36];
    for(int i2=0; i2<36; i2++)
    {
      F[i2]  = 0;
      F1[i2] = 0;
    }
    Transport(daughter, bz, dS, dsdr, m, V, dsdp, F, F1);

//     double V1Tmp[36] = {0.};
//     MultQSQt(F1, fC, V1Tmp, 6);

//     for(int iC=0; iC<21; iC++)
//       V[iC] += V1Tmp[iC];

    double VFT[3][6];
    for(int i=0; i<3; i++)
      for(int j=0; j<6; j++)
      {
        VFT[i][j] = 0;
        for(int k=0; k<3; k++)
        {
          VFT[i][j] +=  p1.fC[IJ(i,k)] * F1[j*6+k];
        }
      }

    double FVFT[6][6];
    for(int i=0; i<6; i++)
      for(int j=0; j<6; j++)
      {
        FVFT[i][j] = 0;
        for(int k=0; k<3; k++)
        {
          FVFT[i][j] += F1[i*6+k] * VFT[k][j];
        }
      }

    for(int i=0; i<3; i++)
      for(int j=0; j<3; j++)
      {
        D[i][j] = 0;
        for(int k=0; k<3; k++)
        {
          D[i][j] += p1.fC[IJ(j,k)] * F1[i*6+k];
        }
      }

    V[0] += FVFT[0][0];
    V[1] += FVFT[1][0];
    V[2] += FVFT[1][1];
    V[3] += FVFT[2][0];
    V[4] += FVFT[2][1];
    V[5] += FVFT[2][2];

//     if(fNDF > 100)
//     {
//       double dx = fP[0] - m[0];
//       double dy = fP[1] - m[1];
//       double dz = fP[2] - m[2];
//       double sigmaS = 3.f*sqrt( (dx*dx + dy*dy + dz*dz) / (m[3]*m[3] + m[4]*m[4] + m[5]*m[5]) );
//
//       double h[3] = { m[3]*sigmaS, m[4]*sigmaS, m[5]*sigmaS };
//       V[0]+= h[0]*h[0];
//       V[1]+= h[1]*h[0];
//       V[2]+= h[1]*h[1];
//       V[3]+= h[2]*h[0];
//       V[4]+= h[2]*h[1];
//       V[5]+= h[2]*h[2];
//     }
  }

  return true;
}

// clang-format on
}  // namespace Tree2Secondaries::Legacy
