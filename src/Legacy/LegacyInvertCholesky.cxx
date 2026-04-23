// Compare with: `original/InvertCholetsky3.txt`

#include <cmath>

#include "Legacy/LegacyFitter.hxx"

namespace Tree2Secondaries::Legacy {
// clang-format off

// Inverts symmetric 3x3 matrix a using modified Choletsky decomposition. The result is stored to the same matrix a.
// \param[in,out] a - 3x3 symmetric matrix
void InvertCholetsky3(double a[6])
{

  double d[3], uud, u[3][3];
  for(int i=0; i<3; i++)
  {
    d[i]=0;
    for(int j=0; j<3; j++)
      u[i][j]=0;
  }

  for(int i=0; i<3 ; i++)
  {
    uud=0;
    for(int j=0; j<i; j++)
      uud += u[j][i]*u[j][i]*d[j];
    uud = a[i*(i+3)/2] - uud;

    if(std::abs(uud)<1.e-12) uud = 1.e-12;

    d[i] = uud/std::abs(uud);
    u[i][i] = std::sqrt(std::abs(uud));

    for(int j=i+1; j<3; j++)
    {
      uud = 0;
      for(int k=0; k<i; k++)
        uud += u[k][i]*u[k][j]*d[k];
      uud = a[j*(j+1)/2+i] - uud;
      u[i][j] = d[i]/u[i][i]*uud;
    }
  }

  double u1[3];

  for(int i=0; i<3; i++)
  {
    u1[i] = u[i][i];
    u[i][i] = 1/u[i][i];
  }
  for(int i=0; i<2; i++)
  {
    u[i][i+1] = - u[i][i+1]*u[i][i]*u[i+1][i+1];
  }
  for(int i=0; i<1; i++)
  {
    u[i][i+2] = u[i][i+1]*u1[i+1]*u[i+1][i+2]-u[i][i+2]*u[i][i]*u[i+2][i+2];
  }

  for(int i=0; i<3; i++)
    a[i+3] = u[i][2]*u[2][2]*d[2];
  for(int i=0; i<2; i++)
    a[i+1] = u[i][1]*u[1][1]*d[1] + u[i][2]*u[1][2]*d[2];
  a[0] = u[0][0]*u[0][0]*d[0] + u[0][1]*u[0][1]*d[1] + u[0][2]*u[0][2]*d[2];
}

// clang-format on
}  // namespace Tree2Secondaries::Legacy
