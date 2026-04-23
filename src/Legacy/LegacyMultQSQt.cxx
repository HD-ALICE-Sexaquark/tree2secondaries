#include "Legacy/LegacyFitter.hxx"

namespace Tree2Secondaries::Legacy {
// clang-format off

// Matrix multiplication SOut = Q*S*Q^T, where Q - square matrix, S - symmetric matrix.
// \param[in] Q - square matrix
// \param[in] S - input symmetric matrix
// \param[out] SOut - output symmetric matrix
// \param[in] kN - dimensionality of the matrices
void MultQSQt( const double Q[], const double S[], double SOut[], const int kN )
{

  double mA[kN*kN];

  for( int i=0, ij=0; i<kN; i++ ){
    for( int j=0; j<kN; j++, ++ij ){
      mA[ij] = 0 ;
      for( int k=0; k<kN; ++k ) mA[ij]+= S[( k<=i ) ? i*(i+1)/2+k :k*(k+1)/2+i] * Q[ j*kN+k];
    }
  }
  for( int i=0; i<kN; i++ ){
    for( int j=0; j<=i; j++ ){
      int ij = ( j<=i ) ? i*(i+1)/2+j :j*(j+1)/2+i;
      SOut[ij] = 0 ;
      for( int k=0; k<kN; k++ )  SOut[ij] += Q[ i*kN+k ] * mA[ k*kN+j ];
    }
  }
}

// clang-format on
}  // namespace Tree2Secondaries::Legacy
