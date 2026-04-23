// Compare this file with `original/Initialize.txt`.

#include <cmath>

#include "Legacy/LegacyParticle.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::Legacy {
// clang-format off

// Initialises the parameters by default:
// 1) all parameters are set to 0;
// 2) all elements of the covariance matrix are set to 0 except Cxx=Cyy=Czz=100;
// 3) Q = 0;
// 4) chi2 is set to 0;
// 5) NDF = -3, since 3 parameters should be fitted: X, Y, Z.
void Particle::Initialize()
{

  for( int i=0; i<8; i++) fP[i] = 0;
  for(int i=0;i<36;++i) fC[i]=0.;
  fC[0] = fC[2] = fC[5] = 100.;
  fC[35] = 1.;
  fNDF  = -3;
  fChi2 =  0.;
  fQ = 0;
  // fSFromDecay = 0;
  // fAtProductionVertex = 0;
  // SumDaughterMass = 0;
  // fMassHypo = -1;
}

// Sets the parameters of the particle:
//
// \param[in] Param[6] = { X, Y, Z, Px, Py, Pz } - position and momentum
// \param[in] Cov[21]  - lower-triangular part of the covariance matrix:
// \verbatim
//           (  0  .  .  .  .  . )
//           (  1  2  .  .  .  . )
// Cov[21] = (  3  4  5  .  .  . )
//           (  6  7  8  9  .  . )
//           ( 10 11 12 13 14  . )
//           ( 15 16 17 18 19 20 )
// \endverbatim
// \param[in] Charge - charge of the particle in elementary charge units
// \param[in] mass - the mass hypothesis
// void Initialize( const float Param[], const float Cov[], int Charge, float Mass )
Particle Particle::FromTrack( const View::Rec::Track &track, double mass )
{

    Particle part;

//   for( int i=0; i<6 ; i++ ) part.fP[i] = Param[i];
    part.fP[0] = track.X();
    part.fP[1] = track.Y();
    part.fP[2] = track.Z();
    part.fP[3] = track.Px();
    part.fP[4] = track.Py();
    part.fP[5] = track.Pz();
//   for( int i=0; i<21; i++ ) part.fC[i] = Cov[i];
    part.fC[IJ(0, 0)] = track.SigmaX2();
    part.fC[IJ(1, 0)] = track.SigmaXY();
    part.fC[IJ(1, 1)] = track.SigmaY2();
    part.fC[IJ(2, 0)] = track.SigmaXZ();
    part.fC[IJ(2, 1)] = track.SigmaYZ();
    part.fC[IJ(2, 2)] = track.SigmaZ2();
    part.fC[IJ(3, 0)] = track.SigmaXPx();
    part.fC[IJ(3, 1)] = track.SigmaYPx();
    part.fC[IJ(3, 2)] = track.SigmaZPx();
    part.fC[IJ(3, 3)] = track.SigmaPx2();
    part.fC[IJ(4, 0)] = track.SigmaXPy();
    part.fC[IJ(4, 1)] = track.SigmaYPy();
    part.fC[IJ(4, 2)] = track.SigmaZPy();
    part.fC[IJ(4, 3)] = track.SigmaPxPy();
    part.fC[IJ(4, 4)] = track.SigmaPy2();
    part.fC[IJ(5, 0)] = track.SigmaXPz();
    part.fC[IJ(5, 1)] = track.SigmaYPz();
    part.fC[IJ(5, 2)] = track.SigmaZPz();
    part.fC[IJ(5, 3)] = track.SigmaPxPz();
    part.fC[IJ(5, 4)] = track.SigmaPyPz();
    part.fC[IJ(5, 5)] = track.SigmaPz2();

  double energy = std::sqrt( mass*mass + part.fP[3]*part.fP[3] + part.fP[4]*part.fP[4] + part.fP[5]*part.fP[5]);
  part.fP[6] = energy;
  part.fP[7] = 0.;
  part.fQ = track.Charge();
  part.fNDF = 0;
  part.fChi2 = 0.;
//   fAtProductionVertex = 0;
//   fSFromDecay = 0;

  double energyInv = 1./energy;
  double
    h0 = part.fP[3]*energyInv,
    h1 = part.fP[4]*energyInv,
    h2 = part.fP[5]*energyInv;

  part.fC[21] = h0*part.fC[ 6] + h1*part.fC[10] + h2*part.fC[15];
  part.fC[22] = h0*part.fC[ 7] + h1*part.fC[11] + h2*part.fC[16];
  part.fC[23] = h0*part.fC[ 8] + h1*part.fC[12] + h2*part.fC[17];
  part.fC[24] = h0*part.fC[ 9] + h1*part.fC[13] + h2*part.fC[18];
  part.fC[25] = h0*part.fC[13] + h1*part.fC[14] + h2*part.fC[19];
  part.fC[26] = h0*part.fC[18] + h1*part.fC[19] + h2*part.fC[20];
  part.fC[27] = ( h0*h0*part.fC[ 9] + h1*h1*part.fC[14] + h2*h2*part.fC[20]
	     + 2*(h0*h1*part.fC[13] + h0*h2*part.fC[18] + h1*h2*part.fC[19] ) );
  for( int i=28; i<36; i++ ) part.fC[i] = 0;
  part.fC[35] = 1.;

//   SumDaughterMass = Mass;
//   fMassHypo = Mass;

  return part;
}

// clang-format on

// Create a `Legacy::Particle`, by setting `fP`, `fC` and `fQ` from a V0 view.
Particle Particle::FromV0(const View::Rec::V0& v0) {

    Particle out{};

    out.fP[0] = v0.X();
    out.fP[1] = v0.Y();
    out.fP[2] = v0.Z();
    out.fP[3] = v0.Px();
    out.fP[4] = v0.Py();
    out.fP[5] = v0.Pz();
    out.fP[6] = v0.Energy();
    out.fP[7] = 0.;

    out.fC[IJ(0, 0)] = v0.SigmaX2();
    out.fC[IJ(1, 0)] = v0.SigmaXY();
    out.fC[IJ(1, 1)] = v0.SigmaY2();
    out.fC[IJ(2, 0)] = v0.SigmaXZ();
    out.fC[IJ(2, 1)] = v0.SigmaYZ();
    out.fC[IJ(2, 2)] = v0.SigmaZ2();
    out.fC[IJ(3, 0)] = v0.SigmaXPx();
    out.fC[IJ(3, 1)] = v0.SigmaYPx();
    out.fC[IJ(3, 2)] = v0.SigmaZPx();
    out.fC[IJ(3, 3)] = v0.SigmaPx2();
    out.fC[IJ(4, 0)] = v0.SigmaXPy();
    out.fC[IJ(4, 1)] = v0.SigmaYPy();
    out.fC[IJ(4, 2)] = v0.SigmaZPy();
    out.fC[IJ(4, 3)] = v0.SigmaPxPy();
    out.fC[IJ(4, 4)] = v0.SigmaPy2();
    out.fC[IJ(5, 0)] = v0.SigmaXPz();
    out.fC[IJ(5, 1)] = v0.SigmaYPz();
    out.fC[IJ(5, 2)] = v0.SigmaZPz();
    out.fC[IJ(5, 3)] = v0.SigmaPxPz();
    out.fC[IJ(5, 4)] = v0.SigmaPyPz();
    out.fC[IJ(5, 5)] = v0.SigmaPz2();

    out.fC[IJ(6, 0)] = v0.SigmaXE();
    out.fC[IJ(6, 1)] = v0.SigmaYE();
    out.fC[IJ(6, 2)] = v0.SigmaZE();
    out.fC[IJ(6, 3)] = v0.SigmaPxE();
    out.fC[IJ(6, 4)] = v0.SigmaPyE();
    out.fC[IJ(6, 5)] = v0.SigmaPzE();
    out.fC[IJ(6, 6)] = v0.SigmaE2();

    out.fC[IJ(7, 7)] = Const::Initial_Css;

    out.fQ = 0;
    out.fNDF = 0;
    out.fChi2 = 0.;

    return out;
}

}  // namespace Tree2Secondaries::Legacy
