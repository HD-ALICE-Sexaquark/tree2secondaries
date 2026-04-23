// Compare with: `original/AddDaughter.txt`

#include "Legacy/LegacyFitter.hxx"
#include "Legacy/LegacyParticle.hxx"

namespace Tree2Secondaries::Legacy {
// clang-format off

// Adds daughter to the current particle. Depending on the selected construction method uses: \n
// 1) Either simplifyed fast mathematics which consideres momentum and energy as
// independent variables and thus ignores constraint on the fixed mass (fConstructMethod = 0).
// In this case the mass of the daughter particle can be corrupted when the constructed vertex
// is added as the measurement and the mass of the output short-lived particle can become
// unphysical - smaller then the threshold. Implemented in the
// AddDaughterWithEnergyFit() function \n
// 2) Or slower but correct mathematics which requires that the masses of daughter particles
// stays fixed in the construction process (fConstructMethod = 2). Implemented in the
// AddDaughterWithEnergyFitMC() function.
// \param[in] Daughter - the daughter particle
void AddDaughter( Particle &part, const Particle& Daughter, double bz )
{

  if( part.fNDF<-1 ){ // first daughter -> just copy
    part.fNDF   = -1;
    part.fQ     =  Daughter.fQ;
    for( int i=0; i<7; i++) part.fP[i] = Daughter.fP[i];
    for( int i=0; i<28; i++) part.fC[i] = Daughter.fC[i];
    // fSFromDecay = 0;
    // fMassHypo = Daughter.fMassHypo;
    // SumDaughterMass = Daughter.SumDaughterMass;
    return;
  }

//   if(static_cast<int>(fConstructMethod) == 0)
    AddDaughterWithEnergyFit(part, Daughter, bz);
//   else if(static_cast<int>(fConstructMethod) == 2)
    // AddDaughterWithEnergyFitMC(Daughter);

//   SumDaughterMass += Daughter.SumDaughterMass;
//   fMassHypo = -1;
}

// clang-format on
}  // namespace Tree2Secondaries::Legacy
