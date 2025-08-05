#ifndef ALICE_CONSTANTS_HXX
#define ALICE_CONSTANTS_HXX

namespace ALICE::Const {

static constexpr double AlmostZero{1.E-13};
static constexpr double Pi{3.1415927};
static constexpr double B2C{-0.299792458E-3};
static constexpr double Max_C0{10'000};           // SigmaY <= 100*100cm
static constexpr double Max_C2{10'000};           // SigmaZ <= 100*100cm
static constexpr double Max_C5{1.};               // SigmaSin <= 1*1
static constexpr double Max_C9{1.};               // SigmaTan <= 1*1
static constexpr double Max_C14{10'000};          // Sigma1/Pt <= 100*100 1/GeV
static constexpr double ResMisAlignPrec{2.5E-7};  // (= 0.0005 * 0.0005) a kind of a residual misalignment precision

static constexpr double CustomV0Finder_DCAmax{10.};  // (d=1.) max DCA between the daughter tracks
static constexpr double CustomV0Finder_Rmin{1.};     // (d=0.9) min radius of the V0 radius
static constexpr double CustomV0Finder_Rmax{200.};   // (d=100) max radius of the V0 radius

}  // namespace ALICE::Const

#endif  // ALICE_CONSTANTS_HXX
