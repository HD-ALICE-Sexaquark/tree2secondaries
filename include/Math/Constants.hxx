#ifndef T2S_CONSTANTS_HXX
#define T2S_CONSTANTS_HXX

namespace Tree2Secondaries {

enum class Channel { A, D, E, H, AntiA, AntiD, AntiE, AntiH };

namespace PdgCode {
constexpr int AntiLambda{-3122};
constexpr int Lambda{3122};
constexpr int KaonZeroShort{310};
constexpr int AntiNeutron{-2112};
constexpr int Neutron{2112};
constexpr int AntiProton{-2212};
constexpr int Proton{2212};
constexpr int NegKaon{-321};
constexpr int PosKaon{321};
constexpr int PiMinus{-211};
constexpr int PiPlus{211};
}  // namespace PdgCode

namespace Const {
constexpr double MassNeutron{0.93956540};  // (GeV/c^2)
constexpr double MassProton{0.93827210};   // (GeV/c^2)
constexpr double MassKaon{0.49367700};     // (GeV/c^2)
constexpr double MassPion{0.13957040};     // (GeV/c^2)
constexpr double Kappa{0.000299792458};    // (GeV/c) / (kG/cm)
constexpr double LocalSmall{1.E-6};
constexpr double AlmostZero{1.E-8};
}  // namespace Const

}  // namespace Tree2Secondaries

#endif  // T2S_CONSTANTS_HXX
