#ifndef T2S_CONSTANTS_HXX
#define T2S_CONSTANTS_HXX

#include <string_view>
#include <vector>

namespace Tree2Secondaries {

enum EReactionChannel { All, A, D, E, H, AntiA, AntiD, AntiE, AntiH };

namespace Name {
static const std::vector<std::string_view> ReactionChannel{"AllChannels", "A", "D", "E", "H", "AntiA", "AntiD", "AntiE", "AntiH"};
static constexpr std::string_view Events{"Events"};
static constexpr std::string_view PackedEvents{"PackedEvents"};
}  // namespace Name

enum EParticle { PiMinus, PiPlus, NegKaon, PosKaon, KaonZeroShort, AntiProton, Proton, AntiNeutron, Neutron, AntiLambda, Lambda };

namespace Particle {
static const std::vector<std::string_view> Acronym{"PM", "PP", "NK", "PK", "K0S", "AP", "P", "AN", "N", "AL", "L"};
static const std::vector<int> PdgCode{-211, 211, -321, 321, 310, -2212, 2212, -2112, 2112, -3122, 3122};

// in (GeV/c^2)
static const std::vector<double> Mass{0.13957040, 0.13957040, 0.49367700, 0.49367700, 0.49761100, 0.93827210,
                                      0.93827210, 0.93956540, 0.93956540, 1.1156830,  1.1156830};
}  // namespace Particle

namespace Const {
static constexpr double Kappa{0.000299792458};  // (GeV/c) / (kG/cm)
static constexpr double AbsAlmostZero{1.E-8};
static constexpr double BigNumber{1.E8};
static constexpr int DummyInt{-1};
static constexpr float DummyFloat{-999.};
static constexpr double DummyDouble{-999.};
static constexpr double StandardSexaquarkMass{1.8};  // (GeV/c^2)
}  // namespace Const

}  // namespace Tree2Secondaries

#endif  // T2S_CONSTANTS_HXX
