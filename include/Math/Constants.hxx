#ifndef T2S_CONSTANTS_HXX
#define T2S_CONSTANTS_HXX

#include <string>
#include <unordered_map>
#include <vector>

#if defined(__AVX512F__)
#define T2S_SIMD_ALIGN 64
#elif defined(__AVX2__)
#define T2S_SIMD_ALIGN 32
#else
#define T2S_SIMD_ALIGN alignof(double)
#endif

namespace Tree2Secondaries {

namespace Name {
static constexpr std::string Events{"Events"};
static constexpr std::string PackedEvents{"PackedEvents"};
static constexpr std::string Injected{"Injected"};
}  // namespace Name

enum EParticle { PiMinus, PiPlus, NegKaon, PosKaon, KaonZeroShort, AntiProton, Proton, AntiNeutron, Neutron, AntiLambda, Lambda };

namespace Particle {
static const std::vector<std::string> Acronym{"PM", "PP", "NK", "PK", "K0S", "AP", "P", "AN", "N", "AL", "L"};
static const std::vector<int> PdgCode{-211, 211, -321, 321, 310, -2212, 2212, -2112, 2112, -3122, 3122};

// in (GeV/c^2)
static const std::vector<double> Mass{0.13957040, 0.13957040, 0.49367700, 0.49367700, 0.49761100, 0.93827210,
                                      0.93827210, 0.93956540, 0.93956540, 1.1156830,  1.1156830};
static const std::vector<int> Charge{-1, +1, -1, +1, 0, -1, +1, 0, 0, 0, 0};
}  // namespace Particle

enum EReactionChannel : char { A = 'A', D = 'D', E = 'E', H = 'H' };
static std::unordered_map<EReactionChannel, EParticle> ReactionNucleonPID{{EReactionChannel::A, EParticle::Neutron},
                                                                          {EReactionChannel::D, EParticle::Proton},
                                                                          {EReactionChannel::E, EParticle::Proton},
                                                                          {EReactionChannel::H, EParticle::Proton}};

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
