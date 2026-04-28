#pragma once

#include <array>
#include <string_view>

#if defined(__AVX512F__)
#define T2S_SIMD_ALIGN 64
#elif defined(__AVX2__)
#define T2S_SIMD_ALIGN 32
#else
#define T2S_SIMD_ALIGN 16
#endif

namespace Tree2Secondaries {

enum PID_StableParticle { PiMinus, PiPlus, NegKaon, PosKaon, AntiProton, Proton, AntiNeutron, Neutron, NStableParticles };
enum PID_V0 { KaonZeroShort, AntiLambda, Lambda, NV0s };
enum EReactionChannel { A, D, E, H, NReactionChannels };

namespace Const {
inline constexpr std::string_view TreeName_Events = "Events";
inline constexpr std::string_view TreeName_PackedEvents = "PackedEvents";
inline constexpr std::string_view TreeName_Injected = "Injected";

inline constexpr std::array<std::string_view, NStableParticles> Particle_Acronym{"PM", "PP", "NK", "PK", "AP", "P", "AN", "N"};
inline constexpr std::array<int, NStableParticles> Particle_PdgCode{-211, 211, -321, 321, -2212, 2212, -2112, 2112};
inline constexpr std::array<double, NStableParticles> Particle_Mass{0.13957040, 0.13957040, 0.49367700, 0.49367700,
                                                                    0.93827210, 0.93827210, 0.93956540, 0.93956540};
inline constexpr std::array<int, NStableParticles> Particle_Charge{-1, +1, -1, +1, -1, +1, 0, 0};

inline constexpr std::array<std::string_view, NV0s> V0_Acronym{"K0S", "AL", "L"};
inline constexpr std::array<int, NV0s> V0_PdgCode{310, -3122, 3122};
inline constexpr std::array<double, NV0s> V0_Mass{0.49761100, 1.1156830, 1.1156830};
inline constexpr std::array<PID_StableParticle, NV0s> V0_NegativePID{PID_StableParticle::PiMinus, PID_StableParticle::AntiProton,
                                                                     PID_StableParticle::PiMinus};
inline constexpr std::array<PID_StableParticle, NV0s> V0_PositivePID{PID_StableParticle::PiPlus, PID_StableParticle::PiPlus,
                                                                     PID_StableParticle::Proton};

inline constexpr std::array<char, NReactionChannels> ReactionChannel_Char{'A', 'D', 'E', 'H'};
inline constexpr std::array<PID_StableParticle, NReactionChannels> ReactionChannel_NucleonPID{PID_StableParticle::Neutron, PID_StableParticle::Proton,
                                                                                              PID_StableParticle::Proton, PID_StableParticle::Proton};

inline constexpr double Kappa = 0.000299792458;  // c in natural units, (GeV/c) / (kG/cm)
inline constexpr double AbsAlmostZero = 1.E-8;
inline constexpr double BigNumber = 1.E8;
inline constexpr int DummyInt = -1;
inline constexpr float DummyFloat = -999.;
inline constexpr double DummyDouble = -999.;
inline constexpr double StandardSexaquarkMass = 1.8;  // (GeV/c^2)
inline constexpr int ReactionID_Offset = 600;
inline constexpr double Initial_Css = 1.;
inline constexpr int NInjectedPerEvent = 20;  // number of injected sexaquark-material reactions per event in dedicated MC
inline constexpr int SignalGeneratorID = 2;

constexpr EReactionChannel ReactionChannelFromChar(char c) {
    for (size_t i = 0; i < NReactionChannels; ++i) {
        if (Const::ReactionChannel_Char[i] == c) return static_cast<EReactionChannel>(i);
    }
    return EReactionChannel::A;  // fallback
}

}  // namespace Const

}  // namespace Tree2Secondaries
