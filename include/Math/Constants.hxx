#pragma once

#include <cstddef>
#include <string_view>

#if defined(__AVX512F__)
#define T2S_SIMD_ALIGN 64
#elif defined(__AVX2__)
#define T2S_SIMD_ALIGN 32
#else
#define T2S_SIMD_ALIGN 16
#endif

namespace Tree2Secondaries::Const {

inline constexpr std::string_view TreeName_Events = "Events";
inline constexpr std::string_view TreeName_PackedEvents = "PackedEvents";
inline constexpr std::string_view TreeName_Injected = "Injected";

inline constexpr double Kappa = 0.000299792458;  // c in natural units, (GeV/c) / (kG/cm)
inline constexpr double AbsAlmostZero = 1.E-8;
inline constexpr double BigNumber = 1.E8;
inline constexpr int DummyInt = -1;
inline constexpr float DummyFloat = -999.;
inline constexpr double DummyDouble = -999.;
inline constexpr double Initial_Css = 1.;

inline constexpr double StandardSexaquarkMass = 1.8;  // (GeV/c^2)
inline constexpr int ReactionID_Offset = 600;
inline constexpr std::size_t NInjectedPerEvent = 20;  // number of injected sexaquark-material reactions per event in dedicated MC
enum EGenerator : char { kHijing, kInjectedAntiNeutron, kInjectedAntiSexaquarkReaction };

}  // namespace Tree2Secondaries::Const
