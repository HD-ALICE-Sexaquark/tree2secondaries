#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string_view>

namespace Skimmer::Process {

// # Process Labels # //

// Backgrounds are per-process and MC-truth-based. The four labels below are deliberately uniform across
// all four analysis channels:
//
//   kSignal         `IsTrueSignal`: every leg is a true daughter of one same injection, under the
//                   hypothesis the candidate was built with
//   kHybrid         neither of the two above -- a mix of signal and non-signal legs, or a true leg of
//                   the wrong species. The most dangerous background, because it peaks where the
//                   signal peaks
//   kCombinatorial  `IsRealBkg`: carries no signal at all, i.e. neither the candidate nor any of its
//                   constituents descends from an injection
//   kReference      the deliberately-wrong-sign combination `t2ds` also reconstructs:
//                   `POD::Sexaquark::IsWrongSignChannel` for the sexaquark, the mixed Lambda-AntiLambda
//                   pair (`Cached::Hdibaryon::IsMixedChannel`) for the h-dibaryon. It is a background
//                   *model*, not a background: keep it out of the sum unless you mean to use it as one.
//
// `kUnknown` only ever appears in real data, where there is no MC truth to label with.

enum EProcess : std::uint8_t {
    kSignal = 0,
    kHybrid = 1,
    kCombinatorial = 2,
    kReference = 3,
    kUnknown = 4,
};

inline constexpr auto Name_Process = std::to_array<std::string_view>({"Signal", "Hybrid", "Combinatorial", "Reference", "Unknown"});
inline constexpr std::size_t NProcesses = Name_Process.size();

[[nodiscard]] constexpr std::string_view Name(EProcess process) { return Name_Process[static_cast<std::size_t>(process)]; }

// The three MC-truth sets are mutually exclusive and exhaustive: `TrueSignal` and `RealBkg` are defined
// positively and `Hybrid` is whatever is left. See `common/docs/MC_LABELS.md`.
//
// Takes the two stored flags rather than a POD, because the sexaquark channels carry them on a
// `POD::Linked::InjectedSexa` and the h-dibaryon on a `POD::Extended::McParticle` -- unrelated types
// that only agree on these two fields.
[[nodiscard]] constexpr EProcess Classify(bool is_true_signal, bool is_real_bkg) {
    if (is_true_signal) return kSignal;
    if (is_real_bkg) return kCombinatorial;
    return kHybrid;
}

}  // namespace Skimmer::Process
