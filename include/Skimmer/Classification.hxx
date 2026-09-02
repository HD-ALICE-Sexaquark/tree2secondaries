#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string_view>

namespace Skimmer::Classification {

// # Classification Labels # //

// The labels are MC-truth-based and deliberately uniform across all four analysis channels:
//
//   kSignal            `IsTrueSignal`: every daughter is a true daughter of one same injection, under the
//                      hypothesis the candidate was built with
//   kRealBkg           `IsRealBkg`: carries no signal at all, i.e. neither the candidate nor any of its
//                      constituents descends from an injection
//   kHybrid            neither of the two above -- a mix of signal and non-signal daughters, or a true daughter of
//                      the wrong species.
//   kReferenceReal     background models:
//                      - wrong-sign channels for the antisexaquark,
//                      - mixed lambda-antilambda pairs for the (anti)h-dibaryon
//                      where none of the (grand)daughters carry signal
//   kReferenceHybrid   same as above, where at least one (grand)daughter carry signal
//   kUnknown           only ever appears in real data (PENDING), where there is no MC truth to label with

enum EClassification : std::uint8_t {
    kSignal = 0,
    kRealBkg = 1,
    kHybrid = 2,
    kReferenceReal = 3,
    kReferenceHybrid = 4,
    kUnknown = 5,
};

inline constexpr auto Name_Classification =
    std::to_array<std::string_view>({"Signal", "RealBkg", "Hybrid", "ReferenceReal", "ReferenceHybrid", "Unknown"});
inline constexpr std::size_t NClassifications = Name_Classification.size();

[[nodiscard]] constexpr std::string_view Name(EClassification classification) {
    return Name_Classification[static_cast<std::size_t>(classification)];
}

// The three MC-truth sets are mutually exclusive and exhaustive: `TrueSignal` and `RealBkg` are defined positively and `Hybrid` is whatever is left.
// See `common/docs/MC_LABELS.md`.
[[nodiscard]] constexpr EClassification Classify(bool is_true_signal, bool is_real_bkg) {
    if (is_true_signal) return kSignal;
    if (is_real_bkg) return kRealBkg;
    return kHybrid;
}

}  // namespace Skimmer::Classification
