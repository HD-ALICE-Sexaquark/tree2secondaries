#pragma once

#include <cstdint>
#include <format>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace Skimmer {

// ## The Config File ## //

// Read a JSON file.
// Every variable referenced here must exist in the compiled registry of `Skimmer/VariableRegistry.hxx` for the declared channel.
// The idea is that future macros/dashboards/optimizers mirror this json structure.

enum class EChannel : std::uint8_t {
    kChannelA,   // antisexaquark + neutron -> antilambda + K0S
    kChannelD,   // antisexaquark + proton -> antilambda + K+
    kChannelH,   // antisexaquark + proton -> K+ + K+ + X
    kHdibaryon,  // (anti)h-dibaryon -> (anti)lambda + (anti)lambda
};

enum class EDirection : std::uint8_t {
    kLower,   // keep x >= value
    kUpper,   // keep x <= value
    kWindow,  // keep value <= x <= value_upper
};

enum class ERole : std::uint8_t {
    kSignal,      // supplies S
    kBackground,  // supplies B
    kBoth,        // one production serving as its own background -- the default
};

struct Sample {
    std::string Path;
    std::string NTuple;
    unsigned int RunNumber{0};
    ERole Role{ERole::kBoth};
    // Injected signals per event in *this* production. 0 means "use the channel trait", the same
    // convention `RunNumber` uses. It exists because `Traits::kNInjectedPerEvent` is a compile-time
    // constant of the channel, not of the production: a sexaquark production skimmed under `Hdibaryon`
    // traits would otherwise claim 100 injected h-dibaryons per event while containing none.
    unsigned int NInjectedPerEvent{0};
};

struct BaselineCut {
    std::string Variable;
    EDirection Direction{EDirection::kLower};
    double Value{0.};
    double ValueUpper{0.};  // only read when `Direction == kWindow`
};

struct CutVariable {
    std::string Name;
    EDirection Direction{EDirection::kLower};
    double RangeMin{0.};
    double RangeMax{0.};
    unsigned int Steps{100};
    double Initial{0.};
    double InitialUpper{0.};  // only read when `Direction == kWindow`
    bool InGrid{false};
};

struct Observable {
    std::string Variable;
    unsigned int Bins{200};
    double Min{0.};
    double Max{5.};
};

// Not used at all by `Skimmer`, but validated here.
struct FOM {
    std::string Formula{"punzi"};    // punzi | asimov | poisson
    double A{3.};                    // target significance, in sigmas (Punzi only)
    double FSyst{0.2};               // assumed relative systematic on the background
    bool HasNSignalExpected{false};  //
    double NSignalExpected{0.};      // expected signal candidates in the full real-data sample at 100% efficiency, `asimov` and `poisson` need it
    bool HasSignalYield{false};
    double DnDy{0.};
    double DeltaY{0.};
    double BranchingRatio{1.};
    unsigned int NInjectedSpecies{1};
    double InteractionProbability{1.};
    std::string YieldSource;
};

struct Guards {
    unsigned int MinRawSignal{20};
    unsigned int MinRawBackground{20};
};

[[nodiscard]] inline std::string_view AsString(EChannel channel) {
    switch (channel) {
        case EChannel::kChannelA:
            return "ChannelA";
        case EChannel::kChannelD:
            return "ChannelD";
        case EChannel::kChannelH:
            return "ChannelH";
        case EChannel::kHdibaryon:
            return "Hdibaryon";
    }
    return "Unknown";
}

[[nodiscard]] inline std::string_view AsString(EDirection direction) {
    switch (direction) {
        case EDirection::kLower:
            return "lower";
        case EDirection::kUpper:
            return "upper";
        case EDirection::kWindow:
            return "window";
    }
    return "unknown";
}

[[nodiscard]] inline std::string_view AsString(ERole role) {
    switch (role) {
        case ERole::kSignal:
            return "signal";
        case ERole::kBackground:
            return "background";
        case ERole::kBoth:
            return "both";
    }
    return "unknown";
}

[[nodiscard]] inline ERole AsRole(std::string_view name) {
    if (name == "signal") return ERole::kSignal;
    if (name == "background") return ERole::kBackground;
    if (name == "both") return ERole::kBoth;
    throw std::runtime_error(std::format("unknown role \"{}\" (expected signal, background or both)", name));
}

[[nodiscard]] inline EChannel AsChannel(std::string_view name) {
    if (name == "ChannelA") return EChannel::kChannelA;
    if (name == "ChannelD") return EChannel::kChannelD;
    if (name == "ChannelH") return EChannel::kChannelH;
    if (name == "Hdibaryon") return EChannel::kHdibaryon;
    throw std::runtime_error(std::format("unknown channel \"{}\" (expected ChannelA, ChannelD, ChannelH or Hdibaryon)", name));
}

[[nodiscard]] inline EDirection AsDirection(std::string_view name) {
    if (name == "lower") return EDirection::kLower;
    if (name == "upper") return EDirection::kUpper;
    if (name == "window") return EDirection::kWindow;
    throw std::runtime_error(std::format("unknown direction \"{}\" (expected lower, upper, window, symm_inside or symm_outside)", name));
}

struct Config {
    EChannel Channel{EChannel::kChannelA};
    std::string Hypothesis;
    double SignalMass{0.};
    std::vector<Sample> Samples;
    Observable Obs;
    std::vector<BaselineCut> Baseline;
    std::vector<CutVariable> Variables;
    FOM FigureOfMerit;
    Guards Guard;
    std::vector<std::string> SentinelOk;  // not used by Skimmer, but validated

    std::string Path;  // where this config was loaded from

    [[nodiscard]] std::vector<std::string> CachedVariables() const;
    void Print() const;
};

Config Load(std::string_view path);

}  // namespace Skimmer
