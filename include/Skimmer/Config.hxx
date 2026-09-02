#pragma once

#include <cstddef>
#include <cstdint>
#include <format>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace Skimmer {

// ## The Skim Config File ## //

// Read a JSON file.
// It holds only what changes the skim's output: which files to read, which variables to cache, and which
// preselection to apply. Every variable named here must exist in the compiled registry of
// `Skimmer/VariableRegistry.hxx` for the declared channel.
//
// Everything downstream needs but the skim does not act on -- the observable's binning, the scan ranges,
// the figure of merit, the guards -- lives in a separate analysis config owned by the consumer, so that a
// re-tune of the optimizer does not imply a re-skim. See `docs/ANALYSIS_CONFIG.md`. The output's
// `CacheSource` record carries a digest of this file -- channel, signal mass, and the path it was loaded
// from -- so a cache says which config produced it, but not what that config said: the baseline cuts live
// only here.

enum class EChannel : std::uint8_t {
    kChannelA,    // antisexaquark + neutron -> antilambda + K0S
    kChannelD,    // antisexaquark + proton -> antilambda + K+
    kChannelH,    // antisexaquark + proton -> K+ + K+ + X
    kLambdaPair,  // (anti)h-dibaryon -> (anti)lambda + (anti)lambda
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

// Every cache row names the sample it came from with a single byte (`SampleIndex`, see `Skimmer/Writers.hxx`),
// which is what bounds how many samples one config may declare.
inline constexpr std::size_t kMaxSamples = 256;

struct Sample {
    std::string Path;
    ERole Role{ERole::kBoth};
};

struct BaselineCut {
    std::string Variable;
    EDirection Direction{EDirection::kLower};
    double Value{0.};
    double ValueUpper{0.};  // only read when `Direction == kWindow`
};

[[nodiscard]] inline std::string_view AsString(EChannel channel) {
    switch (channel) {
        case EChannel::kChannelA:
            return "ChannelA";
        case EChannel::kChannelD:
            return "ChannelD";
        case EChannel::kChannelH:
            return "ChannelH";
        case EChannel::kLambdaPair:
            return "LambdaPair";
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
    if (name == "LambdaPair") return EChannel::kLambdaPair;
    throw std::runtime_error(std::format("unknown channel \"{}\" (expected ChannelA, ChannelD, ChannelH or LambdaPair)", name));
}

[[nodiscard]] inline EDirection AsDirection(std::string_view name) {
    if (name == "lower") return EDirection::kLower;
    if (name == "upper") return EDirection::kUpper;
    if (name == "window") return EDirection::kWindow;
    throw std::runtime_error(std::format("unknown direction \"{}\" (expected lower, upper or window)", name));
}

struct Config {
    EChannel Channel{EChannel::kChannelA};
    double SignalMass{0.};
    std::vector<Sample> Samples;
    std::vector<std::string> Variables;
    std::vector<BaselineCut> Baseline;
    std::string WeightsPt;
    std::string WeightsRadius;

    std::string Path;  // where this config was loaded from

    [[nodiscard]] std::vector<std::string> CachedVariables() const;
    void Print() const;
};

Config Load(std::string_view path);

}  // namespace Skimmer
