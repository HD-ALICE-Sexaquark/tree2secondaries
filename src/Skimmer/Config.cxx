#include <algorithm>
#include <filesystem>
#include <fstream>
#include <utility>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "Utils/Logger.hxx"

#include "Skimmer/Config.hxx"

namespace {

template <typename T>
T Require(const json& node, std::string_view context, std::string_view key) {
    const auto it = node.find(key);
    if (it == node.end()) throw std::runtime_error(std::format("{}: missing required key \"{}\"", context, key));
    try {
        return it->get<T>();
    } catch (const json::exception& exc) {
        throw std::runtime_error(std::format("{}: key \"{}\" has the wrong type ({})", context, key, exc.what()));
    }
}

template <typename T>
T Optional(const json& node, std::string_view key, T fallback) {
    const auto it = node.find(key);
    if (it == node.end()) return fallback;
    try {
        return it->get<T>();
    } catch (const json::exception& exc) {
        throw std::runtime_error(std::format("optional key \"{}\" has the wrong type ({})", key, exc.what()));
    }
}

const json& RequireNode(const json& node, std::string_view context, std::string_view key) {
    const auto it = node.find(key);
    if (it == node.end()) {
        throw std::runtime_error(std::format("{}: missing required key \"{}\"", context, key));
    }
    return *it;
}

const json* FindNode(const json& node, std::string_view key) {
    const auto it = node.find(key);
    return it == node.end() ? nullptr : &*it;
}

const std::string& RequireString(const json& node, std::string_view context, std::string_view key) {
    const json& value = RequireNode(node, context, key);
    if (!value.is_string()) {
        throw std::runtime_error(std::format("{}: key \"{}\" must be a string", context, key));
    }
    return value.get_ref<const std::string&>();
}

std::string_view OptionalString(const json& node, std::string_view key, std::string_view fallback) {
    const json* value = FindNode(node, key);
    if (!value) return fallback;
    if (!value->is_string()) {
        throw std::runtime_error(std::format("optional key \"{}\" must be a string", key));
    }
    return value->get_ref<const std::string&>();
}

// A cut value is a scalar for `lower`/`upper` and a two-element array for `window`.
// Enforced here so that a window with a single number fails loudly instead of defaulting its upper edge to zero.
void ReadCutValue(const json& node, std::string_view context, std::string_view key, Skimmer::EDirection direction, double& low, double& high) {
    const json& value = RequireNode(node, context, key);
    if (direction == Skimmer::EDirection::kWindow) {
        if (!value.is_array() || value.size() != 2)
            throw std::runtime_error(std::format("{}: \"{}\" must be a two-element array [low, high] for a window cut", context, key));
        low = value[0].get<double>();
        high = value[1].get<double>();
        if (high < low) throw std::runtime_error(std::format("{}: \"{}\" has high < low", context, key));
    } else {
        if (!value.is_number()) throw std::runtime_error(std::format("{}: \"{}\" must be a number for a {} cut", context, key, AsString(direction)));
        low = value.get<double>();
        high = low;
    }
}

}  // namespace

namespace Skimmer {

// Every variable the cache must carry: declared and baseline, deduplicated and in that order.
std::vector<std::string> Config::CachedVariables() const {

    std::vector<std::string> names;
    const auto push = [&names](const std::string& name) {
        if (std::ranges::find(names, name) == names.end()) names.push_back(name);
    };

    for (const auto& variable : Variables) push(variable);
    for (const auto& cut : Baseline) push(cut.Variable);

    return names;
}

void Config::Print() const {
    Logger::Info(__FUNCTION__, "{:<12} = {}", "Path", Path);
    Logger::Info(__FUNCTION__, "{:<12} = {}", "Channel", AsString(Channel));
    Logger::Info(__FUNCTION__, "{:<12} = {:.3f}", "SignalMass", SignalMass);
    Logger::Info(__FUNCTION__, "{:<12} = {}", "Output", Output);
    Logger::Info(__FUNCTION__, "{:<12} = {}", "NEventsLimit", NEventsLimit == 0 ? std::string{"all"} : std::format("{}", NEventsLimit));
    Logger::Info(__FUNCTION__, "{:<12} = {}", "KeepRef", KeepReference);
    Logger::Info(__FUNCTION__, "{:<12} = {}", "KeepInjBkg", KeepInjectedBkg);
    if (Channel != EChannel::kLambdaPair) {
        Logger::Info(__FUNCTION__, "{:<12} = {}", "WeightsPt", WeightsPt.empty() ? "(none)" : WeightsPt);
        Logger::Info(__FUNCTION__, "{:<12} = {}", "WeightsRadius", WeightsRadius.empty() ? "(none)" : WeightsRadius);
    }
    Logger::Info(__FUNCTION__, "{:<12} = {} entries", "Samples", Samples.size());
    for (const auto& sample : Samples) {
        Logger::Info(__FUNCTION__, "{:<12}   {} (role: {})", "", sample.Path, AsString(sample.Role));
    }
    Logger::Info(__FUNCTION__, "{:<12} = {} cuts", "Baseline", Baseline.size());
    for (const auto& cut : Baseline) {
        if (cut.Direction == EDirection::kWindow) {
            Logger::Info(__FUNCTION__, "{:<12}   {} window [{}, {}]", "", cut.Variable, cut.Value, cut.ValueUpper);
        } else {
            Logger::Info(__FUNCTION__, "{:<12}   {} {} {}", "", cut.Variable, AsString(cut.Direction), cut.Value);
        }
    }
    const auto cached = CachedVariables();
    Logger::Info(__FUNCTION__, "{:<12} = {} cached", "Variables", cached.size());
    for (const auto& var : cached) {
        Logger::Info(__FUNCTION__, "{:<12}   {}", "", var);
    }
}

Config Load(std::string_view path) {

    std::ifstream stream{std::string(path)};
    if (!stream.is_open()) throw std::runtime_error(std::format("could not open config \"{}\"", path));

    json json_file;
    try {
        stream >> json_file;
    } catch (const json::exception& exc) {
        throw std::runtime_error(std::format("could not parse config \"{}\": {}", path, exc.what()));
    }

    Config config;
    config.Path = std::string(path);
    config.Channel = AsChannel(RequireString(json_file, "config", "channel"));
    config.SignalMass = Require<double>(json_file, "config", "signal_mass");

    // -- where the cache goes

    // A full path rather than a directory: this key is how the consumer finds the cache, and a directory would
    // leave it re-deriving the file name from the channel.
    config.Output = RequireString(json_file, "config", "output");
    if (std::filesystem::path{config.Output}.extension() != ".root") {
        throw std::runtime_error(std::format("config: \"output\" must be a path ending in .root, not \"{}\"", config.Output));
    }

    // -- how many events to read

    // Read as signed: a negative literal would otherwise wrap into a huge limit, i.e. into "all of them".
    const auto n_events_limit = Optional<std::int64_t>(json_file, "n_events_limit", 0);
    if (n_events_limit < 0) throw std::runtime_error(std::format("config: \"n_events_limit\" cannot be negative ({})", n_events_limit));
    config.NEventsLimit = static_cast<std::uint64_t>(n_events_limit);

    // -- what MC truth is allowed in

    config.KeepReference = Optional<bool>(json_file, "keep_reference", false);
    config.KeepInjectedBkg = Optional<bool>(json_file, "keep_injected_bkg", false);

    // -- shape weights

    config.WeightsPt = OptionalString(json_file, "weights_pt", "");
    config.WeightsRadius = OptionalString(json_file, "weights_radius", "");
    if (config.Channel == EChannel::kLambdaPair && (!config.WeightsPt.empty() || !config.WeightsRadius.empty())) {
        throw std::runtime_error(
            "config: \"weights_pt\"/\"weights_radius\" are built to reshape the flat injected antisexaquarks and do not apply to h-dibaryon MC.");
    }

    // -- samples

    // The RNTuple to read is not named here: the channel fixes the schema, and the schema fixes the tuple name
    // (see `Skimmer/Channels.hxx`). A config could otherwise name a tuple its channel cannot read.
    const json& samples = RequireNode(json_file, "config", "samples");
    if (!samples.is_array() || samples.empty()) throw std::runtime_error("config: \"samples\" must be a non-empty array");
    if (samples.size() > kMaxSamples) {
        // `SampleIndex` is one byte per cache row; see `Skimmer/Writers.hxx`.
        throw std::runtime_error(
            std::format("config: \"samples\" holds {} entries, more than the {} a cache row can address", samples.size(), kMaxSamples));
    }
    for (const auto& node : samples) {
        Sample sample;
        sample.Path = RequireString(node, "samples[]", "path");
        sample.Role = AsRole(OptionalString(node, "role", "both"));
        config.Samples.push_back(std::move(sample));
    }

    // -- samples' roles

    // Either every sample is "both" -- one production serving as its own background -- or every sample names a side.
    // Half-declared roles leave it ambiguous which files the background normalization is built from, and a normalization built from the wrong file
    // set is wrong by orders of magnitude without looking wrong anywhere.
    const bool any_both = std::ranges::any_of(config.Samples, [](const Sample& s) { return s.Role == ERole::kBoth; });
    const bool all_both = std::ranges::all_of(config.Samples, [](const Sample& s) { return s.Role == ERole::kBoth; });
    if (any_both && !all_both) {
        throw std::runtime_error(
            "samples[]: a config either declares every sample \"both\", or gives every sample an explicit \"signal\"/\"background\" role. "
            "Mixing the two leaves it ambiguous which files the background normalization is built from.");
    }
    if (!all_both) {
        for (const auto role : {ERole::kSignal, ERole::kBackground}) {
            if (!std::ranges::any_of(config.Samples, [role](const Sample& s) { return s.Role == role; })) {
                throw std::runtime_error(
                    std::format("samples[]: roles are declared but no sample has role \"{}\". Either add one, or declare every sample \"both\".",
                                AsString(role)));
            }
        }
    }

    // -- cached variables

    const json& variables = RequireNode(json_file, "config", "variables");
    if (!variables.is_array() || variables.empty()) {
        throw std::runtime_error("config: \"variables\" must be a non-empty array");
    }
    for (const auto& node : variables) {
        if (!node.is_string()) {
            throw std::runtime_error("variables[]: every entry must be a variable name, as a string");
        }
        const auto& name = node.get_ref<const std::string&>();
        if (std::ranges::find(config.Variables, name) != config.Variables.end()) {
            throw std::runtime_error(std::format("variables[]: \"{}\" is listed twice", name));
        }
        config.Variables.push_back(name);
    }

    // -- baseline preselection

    if (const json* baseline = FindNode(json_file, "baseline")) {
        if (!baseline->is_array()) {
            throw std::runtime_error("config: \"baseline\" must be an array");
        }
        for (const auto& node : *baseline) {
            BaselineCut cut;
            cut.Variable = RequireString(node, "baseline[]", "variable");
            cut.Direction = AsDirection(RequireString(node, "baseline[]", "direction"));
            ReadCutValue(node, "baseline[]", "value", cut.Direction, cut.Value, cut.ValueUpper);
            config.Baseline.push_back(std::move(cut));
        }
    }

    return config;
}

}  // namespace Skimmer
