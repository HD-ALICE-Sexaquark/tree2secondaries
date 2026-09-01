#include <algorithm>
#include <fstream>
#include <utility>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "Utils/Logger.hxx"

#include "Skimmer/Config.hxx"

namespace {

// Reads a required key, and says which key was missing rather than letting nlohmann throw a type error three frames deeper.
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
    return it->get<T>();
}

// A cut value is a scalar for `lower`/`upper` and a two-element array for `window`.
// Enforced here so that a window with a single number fails loudly instead of defaulting its upper edge to zero.
inline void ReadCutValue(const json& node, std::string_view context, std::string_view key, Skimmer::EDirection direction, double& low, double& high) {
    const auto it = node.find(key);
    if (it == node.end()) throw std::runtime_error(std::format("{}: missing required key \"{}\"", context, key));
    if (direction == Skimmer::EDirection::kWindow) {
        if (!it->is_array() || it->size() != 2)
            throw std::runtime_error(std::format("{}: \"{}\" must be a two-element array [low, high] for a window cut", context, key));
        low = (*it)[0].get<double>();
        high = (*it)[1].get<double>();
        if (high < low) throw std::runtime_error(std::format("{}: \"{}\" has high < low", context, key));
    } else {
        if (!it->is_number()) throw std::runtime_error(std::format("{}: \"{}\" must be a number for a {} cut", context, key, AsString(direction)));
        low = it->get<double>();
        high = low;
    }
}

}  // namespace

namespace Skimmer {

// Every variable the cache must carry: the observable, then the scanned variables, then the
// baseline ones, deduplicated and in that order.
// Baseline variables are stored even after applying their cuts.
std::vector<std::string> Config::CachedVariables() const {

    std::vector<std::string> names;
    const auto push = [&names](const std::string& name) {
        if (std::ranges::find(names, name) == names.end()) names.push_back(name);
    };

    push(Obs.Variable);
    for (const auto& variable : Variables) push(variable.Name);
    for (const auto& cut : Baseline) push(cut.Variable);

    return names;
}

void Config::Print() const {
    Logger::Info(__FUNCTION__, "{:<20} = {}", "Path", Path);
    Logger::Info(__FUNCTION__, "{:<20} = {}", "Channel", AsString(Channel));
    Logger::Info(__FUNCTION__, "{:<20} = {}", "Hypothesis", Hypothesis);
    Logger::Info(__FUNCTION__, "{:<20} = {:.2f}", "SignalMass", SignalMass);
    Logger::Info(__FUNCTION__, "{:<20} = {} entries", "Samples", Samples.size());
    for (const auto& sample : Samples)
        Logger::Info(__FUNCTION__, "{:<20}   {} ({}, role {}{})", "", sample.Path, sample.NTuple, AsString(sample.Role),
                     sample.NInjectedPerEvent > 0 ? std::format(", {}/event injected", sample.NInjectedPerEvent) : "");
    Logger::Info(__FUNCTION__, "{:<20} = {} [{}, {}] in {} bins", "Observable", Obs.Variable, Obs.Min, Obs.Max, Obs.Bins);
    Logger::Info(__FUNCTION__, "{:<20} = {} cuts", "Baseline", Baseline.size());
    for (const auto& cut : Baseline) {
        if (cut.Direction == EDirection::kWindow)
            Logger::Info(__FUNCTION__, "{:<20}   {} window [{}, {}]", "", cut.Variable, cut.Value, cut.ValueUpper);
        else
            Logger::Info(__FUNCTION__, "{:<20}   {} {} {}", "", cut.Variable, AsString(cut.Direction), cut.Value);
    }
    Logger::Info(__FUNCTION__, "{:<20} = {} variables", "Scan", Variables.size());
    for (const auto& variable : Variables)
        Logger::Info(__FUNCTION__, "{:<20}   {} ({}, {} steps over [{}, {}]{})", "", variable.Name, AsString(variable.Direction), variable.Steps,
                     variable.RangeMin, variable.RangeMax, variable.InGrid ? ", in grid" : "");
    Logger::Info(__FUNCTION__, "{:<20} = {} (a = {}, f_syst = {})", "Fom", FigureOfMerit.Formula, FigureOfMerit.A, FigureOfMerit.FSyst);
    if (FigureOfMerit.HasSignalYield)
        Logger::Info(__FUNCTION__, "{:<20}   dN/dy {} x dy {} x BR {} x {} species x P(int) {}", "", FigureOfMerit.DnDy, FigureOfMerit.DeltaY,
                     FigureOfMerit.BranchingRatio, FigureOfMerit.NInjectedSpecies, FigureOfMerit.InteractionProbability);
    Logger::Info(__FUNCTION__, "{:<20} = {} signal, {} background", "Guards", Guard.MinRawSignal, Guard.MinRawBackground);
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
    config.Channel = AsChannel(Require<std::string>(json_file, "config", "channel"));
    config.Hypothesis = Require<std::string>(json_file, "config", "hypothesis");
    config.SignalMass = Optional<double>(json_file, "signal_mass", 0.);

    // -- samples

    const auto& samples = Require<json>(json_file, "config", "samples");
    if (!samples.is_array() || samples.empty()) throw std::runtime_error("config: \"samples\" must be a non-empty array");
    for (const auto& node : samples) {
        Sample sample;
        sample.Path = Require<std::string>(node, "samples[]", "path");
        sample.NTuple = Require<std::string>(node, "samples[]", "ntuple");
        sample.RunNumber = Optional<unsigned int>(node, "run_number", 0U);
        sample.Role = AsRole(Optional<std::string>(node, "role", std::string("both")));
        sample.NInjectedPerEvent = Optional<unsigned int>(node, "n_injected_per_event", 0U);
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

    // -- observable

    const auto& observable = Require<json>(json_file, "config", "observable");
    config.Obs.Variable = Require<std::string>(observable, "observable", "variable");
    config.Obs.Bins = Optional<unsigned int>(observable, "bins", 200U);
    const auto& range = Require<json>(observable, "observable", "range");
    if (!range.is_array() || range.size() != 2) throw std::runtime_error("observable: \"range\" must be a two-element array");
    config.Obs.Min = range[0].get<double>();
    config.Obs.Max = range[1].get<double>();

    // -- baseline preselection

    const json baseline = Optional<json>(json_file, "baseline", json::array());
    for (const auto& node : baseline) {
        BaselineCut cut;
        cut.Variable = Require<std::string>(node, "baseline[]", "variable");
        cut.Direction = AsDirection(Require<std::string>(node, "baseline[]", "direction"));
        ReadCutValue(node, "baseline[]", "value", cut.Direction, cut.Value, cut.ValueUpper);
        config.Baseline.push_back(std::move(cut));
    }

    // -- scanned cut variables

    const auto& variables = Require<json>(json_file, "config", "variables");
    if (!variables.is_array() || variables.empty()) throw std::runtime_error("config: \"variables\" must be a non-empty array");
    for (const auto& node : variables) {
        CutVariable variable;
        variable.Name = Require<std::string>(node, "variables[]", "name");
        variable.Direction = AsDirection(Require<std::string>(node, "variables[]", "direction"));
        const auto& scan_range = Require<json>(node, "variables[]", "range");
        if (!scan_range.is_array() || scan_range.size() != 2)
            throw std::runtime_error(std::format("variables[{}]: \"range\" must be a two-element array", variable.Name));
        variable.RangeMin = scan_range[0].get<double>();
        variable.RangeMax = scan_range[1].get<double>();
        if (variable.RangeMax <= variable.RangeMin)
            throw std::runtime_error(std::format("variables[{}]: \"range\" is empty or inverted", variable.Name));
        variable.Steps = Optional<unsigned int>(node, "steps", 100U);
        if (variable.Steps < 2) throw std::runtime_error(std::format("variables[{}]: \"steps\" must be at least 2", variable.Name));
        ReadCutValue(node, std::format("variables[{}]", variable.Name), "initial", variable.Direction, variable.Initial, variable.InitialUpper);
        variable.InGrid = Optional<bool>(node, "in_grid", false);
        config.Variables.push_back(std::move(variable));
    }

    // -- figure of merit

    const auto& fom = Require<json>(json_file, "config", "fom");
    config.FigureOfMerit.Formula = Require<std::string>(fom, "fom", "formula");
    if (config.FigureOfMerit.Formula != "punzi" && config.FigureOfMerit.Formula != "asimov" && config.FigureOfMerit.Formula != "poisson") {
        throw std::runtime_error(std::format("fom: unknown formula \"{}\" (expected punzi, asimov or poisson)", config.FigureOfMerit.Formula));
    }
    config.FigureOfMerit.A = Optional<double>(fom, "a", 3.);
    config.FigureOfMerit.FSyst = Optional<double>(fom, "f_syst", 0.2);
    if (config.FigureOfMerit.FSyst < 0.) {
        throw std::runtime_error("fom: \"f_syst\" must not be negative");
    }
    config.FigureOfMerit.HasNSignalExpected = fom.contains("n_signal_expected");
    config.FigureOfMerit.NSignalExpected = Optional<double>(fom, "n_signal_expected", 0.);

    config.FigureOfMerit.HasSignalYield = fom.contains("signal_yield");
    if (config.FigureOfMerit.HasSignalYield) {
        const auto& yield = Require<json>(fom, "fom", "signal_yield");
        config.FigureOfMerit.DnDy = Require<double>(yield, "fom.signal_yield", "dndy");
        config.FigureOfMerit.DeltaY = Require<double>(yield, "fom.signal_yield", "delta_y");
        config.FigureOfMerit.BranchingRatio = Optional<double>(yield, "branching_ratio", 1.);
        config.FigureOfMerit.NInjectedSpecies = Optional<unsigned int>(yield, "n_injected_species", 1U);
        config.FigureOfMerit.InteractionProbability = Optional<double>(yield, "interaction_probability", 1.);
        config.FigureOfMerit.YieldSource = Optional<std::string>(yield, "source", std::string{});
        if (config.FigureOfMerit.DnDy <= 0. || config.FigureOfMerit.DeltaY <= 0.) {
            throw std::runtime_error("fom.signal_yield: \"dndy\" and \"delta_y\" must be positive");
        }
        if (config.FigureOfMerit.YieldSource.empty()) {
            throw std::runtime_error(
                "fom.signal_yield: \"source\" is required and must not be empty. These yields are pasted from a separate calculation, "
                "not computed here, and the macros that produced them call them provisional -- a number nobody can trace back is one "
                "nobody can check.");
        }
        // The sexaquark MC forces the sexaquark-nucleon interaction to happen, so its efficiency is
        // conditional on a probability that is never claimed. Leaving the factor at 1 quietly asserts the
        // interaction is certain.
        const bool is_sexaquark = config.Channel != EChannel::kHdibaryon;
        if (is_sexaquark && config.FigureOfMerit.Formula != "punzi" && config.FigureOfMerit.InteractionProbability == 1.) {
            throw std::runtime_error(
                std::format("fom.signal_yield: channel {} needs an explicit \"interaction_probability\" -- its MC "
                            "forces the sexaquark-nucleon interaction, so the efficiency is conditional on a "
                            "probability this config would otherwise assert to be 1.",
                            AsString(config.Channel)));
        }
    }

    if (config.FigureOfMerit.HasNSignalExpected && config.FigureOfMerit.HasSignalYield)
        throw std::runtime_error(
            "fom: \"n_signal_expected\" and \"signal_yield\" both given. One is the number, the other is the recipe for it -- keeping "
            "both invites them to disagree, so pick whichever you can defend.");
    if (config.FigureOfMerit.Formula != "punzi" && !config.FigureOfMerit.HasNSignalExpected && !config.FigureOfMerit.HasSignalYield)
        throw std::runtime_error(
            std::format("fom: formula \"{}\" needs an absolute signal yield, so \"n_signal_expected\" or \"signal_yield\" is required. Punzi is the "
                        "default precisely because it does not -- it needs only the signal efficiency, which is what you can actually measure "
                        "without assuming a production rate.",
                        config.FigureOfMerit.Formula));

    // -- guards

    const auto& guards = Optional<json>(json_file, "guards", json::object());
    config.Guard.MinRawSignal = Optional<unsigned int>(guards, "min_raw_signal", 20U);
    config.Guard.MinRawBackground = Optional<unsigned int>(guards, "min_raw_background", 20U);

    // -- acknowledged dummy populations

    const json sentinel_ok = Optional<json>(json_file, "sentinel_ok", json::array());
    for (const auto& node : sentinel_ok) {
        auto name = node.get<std::string>();
        const bool known = std::ranges::any_of(config.Variables, [&name](const CutVariable& v) { return v.Name == name; }) ||
                           std::ranges::any_of(config.Baseline, [&name](const BaselineCut& c) { return c.Variable == name; });
        if (!known)
            throw std::runtime_error(
                std::format("config: \"sentinel_ok\" names \"{}\", which is neither a scanned variable nor a "
                            "baseline cut. It exists to acknowledge a dummy population on a variable the "
                            "selection actually uses.",
                            name));
        config.SentinelOk.push_back(std::move(name));
    }

    return config;
}

}  // namespace Skimmer
