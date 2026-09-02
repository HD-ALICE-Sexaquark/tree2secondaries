#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <format>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <TFile.h>
#include <TH1.h>

#include <Math/Point3D.h>
namespace RMath = ROOT::Math;

#include "Skimmer/Channels.hxx"
#include "Skimmer/Config.hxx"
#include "Skimmer/Writers.hxx"
#include "Utils/Logger.hxx"

#include "Skimmer/Skimmer.hxx"

namespace {

// # The Skim Plan # //

struct BaselineSpec {
    std::size_t FieldIndex;
    Skimmer::EDirection Direction;
    double Value;       // the threshold for `kLower`/`kUpper`, the lower edge for `kWindow`
    double ValueUpper;  // only meaningful for `kWindow`
};

struct SkimPlan {
    std::vector<std::string> Fields;
    std::vector<std::size_t> RegistryIndex;  // parallel to `Fields`
    std::vector<BaselineSpec> Baseline;
};

// Check if every variable named in the config exists in the registry for its channel.
template <typename Traits>
SkimPlan BuildPlan(const Skimmer::Config& config) {

    SkimPlan plan;
    const std::vector<std::string> names = config.CachedVariables();

    for (const auto& name : names) {
        // Checked before the registry lookup so that the reserved name gets the specific message,
        // whether it arrived from the config or from a future registry entry.
        if (std::ranges::find(Skimmer::kReservedFields, name) != Skimmer::kReservedFields.end()) {
            throw std::runtime_error(std::format("variable \"{}\" collides with a reserved cache field (see Skimmer/Writers.hxx)", name));
        }
        const auto index = Skimmer::Variables::FindIndex(Traits::kVariables, name);
        if (!index) {
            throw std::runtime_error(
                std::format("variable \"{}\" does not exist in the registry for channel {} (see Skimmer/VariableRegistry.hxx)",  //
                            name, Skimmer::AsString(config.Channel)));
        }
        plan.RegistryIndex.push_back(*index);
        plan.Fields.push_back(name);
    }

    for (const auto& cut : config.Baseline) {
        const auto it = std::ranges::find(names, cut.Variable);
        if (it == names.end()) {
            // unreachable via `CachedVariables()`, which includes every baseline variable by construction
            throw std::runtime_error(std::format("baseline variable \"{}\" was not registered in the cache", cut.Variable));
        }
        plan.Baseline.push_back({static_cast<std::size_t>(it - names.begin()), cut.Direction, cut.Value, cut.ValueUpper});
    }

    return plan;
}

[[nodiscard]] bool PassesBaseline(const SkimPlan& plan, std::span<const double> values) {
    for (const auto& cut : plan.Baseline) {
        const double value = values[cut.FieldIndex];
        switch (cut.Direction) {
            case Skimmer::EDirection::kLower:
                if (value < cut.Value) return false;
                break;
            case Skimmer::EDirection::kUpper:
                if (value > cut.Value) return false;
                break;
            case Skimmer::EDirection::kWindow:
                if (value < cut.Value || value > cut.ValueUpper) return false;
                break;
        }
    }
    return true;
}

// # Where the Output Goes # //

std::filesystem::path ResolveOutput(const Skimmer::Config& config, std::string_view output_dir) {
    std::filesystem::path output{output_dir};
    if (output.extension() == ".root") {
        std::filesystem::create_directories(output.parent_path().empty() ? std::filesystem::path{"."} : output.parent_path());
        return output;
    }
    std::filesystem::create_directories(output);
    return output / std::format("Cache_{}.root", Skimmer::AsString(config.Channel));
}

// # The Event Count # //

// `t2ds` writes an `N_Events` counter beside each RNTuple, filled once per event at the top of `Finder::ProcessEvent`, i.e. before any candidate is
// built.
// NOTE: it's the correct event denominator, because `t2ds` drops events that produced neither a candidate nor an injection; so the entry count
// of the RNTuple is always the smaller number.

[[nodiscard]] std::uint64_t ReadNEvents(std::string_view path) {

    const std::unique_ptr<TFile> file{TFile::Open(std::string(path).c_str(), "READ")};
    if (!file || file->IsZombie()) {
        throw std::runtime_error(std::format("could not open \"{}\"", path));
    }

    const auto* histogram = file->Get<TH1>("N_Events");
    if (!histogram) {
        // Returning 0 here would zero this sample's denominator without zeroing its candidates, which no
        // consumer can detect afterwards.
        throw std::runtime_error(
            std::format("\"{}\" has no N_Events histogram -- it was not produced by `t2ds`, so it has no normalization denominator", path));
    }
    return static_cast<std::uint64_t>(histogram->GetEntries());
}

// # Implementation per Channel Trait # //

template <typename Traits>
void RunChannel(const Skimmer::Config& config, std::string_view output_dir, std::uint64_t n_events_limit) {

    const SkimPlan plan = BuildPlan<Traits>(config);

    const std::filesystem::path output_path = ResolveOutput(config, output_dir);
    const std::unique_ptr<TFile> output_file{TFile::Open(output_path.c_str(), "RECREATE")};
    if (!output_file || output_file->IsZombie()) {
        throw std::runtime_error(std::format("could not create \"{}\"", output_path.string()));
    }

    std::uint64_t n_events_total{0};
    std::uint64_t n_events_read_total{0};
    std::uint64_t n_read_total{0};
    std::uint64_t n_written_total{0};
    std::uint64_t n_injected_total{0};

    // NOTE: the writers have to be destroyed before the file is closed -- an `RNTupleWriter` commits its
    //       final cluster and footer in its destructor. Hence the scope.
    {
        Skimmer::CacheWriter cache{Traits::kName_OutputRNT, *output_file, plan.Fields};
        Skimmer::MetaWriter meta{*output_file};

        std::vector<double> values(plan.Fields.size(), 0.);

        for (std::size_t s = 0; s < config.Samples.size(); ++s) {

            const Skimmer::Sample& sample = config.Samples[s];
            const auto sample_index = static_cast<std::uint8_t>(s);

            const std::uint64_t n_events = ReadNEvents(sample.Path);

            typename Traits::Schema schema;
            Framework::Reader reader{schema.CreateModel(true), Traits::kName_InputRNT, sample.Path};

            const auto n_entries = reader.Iter()->GetNEntries();
            std::uint64_t n_events_read{0};
            std::uint64_t n_injected{0};
            std::uint64_t n_read{0};
            std::uint64_t n_written{0};

            for (ROOT::NTupleSize_t entry = 0; entry < n_entries; ++entry) {

                if (n_events_limit > 0 && n_events_read >= n_events_limit) break;
                reader.Load(entry);
                ++n_events_read;

                n_injected += schema.Injected.size();

                // `Traits::Label` indexes the MC vector with an index bounded by the reconstructed one.
                if (Traits::McSize(schema) != Traits::Size(schema)) {
                    throw std::runtime_error(std::format("{}: entry {} has {} candidates but {} linked MC entries -- the vectors are not parallel",
                                                         sample.Path, entry, Traits::Size(schema), Traits::McSize(schema)));
                }

                const RMath::XYZPoint pv{schema.Event.PV_X, schema.Event.PV_Y, schema.Event.PV_Z};

                const std::size_t n_candidates = Traits::Size(schema);
                n_read += n_candidates;

                for (std::size_t i = 0; i < n_candidates; ++i) {

                    const typename Traits::Cached cached = Traits::Build(schema, i, pv);
                    for (std::size_t f = 0; f < plan.RegistryIndex.size(); ++f) values[f] = Traits::kVariables[plan.RegistryIndex[f]].Extract(cached);
                    if (!PassesBaseline(plan, values)) continue;

                    cache.Fill(sample_index, Traits::Label(schema, i, cached), schema.Event, 1.F, values);
                    ++n_written;
                }
            }

            meta.Fill({.SampleIndex = sample_index,
                       .Path = sample.Path,
                       .Role = sample.Role,
                       .NInjected = n_injected,
                       .NEvents = n_events,
                       .NEventsRead = n_events_read,
                       .NCandidatesRead = n_read,
                       .NCandidatesWritten = n_written});

            Logger::Info(__FUNCTION__, "{} ({}): {} events, {} candidates read, {} written ({:.2f}%)", sample.Path, Skimmer::AsString(sample.Role),
                         n_events_read, n_read, n_written, n_read > 0 ? 100. * static_cast<double>(n_written) / static_cast<double>(n_read) : 0.);

            n_events_total += n_events;
            n_events_read_total += n_events_read;
            n_read_total += n_read;
            n_written_total += n_written;
            n_injected_total += n_injected;
        }
    }

    Skimmer::WriteCacheSource(*output_file, config, Traits::kName_OutputRNT, plan.Fields, n_events_limit);

    output_file->Close();

    Logger::Info(__FUNCTION__, "wrote {} rows to \"{}\" ({} candidates read from {} events)", n_written_total, output_path.string(), n_read_total,
                 n_events_read_total);
    Logger::Info(__FUNCTION__, "the inputs hold {} events; {} read, holding {} injected signals",  //
                 n_events_total, n_events_read_total, n_injected_total);

    if (n_events_limit > 0) {
        Logger::Info(__FUNCTION__,
                     "this is a partial skim ({} events per sample at most); it is flagged as such and must not be used for normalization",
                     n_events_limit);
    }
}

}  // namespace

// # Entry Point # //

namespace Skimmer {

void Run(const Config& config, std::string_view output_dir, std::uint64_t n_events_limit) {
    switch (config.Channel) {
        case EChannel::kChannelA:
            RunChannel<TraitsChannelA>(config, output_dir, n_events_limit);
            break;
        case EChannel::kChannelD:
            RunChannel<TraitsChannelD>(config, output_dir, n_events_limit);
            break;
        case EChannel::kChannelH:
            RunChannel<TraitsChannelH>(config, output_dir, n_events_limit);
            break;
        case EChannel::kLambdaPair:
            RunChannel<TraitsLambdaPair>(config, output_dir, n_events_limit);
            break;
    }
}

}  // namespace Skimmer
