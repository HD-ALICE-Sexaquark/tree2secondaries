#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <TFile.h>

#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleWriteOptions.hxx>

#include "common/Constants.hpp"
#include "common/Framework.hpp"

#include "Skimmer/Config.hxx"
#include "Skimmer/Writers.hxx"

namespace Skimmer {

// # Flat Cache # //

CacheWriter::CacheWriter(std::string_view ntuple_name, TFile& file, std::span<const std::string> fields) {

    // NOTE: the addresses of `fValues`' elements are handed to the entry, so this must be sized once up front and never resized again.
    fValues.assign(fields.size(), 0.F);

    Framework::Model model;
    model.RegisterField<std::uint8_t>(&fSampleIndex, "SampleIndex");
    model.RegisterField<std::uint8_t>(&fProcess, "Process");
    model.RegisterField<unsigned int>(&fRunNumber, "RunNumber");
    model.RegisterField<unsigned int>(&fEventNumber, "EventNumber");
    model.RegisterField<float>(&fWeight, "Weight");
    for (std::size_t i = 0; i < fields.size(); ++i) {
        model.RegisterField<float>(&fValues[i], fields[i]);
    }

    fWriter = std::make_unique<Framework::Writer>(std::move(model), ntuple_name, file, ROOT::RNTupleWriteOptions());
}

void CacheWriter::Fill(std::uint8_t sample_index, Process::EProcess process, unsigned int run_number, unsigned int event_number, float weight,
                       std::span<const double> values) {
    fSampleIndex = sample_index;
    fProcess = static_cast<std::uint8_t>(process);
    fRunNumber = run_number;
    fEventNumber = event_number;
    fWeight = weight;
    for (std::size_t i = 0; i < fValues.size(); ++i) fValues[i] = static_cast<float>(values[i]);
    fWriter->Fill();
}

// # Bookkeeping # //

MetaWriter::MetaWriter(TFile& file) {
    Framework::Model model;
    model.RegisterField<std::uint8_t>(&fSampleIndex, "SampleIndex");
    model.RegisterField<std::string>(&fPath, "Path");
    model.RegisterField<unsigned int>(&fRunNumber, "RunNumber");
    model.RegisterField<std::uint8_t>(&fRole, "Role");
    model.RegisterField<std::uint32_t>(&fNInjectedPerEvent, "NInjectedPerEvent");
    model.RegisterField<std::uint64_t>(&fNEvents, "NEvents");
    model.RegisterField<std::uint64_t>(&fNEventsRead, "NEventsRead");
    model.RegisterField<std::uint64_t>(&fNCandidatesRead, "NCandidatesRead");
    model.RegisterField<std::uint64_t>(&fNCandidatesWritten, "NCandidatesWritten");
    fWriter = std::make_unique<Framework::Writer>(std::move(model), "Meta", file, ROOT::RNTupleWriteOptions());
}

void MetaWriter::Fill(const MetaRow& row) {
    fSampleIndex = row.SampleIndex;
    fPath = std::string(row.Path);
    fRunNumber = row.RunNumber;
    fRole = static_cast<std::uint8_t>(row.Role);
    fNInjectedPerEvent = row.NInjectedPerEvent;
    fNEvents = row.NEvents;
    fNEventsRead = row.NEventsRead;
    fNCandidatesRead = row.NCandidatesRead;
    fNCandidatesWritten = row.NCandidatesWritten;
    fWriter->Fill();
}

void WriteProvenance(TFile& file, const Skimmer::Config& config, bool is_partial) {

    std::string channel{AsString(config.Channel)};
    std::string hypothesis{config.Hypothesis};
    std::string config_path{config.Path};
    std::string observable{config.Obs.Variable};
    double signal_mass{config.SignalMass};
    std::uint64_t expected_real_events{E2T::NExpectedEventsInRealData};

    Framework::Model model;
    model.RegisterField<std::string>(&channel, "Channel");
    model.RegisterField<std::string>(&hypothesis, "Hypothesis");
    model.RegisterField<std::string>(&config_path, "ConfigPath");
    model.RegisterField<std::string>(&observable, "Observable");
    model.RegisterField<double>(&signal_mass, "SignalMass");
    model.RegisterField<std::uint64_t>(&expected_real_events, "NExpectedEventsInRealData");
    model.RegisterField<bool>(&is_partial, "IsPartial");

    Framework::Writer writer{std::move(model), "Provenance", file, ROOT::RNTupleWriteOptions()};
    writer.Fill();
}

}  // namespace Skimmer
