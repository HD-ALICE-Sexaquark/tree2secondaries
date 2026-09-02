#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <TFile.h>
#include <TObjString.h>

#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleWriteOptions.hxx>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

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
    model.RegisterField<std::uint8_t>(&fClassification, "Classification");
    model.RegisterField<unsigned int>(&fRunNumber, "RunNumber");
    model.RegisterField<unsigned int>(&fDirNumber, "DirNumber");
    model.RegisterField<unsigned int>(&fEventNumber, "EventNumber");
    model.RegisterField<float>(&fWeight, "Weight");
    for (std::size_t i = 0; i < fields.size(); ++i) {
        model.RegisterField<float>(&fValues[i], fields[i]);
    }

    fWriter = std::make_unique<Framework::Writer>(std::move(model), ntuple_name, file, ROOT::RNTupleWriteOptions());
}

void CacheWriter::Fill(std::uint8_t sample_index, Classification::EClassification classification, const POD::Event& event, float weight,
                       std::span<const double> values) {
    fSampleIndex = sample_index;
    fClassification = static_cast<std::uint8_t>(classification);
    fRunNumber = event.RunNumber;
    fDirNumber = event.DirNumber;
    fEventNumber = event.EventNumber;
    fWeight = weight;
    for (std::size_t i = 0; i < fValues.size(); ++i) fValues[i] = static_cast<float>(values[i]);
    fWriter->Fill();
}

// # Bookkeeping # //

MetaWriter::MetaWriter(TFile& file) {
    Framework::Model model;
    model.RegisterField<std::uint8_t>(&fSampleIndex, "SampleIndex");
    model.RegisterField<std::string>(&fPath, "Path");
    model.RegisterField<std::uint8_t>(&fRole, "Role");
    model.RegisterField<std::uint64_t>(&fNInjected, "NInjected");
    model.RegisterField<std::uint64_t>(&fNEvents, "NEvents");
    model.RegisterField<std::uint64_t>(&fNEventsRead, "NEventsRead");
    model.RegisterField<std::uint64_t>(&fNCandidatesRead, "NCandidatesRead");
    model.RegisterField<std::uint64_t>(&fNCandidatesWritten, "NCandidatesWritten");
    fWriter = std::make_unique<Framework::Writer>(std::move(model), "Meta", file, ROOT::RNTupleWriteOptions());
}

void MetaWriter::Fill(const MetaRow& row) {
    fSampleIndex = row.SampleIndex;
    fPath = std::string(row.Path);
    fRole = static_cast<std::uint8_t>(row.Role);
    fNInjected = row.NInjected;
    fNEvents = row.NEvents;
    fNEventsRead = row.NEventsRead;
    fNCandidatesRead = row.NCandidatesRead;
    fNCandidatesWritten = row.NCandidatesWritten;
    fWriter->Fill();
}

// # CacheSource # //

void WriteCacheSource(TFile& file, const Skimmer::Config& config, std::string_view ntuple_name, std::span<const std::string> fields,
                      std::uint64_t n_events_limit) {

    json record;
    record["channel"] = std::string(AsString(config.Channel));
    record["signal_mass"] = config.SignalMass;
    record["ntuple"] = std::string(ntuple_name);
    record["fields"] = std::vector<std::string>(fields.begin(), fields.end());
    record["config_path"] = config.Path;
    record["is_partial"] = n_events_limit > 0;
    record["n_events_limit"] = n_events_limit;

    TObjString text{record.dump(2).c_str()};
    file.WriteObject(&text, "CacheSource");
}

}  // namespace Skimmer
