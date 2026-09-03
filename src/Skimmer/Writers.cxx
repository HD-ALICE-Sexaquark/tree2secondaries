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

#include "common/Framework.hpp"

#include "Skimmer/Writers.hxx"

namespace Skimmer {

// # Flat Cache # //

CacheWriter::CacheWriter(std::string_view ntuple_name, TFile& file, std::span<const std::string> fields) {

    // NOTE: the addresses of `fValues`' elements are handed to the entry, so this must be sized once up front and never resized again.
    fValues.assign(fields.size(), 0.F);

    Framework::Model model;
    model.RegisterField<std::uint8_t>(&fSampleIndex, "SampleIndex");
    model.RegisterField<std::uint8_t>(&fClassification, "Classification");
    model.RegisterField<std::uint8_t>(&fGeneratorMask, "GeneratorMask");
    model.RegisterField<unsigned int>(&fRunNumber, "RunNumber");
    model.RegisterField<unsigned int>(&fDirNumber, "DirNumber");
    model.RegisterField<unsigned int>(&fDirNumberB, "DirNumberB");
    model.RegisterField<unsigned int>(&fEventNumber, "EventNumber");
    model.RegisterField<float>(&fWeight, "Weight");
    for (std::size_t i = 0; i < fields.size(); ++i) {
        model.RegisterField<float>(&fValues[i], fields[i]);
    }

    fWriter = std::make_unique<Framework::Writer>(std::move(model), ntuple_name, file, ROOT::RNTupleWriteOptions());
}

void CacheWriter::Fill(std::uint8_t sample_index, Classification::EClassification classification, std::uint8_t generator_mask,
                       const POD::Event& event, float weight, std::span<const double> values) {
    fSampleIndex = sample_index;
    fClassification = static_cast<std::uint8_t>(classification);
    fGeneratorMask = generator_mask;
    fRunNumber = event.RunNumber;
    fDirNumber = event.DirNumber;
    fDirNumberB = event.DirNumberB;
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
    model.RegisterField<std::uint64_t>(&fNDropped_Truth, "NDropped_Truth");
    model.RegisterField<std::uint64_t>(&fNDropped_Origin, "NDropped_Origin");
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
    fNDropped_Truth = row.NDropped_Truth;
    fNDropped_Origin = row.NDropped_Origin;
    fWriter->Fill();
}

}  // namespace Skimmer
