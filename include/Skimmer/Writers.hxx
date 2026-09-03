#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include <TFile.h>

#include "common/Framework.hpp"
#include "common/POD_Event.hpp"

#include "Skimmer/Classification.hxx"
#include "Skimmer/Config.hxx"

namespace Skimmer {

// # The Flat Cache # //

// The bookkeeping columns every cache carries, whatever the channel. A config variable sharing one of
// these names would register a duplicate RNTuple field, so `BuildPlan` rejects it at startup rather
// than letting the writer fail halfway through a production.
inline constexpr auto kReservedFields = std::to_array<std::string_view>(
    {"SampleIndex", "Classification", "GeneratorMask", "RunNumber", "DirNumber", "DirNumberB", "EventNumber", "Weight"});

class CacheWriter {
   public:
    CacheWriter(const CacheWriter&) = delete;
    CacheWriter(CacheWriter&&) = delete;
    CacheWriter& operator=(const CacheWriter&) = delete;
    CacheWriter& operator=(CacheWriter&&) = delete;
    ~CacheWriter() = default;

    CacheWriter(std::string_view ntuple_name, TFile& file, std::span<const std::string> fields);

    // `sample_index` is the row's position in the config's `samples[]`, i.e. the `Meta` row it belongs to.
    // Without it a row cannot be traced back to its file, since a file may span several runs.
    // The event is taken whole rather than field by field, so that the identity columns and `kReservedFields`
    // only ever have to agree with each other in one place.
    // `generator_mask` is stored rather than consumed, so that a consumer can still cut on the generator of origin
    // without re-running the skim -- and so that an origin filter applied here stays visible after the fact.
    void Fill(std::uint8_t sample_index, Classification::EClassification classification, std::uint8_t generator_mask, const POD::Event& event,
              float weight, std::span<const double> values);

   private:
    std::vector<float> fValues;
    std::uint8_t fSampleIndex{};
    std::uint8_t fClassification{};
    std::uint8_t fGeneratorMask{};
    unsigned int fRunNumber{};
    unsigned int fDirNumber{};
    unsigned int fDirNumberB{};
    unsigned int fEventNumber{};
    float fWeight{};
    std::unique_ptr<Framework::Writer> fWriter;
};

// # Bookkeeping # //

struct MetaRow {
    std::uint8_t SampleIndex{};
    std::string_view Path;
    ERole Role{ERole::kBoth};
    std::uint64_t NInjected{};
    std::uint64_t NEvents{};
    std::uint64_t NEventsRead{};
    std::uint64_t NCandidatesRead{};
    std::uint64_t NCandidatesWritten{};
    std::uint64_t NDropped_Truth{};   // dropped by `Classification::ShouldBeKept`: hybrids always, reference unless kept
    std::uint64_t NDropped_Origin{};  // dropped by the injected-background veto
};

// Metadata: one row per input file.
class MetaWriter {
   public:
    MetaWriter(const MetaWriter&) = delete;
    MetaWriter(MetaWriter&&) = delete;
    MetaWriter& operator=(const MetaWriter&) = delete;
    MetaWriter& operator=(MetaWriter&&) = delete;
    ~MetaWriter() = default;

    explicit MetaWriter(TFile& file);

    void Fill(const MetaRow& row);

   private:
    std::uint8_t fSampleIndex{};
    std::string fPath;
    std::uint8_t fRole{};
    std::uint64_t fNInjected{};
    std::uint64_t fNEvents{};
    std::uint64_t fNEventsRead{};
    std::uint64_t fNCandidatesRead{};
    std::uint64_t fNCandidatesWritten{};
    std::uint64_t fNDropped_Truth{};
    std::uint64_t fNDropped_Origin{};
    std::unique_ptr<Framework::Writer> fWriter;
};

}  // namespace Skimmer
