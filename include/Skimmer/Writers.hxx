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

#include "Skimmer/Config.hxx"
#include "Skimmer/Process.hxx"

namespace Skimmer {

// # The Flat Cache # //

// The bookkeeping columns every cache carries, whatever the channel. A config variable sharing one of
// these names would register a duplicate RNTuple field, so `BuildPlan` rejects it at startup rather
// than letting the writer fail halfway through a production.
inline constexpr auto kReservedFields = std::to_array<std::string_view>({"SampleIndex", "Process", "RunNumber", "EventNumber", "Weight"});

class CacheWriter {
   public:
    CacheWriter(const CacheWriter&) = delete;
    CacheWriter(CacheWriter&&) = delete;
    CacheWriter& operator=(const CacheWriter&) = delete;
    CacheWriter& operator=(CacheWriter&&) = delete;
    ~CacheWriter() = default;

    CacheWriter(std::string_view ntuple_name, TFile& file, std::span<const std::string> fields);

    // `sample_index` is the row's position in the config's `samples[]`, i.e. the `Meta` row it belongs to.
    // Without it a row cannot be traced back to its file, since two samples may share (or omit) a run number.
    void Fill(std::uint8_t sample_index, Process::EProcess process, unsigned int run_number, unsigned int event_number, float weight,
              std::span<const double> values);

   private:
    std::vector<float> fValues;
    std::uint8_t fSampleIndex{};
    std::uint8_t fProcess{};
    unsigned int fRunNumber{};
    unsigned int fEventNumber{};
    float fWeight{};
    std::unique_ptr<Framework::Writer> fWriter;
};

// # Bookkeeping # //

struct MetaRow {
    std::uint8_t SampleIndex{};
    std::string_view Path;
    unsigned int RunNumber{};
    ERole Role{ERole::kBoth};
    unsigned int NInjectedPerEvent{};
    std::uint64_t NEvents{};
    std::uint64_t NEventsRead{};
    std::uint64_t NCandidatesRead{};
    std::uint64_t NCandidatesWritten{};
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
    unsigned int fRunNumber{};
    std::uint8_t fRole{};
    std::uint32_t fNInjectedPerEvent{};
    std::uint64_t fNEvents{};
    std::uint64_t fNEventsRead{};
    std::uint64_t fNCandidatesRead{};
    std::uint64_t fNCandidatesWritten{};
    std::unique_ptr<Framework::Writer> fWriter;
};

void WriteProvenance(TFile& file, const Config& config, bool is_partial);

}  // namespace Skimmer
