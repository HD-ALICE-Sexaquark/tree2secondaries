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

#include "Skimmer/Process.hxx"

namespace Skimmer {

struct Config;

// # The Flat Cache # //

// The bookkeeping columns every cache carries, whatever the channel. A config variable sharing one of
// these names would register a duplicate RNTuple field, so `BuildPlan` rejects it at startup rather
// than letting the writer fail halfway through a production.
inline constexpr auto kReservedFields = std::to_array<std::string_view>({"Process", "RunNumber", "EventNumber", "Weight"});

class CacheWriter {
   public:
    CacheWriter(const CacheWriter&) = delete;
    CacheWriter(CacheWriter&&) = delete;
    CacheWriter& operator=(const CacheWriter&) = delete;
    CacheWriter& operator=(CacheWriter&&) = delete;
    ~CacheWriter() = default;

    CacheWriter(std::string_view ntuple_name, TFile& file, std::span<const std::string> fields);

    void Fill(Process::EProcess process, unsigned int run_number, unsigned int event_number, float weight, std::span<const double> values);

   private:
    std::vector<float> fValues;
    std::uint8_t fProcess{};
    unsigned int fRunNumber{};
    unsigned int fEventNumber{};
    float fWeight{};
    std::unique_ptr<Framework::Writer> fWriter;
};

// # Bookkeeping # //

// Metadata: one row per input file.
class MetaWriter {
   public:
    MetaWriter(const MetaWriter&) = delete;
    MetaWriter(MetaWriter&&) = delete;
    MetaWriter& operator=(const MetaWriter&) = delete;
    MetaWriter& operator=(MetaWriter&&) = delete;
    ~MetaWriter() = default;

    explicit MetaWriter(TFile& file);

    void Fill(std::string_view path, unsigned int run_number, std::uint64_t n_events, std::uint64_t n_events_read, std::uint64_t n_read,
              std::uint64_t n_written);

   private:
    std::string fPath;
    unsigned int fRunNumber{};
    std::uint64_t fNEvents{};
    std::uint64_t fNEventsRead{};
    std::uint64_t fNCandidatesRead{};
    std::uint64_t fNCandidatesWritten{};
    std::unique_ptr<Framework::Writer> fWriter;
};

void WriteProvenance(TFile& file, const Config& config, unsigned int n_injected_per_event, bool is_partial);

}  // namespace Skimmer
