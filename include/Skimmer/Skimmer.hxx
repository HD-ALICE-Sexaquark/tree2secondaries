#pragma once

#include <cstdint>
#include <string_view>

#include "Skimmer/Config.hxx"

namespace Skimmer {

void Run(const Config& config, std::string_view output_dir, std::uint64_t n_events_limit);

}
