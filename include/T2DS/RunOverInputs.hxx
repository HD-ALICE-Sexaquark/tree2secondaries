#pragma once

#include <algorithm>
#include <limits>

#include "T2DS/Settings.hxx"

namespace T2DS {

// Read every event of every input file, calling `process_event` on each.
// Possible `Worker`: `Finder`, `Verifier`.
template <typename Worker, typename ProcessEvent>
void RunOverInputs(Worker &worker, const T2DS::Settings &settings, const ProcessEvent &process_event) {

    // `remaining` clamps the total number of events to be read from *all* files, not per file
    auto remaining = settings.LimitToNEvents.has_value() ? settings.LimitToNEvents.value() : std::numeric_limits<unsigned long>::max();

    for (const auto &input_path : settings.PathInputFiles) {
        if (remaining == 0) break;
        if (!worker.OpenInput(input_path)) continue;  // unreadable files are skipped

        const auto n_events = std::min(worker.NumberEventsToRead(), remaining);
        for (unsigned long id_event = 0; id_event < n_events; ++id_event) {
            worker.Load(id_event);
            process_event();
        }
        remaining -= n_events;
    }
}

}  // namespace T2DS
