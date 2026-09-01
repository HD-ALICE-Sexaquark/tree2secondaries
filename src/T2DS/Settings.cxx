#include <array>
#include <optional>
#include <vector>

#include "Utils/Logger.hxx"

#include "T2DS/Settings.hxx"

namespace T2DS {

void Settings::Print() const {
    Logger::Info("Settings", "Mode            = {}", Name_ProgramMode[Mode]);
    Logger::Info("Settings", "Data Type       = {}", IsMC ? "MC" : "RD");
    if (Mode == EProgramMode::FINDER && IsMC) {
        Logger::Info("Settings", "SexaquarkMass   = {:.2f}", SexaquarkMass);
    }
    if (PathInputFiles.size() == 1) {
        Logger::Info("Settings", "InputFile       = {}", PathInputFiles.front());
    } else {
        Logger::Info("Settings", "InputFiles      = {} files", PathInputFiles.size());
        for (const auto &file : PathInputFiles) {
            Logger::Info("Settings", "-- {}", file);
        }
    }
    Logger::Info("Settings", "OutputFile      = {}", PathOutputFile);
    if (LimitToNEvents.has_value()) {
        Logger::Info("Settings", "LimitToNEvents  = {}", LimitToNEvents.value());
    } else {
        Logger::Info("Settings", "LimitToNEvents  = --");
    }
}

}  // namespace T2DS
