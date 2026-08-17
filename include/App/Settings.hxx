#pragma once

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "App/Logger.hxx"

namespace T2DS {

enum EProgramMode : std::uint8_t { FINDER, PACKAGER, VERIFIER };
inline constexpr std::array<const char *, 3> Name_ProgramMode{"FINDER", "PACKAGER", "VERIFIER"};

struct Settings {
    void Print() const {
        Logger::Info("Settings", "Mode            = {}", Name_ProgramMode[Mode]);
        Logger::Info("Settings", "Data Kind       = {}", IsMC ? "MC" : "RD");
        if (IsMC && Mode == EProgramMode::PACKAGER) {
            Logger::Info("Settings", "ReactionChannel = {}", ReactionChannel);
            Logger::Info("Settings", "SexaquarkMass   = {:.2f}", SexaquarkMass);
        } else if (Mode == EProgramMode::FINDER) {
            Logger::Info("Settings", "ReactionChannel = {}", ReactionChannel);
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

    std::vector<std::string> PathInputFiles;
    std::string PathOutputFile;
    std::optional<unsigned long> LimitToNEvents;
    double SexaquarkMass{};
    EProgramMode Mode{EProgramMode::PACKAGER};
    char ReactionChannel{};
    bool IsMC{false};
};

}  // namespace T2DS
