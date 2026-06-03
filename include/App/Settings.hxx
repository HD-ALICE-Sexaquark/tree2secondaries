#pragma once

#include <cstdint>
#include <optional>
#include <string>

#include "App/Logger.hxx"

namespace R2DS {

enum EProgramMode : std::uint8_t { FINDER, PACKAGER, VERIFIER };
const std::array Name_ProgramMode{"FINDER", "PACKAGER", "VERIFIER"};

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
        Logger::Info("Settings", "InputFile       = {}", PathInputFile);
        Logger::Info("Settings", "OutputFile      = {}", PathOutputFile);
        if (LimitToNEvents.has_value()) {
            Logger::Info("Settings", "LimitToNEvents  = {}", LimitToNEvents.value());
        } else {
            Logger::Info("Settings", "LimitToNEvents  = --");
        }
    }

    std::string PathInputFile;
    std::string PathOutputFile;
    std::optional<unsigned long> LimitToNEvents;
    double SexaquarkMass{};
    EProgramMode Mode{EProgramMode::PACKAGER};
    char ReactionChannel{};
    bool IsMC;
};

}  // namespace R2DS
