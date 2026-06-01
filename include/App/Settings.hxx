#pragma once

#include <optional>
#include <string>

#include "common/DB_ReactionChannels.hpp"

#include "App/Logger.hxx"

namespace R2DS {

enum class EProgramMode { FINDER, PACKAGER, VERIFIER };

struct Settings {
    void Print() const {
        Logger::Info("Settings", "Mode            = {}",
                     Mode == EProgramMode::FINDER ? "FINDER" : (Mode == EProgramMode::PACKAGER ? "PACKAGER" : "VERIFIER"));
        if (ReactionChannel.has_value()) {
            Logger::Info("Settings", "ReactionChannel = {}", ReactionChannel.value().name);
        }
        Logger::Info("Settings", "InputFile       = {}", PathInputFile);
        Logger::Info("Settings", "OutputFile      = {}", PathOutputFile);
        if (SexaquarkMass.has_value()) {
            Logger::Info("Settings", "SexaquarkMass   = {:.2f}", SexaquarkMass.value());
        }
        if (LimitToNEvents.has_value()) {
            Logger::Info("Settings", "LimitToNEvents  = {}", LimitToNEvents.value());
        } else {
            Logger::Info("Settings", "LimitToNEvents  = --");
        }
    }

    std::string PathInputFile;
    std::string PathOutputFile;
    std::optional<DB::ReactionChannels::Definition> ReactionChannel;
    std::optional<double> SexaquarkMass;
    std::optional<unsigned long> LimitToNEvents;
    EProgramMode Mode{EProgramMode::PACKAGER};
};

}  // namespace R2DS
