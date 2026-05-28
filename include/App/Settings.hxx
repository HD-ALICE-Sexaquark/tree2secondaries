#pragma once

#include <string>
#include <vector>

#include "common/Constants.hpp"
#include "common/DB_ReactionChannels.hpp"

#include "App/Logger.hxx"

namespace R2DS {

enum class EProgramMode { FINDER, PACKAGER, VERIFIER };

struct Settings {
    void Print() const {
        Logger::Info("Settings", "Mode            = {}",
                     Mode == EProgramMode::FINDER ? "FINDER" : (Mode == EProgramMode::PACKAGER ? "PACKAGER" : "VERIFIER"));
        Logger::Info("Settings", "ReactionChannel = {}", ReactionChannel.name);
        Logger::Info("Settings", "InputFiles      = ");
        for (const auto& path : PathInputFiles) {
            Logger::Info("Settings", "- {}", path);
        }
        Logger::Info("Settings", "IsMC            = {}", IsMC);
        Logger::Info("Settings", "OutputFile      = {}", PathOutputFile);
        Logger::Info("Settings", "SexaquarkMass   = {}", SexaquarkMass);
        Logger::Info("Settings", "LimitToNEvents  = {}", LimitToNEvents);
    }

    std::vector<std::string> PathInputFiles;
    std::string PathOutputFile;
    DB::ReactionChannels::Definition ReactionChannel;
    double SexaquarkMass{Common::DummyDouble};
    long long LimitToNEvents{0};
    EProgramMode Mode{EProgramMode::PACKAGER};
    bool IsMC{false};
};

}  // namespace R2DS
