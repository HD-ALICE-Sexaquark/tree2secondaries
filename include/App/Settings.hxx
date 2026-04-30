#pragma once

#include <string>
#include <vector>

#include "App/DB_ReactionChannels.hxx"
#include "App/Logger.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries {

struct Settings {
    void Print() const {
        Logger::Info("Settings", "Mode            = {}", DoTheSearch ? "FINDER" : "PACKAGER");
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
    double SexaquarkMass{Const::StandardSexaquarkMass};
    long long LimitToNEvents{0};
    ReactionChannels::Definition ReactionChannel;
    bool IsMC{false};
    bool DoTheSearch{false};
};

}  // namespace Tree2Secondaries
