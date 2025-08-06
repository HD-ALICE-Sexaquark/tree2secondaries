#ifndef T2S_SETTINGS_HXX
#define T2S_SETTINGS_HXX

#include <cctype>
#include <iostream>
#include <string>
#include <vector>

#include "App/Logger.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries {

struct Settings {
    std::string StrReactionChannel() const {
        if (Channel == ReactionChannel::All) return "AllChannels";
        std::string prefix{std::islower(static_cast<char>(Channel)) ? "Anti" : "Channel"};
        std::string suffix{static_cast<char>(std::toupper(static_cast<char>(Channel)))};
        return prefix + suffix;
    }
    void Print() const {
        Logger::Info("Settings", "Mode             = {}", (DoTheSearch ? "FINDER" : "PACKAGER"));
        Logger::Info("Settings", "ReactionChannel  = {}", StrReactionChannel());
        Logger::Info("Settings", "InputFiles       =");
        for (const auto& path : PathInputFiles) {
            Logger::Info("Settings", "- {}", path);
        }
        Logger::Info("Settings", "OutputFile       = {}", PathOutputFile);
        Logger::Info("Settings", "IsMC             = {}", IsMC);
        Logger::Info("Settings", "IsSignalMC       = {}", IsSignalMC);
        if (IsSignalMC) Logger::Info("Settings", "SexaquarkMass = {}", SexaquarkMass);
        Logger::Info("Settings", "LimitToNEvents   = {}", LimitToNEvents);
    }

    std::vector<std::string> PathInputFiles;
    std::string PathOutputFile;
    int LimitToNEvents{0};
    double SexaquarkMass{Const::StandardSexaquarkMass};
    ReactionChannel Channel{ReactionChannel::All};
    bool IsMC{false};
    bool IsSignalMC{false};
    bool DoTheSearch{false};
};

}  // namespace Tree2Secondaries

#endif  // T2S_SETTINGS_HXX
