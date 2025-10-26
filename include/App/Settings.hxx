#ifndef T2S_SETTINGS_HXX
#define T2S_SETTINGS_HXX

#include <cctype>
#include <string>
#include <vector>

#include "App/Logger.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries {

struct Settings {
    void Print() const {
        Logger::Info("Settings", "Mode            = {}", (DoTheSearch ? "FINDER" : "PACKAGER"));
        Logger::Info("Settings", "ReactionChannel = {}", static_cast<char>(ReactionChannel));
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
    int LimitToNEvents{0};
    double SexaquarkMass{Const::StandardSexaquarkMass};
    EReactionChannel ReactionChannel{EReactionChannel::A};
    bool IsMC{false};
    bool DoTheSearch{false};
};

}  // namespace Tree2Secondaries

#endif  // T2S_SETTINGS_HXX
