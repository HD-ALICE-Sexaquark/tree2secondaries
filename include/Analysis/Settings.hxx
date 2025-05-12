#ifndef T2S_SETTINGS_HXX
#define T2S_SETTINGS_HXX

#include <cctype>
#include <string>

#include "Math/Constants.hxx"
#include "Utilities/Logger.hxx"

namespace Tree2Secondaries {

struct Settings {
    [[nodiscard]] std::string StrReactionChannel() const {
        std::string prefix{std::islower(static_cast<char>(Channel)) ? "Anti" : ""};
        std::string suffix{static_cast<char>(std::toupper(static_cast<char>(Channel)))};
        return prefix + suffix;
    }
    void Print() const {
        INFO("SETTINGS");
        INFO("========");
        INFO("ReactionChannel = %s", StrReactionChannel().c_str());
        INFO("InputFiles      = %s", PathInputFiles.c_str());
        INFO("OutputFile      = %s", PathOutputFile.c_str());
        INFO("IsMC            = %i", IsMC);
        INFO("IsSignalMC      = %i", IsSignalMC);
        INFO("LimitToNEvents  = %lld", LimitToNEvents);
    }

    std::string PathInputFiles;
    std::string PathOutputFile;
    long long LimitToNEvents{0};
    ReactionChannel Channel{ReactionChannel::A};
    bool IsMC{false};
    bool IsSignalMC{false};
};

}  // namespace Tree2Secondaries

#endif  // T2S_SETTINGS_HXX
