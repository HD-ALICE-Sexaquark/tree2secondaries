#ifndef T2S_SETTINGS_HXX
#define T2S_SETTINGS_HXX

#include <cctype>
#include <iostream>
#include <string>
#include <vector>

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
        std::cout << "SETTINGS" << '\n';
        std::cout << "========" << '\n';
        std::cout << "Mode             = " << (DoTheSearch ? "FINDER" : "PACKAGER") << '\n';
        std::cout << "ReactionChannel  = " << StrReactionChannel() << '\n';
        std::cout << "InputFiles       =" << '\n';
        for (const auto& path : PathInputFiles) {
            std::cout << "- " << path << '\n';
        }
        std::cout << "OutputFile       = " << PathOutputFile << '\n';
        std::cout << "IsMC             = " << IsMC << '\n';
        std::cout << "IsSignalMC       = " << IsSignalMC << '\n';
        if (IsSignalMC) std::cout << "   SexaquarkMass = " << SexaquarkMass << '\n';
        std::cout << "LimitToNEvents   = " << LimitToNEvents << '\n';
    }

    std::vector<std::string> PathInputFiles;
    std::string PathOutputFile;
    long long LimitToNEvents{0};
    double SexaquarkMass{Const::StandardSexaquarkMass};
    ReactionChannel Channel{ReactionChannel::All};
    bool IsMC{false};
    bool IsSignalMC{false};
    bool DoTheSearch{false};
};

}  // namespace Tree2Secondaries

#endif  // T2S_SETTINGS_HXX
