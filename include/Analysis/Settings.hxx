#ifndef T2S_SETTINGS_HXX
#define T2S_SETTINGS_HXX

#include <string>

#include "Utilities/Logger.hxx"

namespace Tree2Secondaries {

struct Settings {
    void Print() const {
        INFO("SETTINGS");
        INFO("========");
        INFO("IsMC           = %i", IsMC);
        INFO("IsSignalMC     = %i", IsSignalMC);
        INFO("InputFiles     = %s", PathInputFiles.c_str());
        INFO("OutputFile     = %s", PathOutputFile.c_str());
        INFO("LimitToNEvents = %lld", LimitToNEvents);
    }

    std::string PathInputFiles;
    std::string PathOutputFile;
    bool IsMC;
    bool IsSignalMC;
    long long LimitToNEvents;
};

}  // namespace Tree2Secondaries

#endif  // T2S_SETTINGS_HXX
