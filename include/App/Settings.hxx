#ifndef T2S_SETTINGS_HXX
#define T2S_SETTINGS_HXX

#include <string>

#include "App/Logger.hxx"

namespace Tree2Secondaries {

struct Settings {
    void Print() const {
        INFO("SETTINGS");
        INFO("========");
        INFO("IsMC           = %i", isMC);
        INFO("IsSignalMC     = %i", isSignalMC);
        INFO("InputFile      = %s", pathInputFile.c_str());
        INFO("OutputFile     = %s", pathOutputFile.c_str());
        INFO("LimitToNEvents = %lld", limitToNEvents);
    }

    std::string pathInputFile;
    std::string pathOutputFile;
    bool isMC;
    bool isSignalMC;
    long long limitToNEvents;
};

}  // namespace Tree2Secondaries

#endif  // T2S_SETTINGS_HXX
