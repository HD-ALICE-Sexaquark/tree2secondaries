#ifndef T2S_PARSER_HXX
#define T2S_PARSER_HXX

#include <string>

#include "CLI11.hpp"

#include "Analysis/Settings.hxx"

namespace Tree2Secondaries {

class Parser {

   public:
    Parser(int argc, char* argv[], Settings& settings, std::string app_description)
        : ExitCode(0), HelpOrError(false), CLI_APP(std::move(app_description)) {
        AddOptions(settings);
        Parse(argc, argv, settings);
    }
    ~Parser() = default;

    void AddOptions(Settings& settings) {
        CLI_APP.add_option("-i,--input", settings.PathInputFiles, "Path of input files");
        CLI_APP.add_option("-o,--output", settings.PathOutputFile, "Path of output file");
        CLI_APP.add_flag("-m,--mc", settings.IsMC, "Flag to process MC");
        CLI_APP
            .add_flag("-s,--signal", settings.IsSignalMC, "Flag to process Signal MC")  //
            ->needs("-m");
        CLI_APP
            .add_option("-n,--nevents", settings.LimitToNEvents, "Limit to N events")  //
            ->check(CLI::NonNegativeNumber & CLI::TypeValidator<int>());
    }

    int Parse(int argc, char* argv[], Settings& settings) {
        argv = CLI_APP.ensure_utf8(argv);
        try {
            CLI_APP.parse(argc, argv);
        } catch (const CLI::ParseError& e) {
            ExitCode = e.get_exit_code();
            HelpOrError = e.get_name() == "CallForHelp" || (ExitCode != 0);
            return CLI_APP.exit(e);
        }
        // default options //
        if (settings.PathOutputFile.empty()) {
            if (settings.IsMC) {
                if (settings.IsSignalMC) {
                    settings.PathOutputFile = "../files/AnalysisResults.root";
                } else {
                    settings.PathOutputFile = "../files/AnalysisResults_what.root";
                }
            } else {
                settings.PathOutputFile = "../files/AnalysisResults_data_15o_245554_001.root";
            }
        }
        if (settings.PathOutputFile.empty()) {
            if (settings.IsMC) {
                if (settings.IsSignalMC) {
                    settings.PathOutputFile = "../files/SecondariesResults_BkgMC.root";
                } else {
                    settings.PathOutputFile = "../files/SecondariesResults_SignalMC.root";
                }
            } else {
                settings.PathOutputFile = "../files/SecondariesResults_Data.root";
            }
        }
        return 0;
    }
    int ExitCode;
    bool HelpOrError;

   protected:
    void AddOptions();
    CLI::App CLI_APP;
};

}  // namespace Tree2Secondaries

#endif  // T2S_PARSER_HXX
