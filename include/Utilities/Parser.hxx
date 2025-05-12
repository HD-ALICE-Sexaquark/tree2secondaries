#ifndef T2S_PARSER_HXX
#define T2S_PARSER_HXX

#include <cctype>
#include <set>
#include <string>

#include "CLI11.hpp"

#include "Analysis/Settings.hxx"
#include "Math/Constants.hxx"

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
        // input/output paths //
        CLI_APP.add_option("-i,--input", settings.PathInputFiles, "Path of input files")->required();
        CLI_APP.add_option("-o,--output", settings.PathOutputFile, "Path of output file")->expected(1);
        // data type flags //
        CLI_APP.add_flag("-m,--mc", settings.IsMC, "Process MC");
        CLI_APP
            .add_flag("-s,--signal", settings.IsSignalMC, "Process Signal MC")  //
            ->needs("-m");
        //
        CLI_APP
            .add_option("-n,--nevents", settings.LimitToNEvents, "Limit to N events")  //
            ->check(CLI::NonNegativeNumber);
        // reaction channel //
        const std::set<char> allowed_channels{'A', 'D', 'E', 'H'};
        CLI_APP
            .add_option_function<char>(
                "-c,--channel", [&settings](const char& val) { settings.Channel = static_cast<ReactionChannel>(val); },
                "Process a standard reaction channel")
            ->expected(1)
            ->check(CLI::IsMember(allowed_channels));
        CLI_APP
            .add_option_function<char>(
                "-a,--anti", [&settings](const char& val) { settings.Channel = static_cast<ReactionChannel>(std::tolower(val)); },
                "Process an anti reaction channel")
            ->expected(1)
            ->excludes("-c")
            ->check(CLI::IsMember(allowed_channels));
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
                    settings.PathOutputFile = "Packed_Channel" + settings.StrReactionChannel() + "_BkgMC.root";
                } else {
                    settings.PathOutputFile = "Packed_Channel" + settings.StrReactionChannel() + "_SignalMC.root";
                }
            } else {
                settings.PathOutputFile = "Packed_Channel" + settings.StrReactionChannel() + "_Data.root";
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
