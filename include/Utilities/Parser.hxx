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
    Parser(int argc, char* argv[], Settings& settings, std::string app_description) : CLI_APP{std::move(app_description)} {
        AddOptions(settings);
        Parse(argc, argv, settings);
    }
    ~Parser() = default;

    void AddOptions(Settings& settings) {
        // input/output paths //
        CLI_APP
            .add_option("-i,--input", settings.PathInputFiles, "Path(s) of input file(s)")  //
            ->required()
            ->expected(1, -1);
        CLI_APP
            .add_option("-o,--output", settings.PathOutputFile, "Path of output file")  //
            ->expected(1);
        // n events //
        CLI_APP
            .add_option("-n,--nevents", settings.LimitToNEvents, "Limit to N events")  //
            ->expected(1)
            ->check(CLI::NonNegativeNumber);

        // MC //
        auto* mc_cmd = CLI_APP.add_subcommand("mc", "Process MC");

        // is signal MC //
        mc_cmd->add_flag("-s,--signal", settings.IsSignalMC, "Process Signal MC");

        // reaction channel //
        auto* ch_opt_group = mc_cmd->add_option_group("channels");
        const std::set<char> allowed_channels{'A', 'D', 'E', 'H'};
        ch_opt_group
            ->add_option_function<char>(
                "-c,--channel", [&settings](char val) { settings.Channel = static_cast<ReactionChannel>(val); },
                "Process a standard reaction channel")
            ->expected(1)
            ->check(CLI::IsMember(allowed_channels));
        ch_opt_group
            ->add_option_function<char>(
                "-a,--anti", [&settings](char val) { settings.Channel = static_cast<ReactionChannel>(std::tolower(val)); },
                "Process an anti-reaction channel")
            ->expected(1)
            ->check(CLI::IsMember(allowed_channels));
        ch_opt_group->require_option(1);

        // data //
        CLI_APP.add_subcommand("data", "Process data");

        CLI_APP.require_subcommand(1);
    }

    int Parse(int argc, char* argv[], Settings& settings) {
        argv = CLI_APP.ensure_utf8(argv);
        try {
            CLI_APP.parse(argc, argv);
            settings.IsMC = CLI_APP.got_subcommand("mc");
        } catch (const CLI::ParseError& e) {
            ExitCode = e.get_exit_code();
            HelpOrError = (e.get_name() == "CallForHelp") || (ExitCode != 0);
            return CLI_APP.exit(e);
        }
        // default options //
        if (settings.PathOutputFile.empty()) {
            if (settings.IsMC) {
                if (settings.IsSignalMC) {
                    settings.PathOutputFile = "Packed_Channel" + settings.StrReactionChannel() + "_SignalMC.root";
                } else {
                    settings.PathOutputFile = "Packed_Channel" + settings.StrReactionChannel() + "_BkgMC.root";
                }
            } else {
                settings.PathOutputFile = "Packed_Data.root";
            }
        }
        return 0;
    }
    int ExitCode{0};
    bool HelpOrError{false};

   protected:
    void AddOptions();
    CLI::App CLI_APP;
};

}  // namespace Tree2Secondaries

#endif  // T2S_PARSER_HXX
