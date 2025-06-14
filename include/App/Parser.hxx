#ifndef T2S_PARSER_HXX
#define T2S_PARSER_HXX

#include <cctype>
#include <set>
#include <string>

#include <CLI/CLI.hpp>

#include "App/Settings.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries {

class Parser {
   public:
    Parser(const Parser&) = delete;
    Parser(Parser&&) = delete;
    Parser& operator=(const Parser&) = delete;
    Parser& operator=(Parser&&) = delete;

    explicit Parser(const std::string& app_name) : CLI_APP{app_name} { AddOptions(); }
    ~Parser() = default;

    void AddOptions() {

        const std::set<char> allowed_channels{'A', 'D', 'E', 'H'};

        auto add_channels_opt = [allowed_channels](CLI::App* subcmd) {
            auto* ch_opt_group = subcmd->add_option_group("channels");
            ch_opt_group->add_option("-c,--channel", "Process a standard reaction channel")->expected(1)->check(CLI::IsMember(allowed_channels));
            ch_opt_group->add_option("-a,--anti", "Process an anti-reaction channel")->expected(1)->check(CLI::IsMember(allowed_channels));
            ch_opt_group->require_option(1);
        };

        CLI_APP
            .add_option("-i,--input", "Path(s) of input file(s)")  //
            ->required()
            ->expected(1, -1);
        CLI_APP
            .add_option("-o,--output", "Path of output file")  //
            ->expected(1);
        CLI_APP
            .add_option("-n,--nevents", "Limit to N events")  //
            ->expected(1)
            ->check(CLI::NonNegativeNumber);

        // # package mode //

        auto* package_cmd = CLI_APP.add_subcommand("pack", "Package V0s and necessary tracks");
        auto* package_mc = package_cmd->add_subcommand("mc", "Process MC");
        package_mc->add_flag("-s,--signal", "Process Signal MC");
        add_channels_opt(package_mc);
        package_cmd->add_subcommand("data", "Process data");
        package_cmd->require_subcommand(1);

        // # search mode //

        auto* search_cmd = CLI_APP.add_subcommand("search", "Search for anti-sexaquark reactions");
        auto* search_mc = search_cmd->add_subcommand("mc", "Process MC");
        search_mc->add_flag("-s,--signal", "Process Signal MC");
        add_channels_opt(search_mc);
        auto* search_data = search_cmd->add_subcommand("data", "Process data");
        add_channels_opt(search_data);
        search_cmd->require_subcommand(1);

        CLI_APP.require_subcommand(1);
    }

    int Parse(int argc, char* argv[]) {
        argv = CLI_APP.ensure_utf8(argv);
        try {
            CLI_APP.parse(argc, argv);
        } catch (const CLI::ParseError& e) {
            ExitCode = e.get_exit_code();
            HelpOrError = (e.get_name() == "CallForHelp") || (ExitCode != 0);
            return CLI_APP.exit(e);
        }
        return 0;  // = success
    }

    void Assign(Settings& settings) {
        // mode //
        settings.DoTheSearch = CLI_APP.got_subcommand("search");
        auto* mode_cmd = settings.DoTheSearch ? CLI_APP.get_subcommand("search") : CLI_APP.get_subcommand("pack");
        // mc vs data //
        settings.IsMC = mode_cmd->got_subcommand("mc");
        auto* type_cmd = settings.IsMC ? mode_cmd->get_subcommand("mc") : mode_cmd->get_subcommand("data");
        if (settings.IsMC) settings.IsSignalMC = type_cmd->get_option("-s")->as<bool>();
        // reaction channels //
        if (settings.DoTheSearch || settings.IsMC) {
            auto* opt_grp = type_cmd->get_option_group("channels");
            auto* std_channel = opt_grp->get_option("-c");
            auto* anti_channel = opt_grp->get_option("-a");
            if (std_channel->count()) {
                settings.Channel = static_cast<ReactionChannel>(std_channel->as<char>());
            } else {
                settings.Channel = static_cast<ReactionChannel>(std::tolower(anti_channel->as<char>()));
            }
        }
        // input path //
        settings.PathInputFiles = CLI_APP.get_option("-i")->as<std::vector<std::string>>();
        // n events limit //
        settings.LimitToNEvents = CLI_APP.get_option("-n")->as<long long>();
        // output path //
        if (settings.PathOutputFile.empty()) {
            std::string prefix_file{settings.DoTheSearch ? "Searched" : "Packed"};
            std::string suffix_file{".root"};
            if (settings.IsMC && settings.IsSignalMC)
                suffix_file = "SignalMC" + suffix_file;
            else if (settings.IsMC && !settings.IsSignalMC)
                suffix_file = "BkgMC" + suffix_file;
            else
                suffix_file = "Data" + suffix_file;
            settings.PathOutputFile = prefix_file + "_" + settings.StrReactionChannel() + "_" + suffix_file;
        }
    }

    // member vars //
    int ExitCode{0};
    bool HelpOrError{false};

   protected:
    CLI::App CLI_APP;
};

}  // namespace Tree2Secondaries

#endif  // T2S_PARSER_HXX
