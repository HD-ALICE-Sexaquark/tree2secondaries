#ifndef T2S_PARSER_HXX
#define T2S_PARSER_HXX

#include <cctype>
#include <set>
#include <string>
#include <vector>

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
        const std::set<float> allowed_masses{1.73, 1.8, 1.87, 1.94, 2.01};

        auto add_channels_opt = [allowed_channels](CLI::App* subcmd) {
            subcmd->add_option("-c,--channel", "Process a standard reaction channel")
                ->expected(1)
                ->check(CLI::IsMember(allowed_channels))
                ->required();
        };

        auto add_mass_opt = [allowed_masses](CLI::App* subcmd) {
            subcmd
                ->add_option("-m,--mass", "Assign Injected Sexaquark Mass")  //
                ->expected(1)
                ->check(CLI::IsMember(allowed_masses))
                ->required();
        };

        CLI_APP
            .add_option("-i,--input", InputFiles, "Path(s) of input file(s)")  //
            ->required();
        CLI_APP
            .add_option("-o,--output", "Path of output file")  //
            ->expected(1);
        CLI_APP
            .add_option("-n,--nevents", "Limit to N events")  //
            ->expected(1)
            ->check(CLI::NonNegativeNumber);

        // # package mode //

        auto* package_cmd = CLI_APP.add_subcommand("pack", "Package V0s and necessary tracks");
        auto* package_mc_cmd = package_cmd->add_subcommand("mc", "Process MC");
        add_channels_opt(package_mc_cmd);
        add_mass_opt(package_mc_cmd);
        auto* package_data = package_cmd->add_subcommand("data", "Process data");
        add_channels_opt(package_data);
        package_cmd->require_subcommand(1);

        // # search mode //

        auto* search_cmd = CLI_APP.add_subcommand("search", "Search for anti-sexaquark reactions");
        auto* search_mc_cmd = search_cmd->add_subcommand("mc", "Process MC");
        add_channels_opt(search_mc_cmd);
        add_mass_opt(search_mc_cmd);
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
        // reaction channels //
        auto* opt_channel = type_cmd->get_option("-c");
        settings.ReactionChannel = static_cast<EReactionChannel>(opt_channel->as<char>());
        // sexaquark mass //
        if (settings.IsMC) settings.SexaquarkMass = type_cmd->get_option("-m")->as<double>();
        // input path //
        settings.PathInputFiles = InputFiles;
        // n events limit //
        settings.LimitToNEvents = CLI_APP.get_option("-n")->as<int>();
        // output path //
        settings.PathOutputFile = CLI_APP.get_option("-o")->as<std::string>();
        if (settings.PathOutputFile.empty()) {
            std::string filename_prefix{settings.DoTheSearch ? "Found" : "Packed"};
            std::string filename_suffix{};
            if (settings.IsMC)
                filename_suffix = std::format("MC_{}{:.2f}", static_cast<char>(settings.ReactionChannel), settings.SexaquarkMass);
            else
                filename_suffix = std::format("Data_Channel{}", static_cast<char>(settings.ReactionChannel));
            settings.PathOutputFile = std::format("{}_{}.root", filename_prefix, filename_suffix);
        }
    }

    // member vars //
    int ExitCode{0};
    bool HelpOrError{false};

   protected:
    CLI::App CLI_APP;
    std::vector<std::string> InputFiles;
};

}  // namespace Tree2Secondaries

#endif  // T2S_PARSER_HXX
