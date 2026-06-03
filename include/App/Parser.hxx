#pragma once

#include <array>
#include <string>

#include <CLI/CLI.hpp>

#include "App/Settings.hxx"

namespace R2DS {

class Parser {
   public:
    Parser(const Parser&) = delete;
    Parser(Parser&&) = delete;
    Parser& operator=(const Parser&) = delete;
    Parser& operator=(Parser&&) = delete;

    explicit Parser(const std::string& app_name) : CLI_APP{app_name} { AddOptions(); }
    ~Parser() = default;

    int Parse(int argc, char* argv[]) {
        argv = CLI_APP.ensure_utf8(argv);
        try {
            CLI_APP.parse(argc, argv);
        } catch (const CLI::ParseError& e) {
            ExitCode = e.get_exit_code();
            HelpOrError = (e.get_name() == "CallForHelp") || (ExitCode != 0);
            return CLI_APP.exit(e);
        }
        return 0;
    }

    void Assign(Settings& settings) {
        // -- mode
        CLI::App* mode_cmd = nullptr;
        if (CLI_APP.got_subcommand("pack")) {
            settings.Mode = EProgramMode::PACKAGER;
            mode_cmd = CLI_APP.get_subcommand("pack");
        } else if (CLI_APP.got_subcommand("search")) {
            settings.Mode = EProgramMode::FINDER;
            mode_cmd = CLI_APP.get_subcommand("search");
        } else if (CLI_APP.got_subcommand("verify")) {
            settings.Mode = EProgramMode::VERIFIER;
            mode_cmd = CLI_APP.get_subcommand("verify");
        }

        // -- data kind
        settings.IsMC = mode_cmd->got_subcommand("mc");
        auto* data_kind_cmd = settings.IsMC ? mode_cmd->get_subcommand("mc") : mode_cmd->get_subcommand("data");

        // -- input path
        settings.PathInputFile = InputFile;

        // -- n events limit
        auto* opt_n = CLI_APP.get_option("-n");
        if (opt_n->count() > 0) {
            settings.LimitToNEvents = opt_n->as<long long>();
        }

        // -- output path, reaction channel & mass
        settings.PathOutputFile = CLI_APP.get_option("-o")->as<std::string>();
        if (settings.PathOutputFile.empty()) {
            if (settings.Mode == EProgramMode::PACKAGER) {
                if (settings.IsMC) {
                    settings.ReactionChannel = data_kind_cmd->get_option("-c")->as<char>();
                    settings.SexaquarkMass = data_kind_cmd->get_option("-m")->as<double>();
                    settings.PathOutputFile = std::format("PackedRNT_{}{:.2f}.root", settings.ReactionChannel, settings.SexaquarkMass);
                } else {
                    settings.PathOutputFile = "PackedRNT.root";
                }
            }
            if (settings.Mode == EProgramMode::FINDER) {
                settings.ReactionChannel = data_kind_cmd->get_option("-c")->as<char>();
                settings.PathOutputFile = std::format("FoundRNT_Channel{}.root", settings.ReactionChannel);
            }
            if (settings.Mode == EProgramMode::VERIFIER) {
                settings.PathOutputFile = "VerifiedRNT.root";
            }
        }
    }

    int ExitCode{0};
    bool HelpOrError{false};

   protected:
    void AddOptions() {

        CLI_APP.add_option("-i,--input", InputFile, "Path of input file")->required();
        CLI_APP.add_option("-o,--output", "Path of output file")->expected(1);
        CLI_APP.add_option("-n,--nevents", "Limit to N events")->expected(1)->check(CLI::PositiveNumber);

        constexpr std::array<char, 3> allowed_channels{'A', 'D', 'H'};
        constexpr std::array<double, 5> allowed_masses{1.73, 1.8, 1.87, 1.94, 2.01};

        auto add_channels_opt = [allowed_channels](CLI::App* subcmd) {
            subcmd->add_option("-c,--channel", "Process a standard reaction channel")
                ->expected(1)
                ->check(CLI::IsMember(allowed_channels))
                ->required();
        };
        auto add_mass_opt = [allowed_masses](CLI::App* subcmd) {
            subcmd
                ->add_option("-m,--mass", "Assign injected sexaquark mass")  //
                ->expected(1)
                ->check(CLI::IsMember(allowed_masses))
                ->required();
        };

        // package mode
        auto* package_cmd = CLI_APP.add_subcommand("pack", "Read \"EventsRNT.root\" to package injected anti-sexaquarks, tracks and V0s");
        // -- mc
        auto* package_mc_cmd = package_cmd->add_subcommand("mc", "Process MC");
        add_channels_opt(package_mc_cmd);
        add_mass_opt(package_mc_cmd);
        // -- rd
        package_cmd->add_subcommand("data", "Process data");
        package_cmd->require_subcommand(1);

        // -- search mode
        auto* search_cmd = CLI_APP.add_subcommand("search", "Read \"PackedRNT.root\" files to search for anti-sexaquark reactions");
        // -- mc
        auto* search_mc_cmd = search_cmd->add_subcommand("mc", "Process MC");
        add_channels_opt(search_mc_cmd);
        // -- rd
        auto* search_data_cmd = search_cmd->add_subcommand("data", "Process data");
        add_channels_opt(search_data_cmd);
        search_cmd->require_subcommand(1);

        // -- verify mode
        auto* verify_cmd = CLI_APP.add_subcommand("verify", "Read \"EventsRNT.root\" to verify the existence of h-dibaryons");
        verify_cmd->add_subcommand("mc", "Process MC");
        verify_cmd->add_subcommand("data", "Process data");
        verify_cmd->require_subcommand(1);

        CLI_APP.require_subcommand(1);
    }

    CLI::App CLI_APP;
    std::string InputFile;
};

}  // namespace R2DS
