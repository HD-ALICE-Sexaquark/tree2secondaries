#pragma once

#include <array>
#include <string>

#include <CLI/CLI.hpp>

#include "common/DB_ReactionChannels.hpp"

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
        if (CLI_APP.get_subcommand("search")) {
            settings.Mode = EProgramMode::FINDER;
            mode_cmd = CLI_APP.get_subcommand("search");
        } else if (CLI_APP.get_subcommand("pack")) {
            settings.Mode = EProgramMode::PACKAGER;
            mode_cmd = CLI_APP.get_subcommand("pack");
        } else if (CLI_APP.get_subcommand("verify")) {
            settings.Mode = EProgramMode::VERIFIER;
            mode_cmd = CLI_APP.get_subcommand("verify");
        }

        // -- reaction channel & mass (not applicable to verifier)
        if (settings.Mode == EProgramMode::PACKAGER) {
            auto* opt_channel = mode_cmd->get_option("-c");
            settings.ReactionChannel = DB::ReactionChannels::FindReactionChannel(opt_channel->as<char>());
            settings.SexaquarkMass = mode_cmd->get_option("-m")->as<double>();
        }

        // -- input path
        settings.PathInputFile = InputFile;

        // -- n events limit
        auto* opt_n = CLI_APP.get_option("-n");
        if (opt_n->count() > 0) {
            settings.LimitToNEvents = opt_n->as<long long>();
        }

        // -- output path
        settings.PathOutputFile = CLI_APP.get_option("-o")->as<std::string>();
        if (settings.PathOutputFile.empty()) {
            if (settings.Mode == EProgramMode::VERIFIER) {
                if (settings.ReactionChannel.has_value() && settings.SexaquarkMass.has_value()) {
                    settings.PathOutputFile =
                        std::format("Packed_{}{:.2f}.root", settings.ReactionChannel.value().name, settings.SexaquarkMass.value());
                } else {
                    settings.PathOutputFile = "Packed.root";
                }
                if (settings.Mode == EProgramMode::VERIFIER) settings.PathOutputFile = "Verified.root";
                if (settings.Mode == EProgramMode::FINDER) settings.PathOutputFile = "Found.root";
            }
        }
    }

    int ExitCode{0};
    bool HelpOrError{false};

   protected:
    void AddOptions() {

        constexpr std::array<char, 3> allowed_channels{'A', 'D', 'H'};
        constexpr std::array<double, 5> allowed_masses{1.73, 1.8, 1.87, 1.94, 2.01};

        CLI_APP.add_option("-i,--input", InputFile, "Path of input file")->required();
        CLI_APP.add_option("-o,--output", "Path of output file")->expected(1);
        CLI_APP.add_option("-n,--nevents", "Limit to N events")->expected(1)->check(CLI::PositiveNumber);

        // -- package mode
        auto* package_cmd = CLI_APP.add_subcommand("pack", "Read \"AnalysisResults.root\" to package injected anti-sexaquarks, tracks and V0s");
        package_cmd  //
            ->add_option("-c,--channel", "(Required when reading MC) Process a standard reaction channel")
            ->expected(1)
            ->check(CLI::IsMember(allowed_channels));
        package_cmd  //
            ->add_option("-m,--mass", "(Required when reading MC) Assign injected sexaquark mass")
            ->expected(1)
            ->check(CLI::IsMember(allowed_masses));

        // -- search mode
        CLI_APP.add_subcommand("search", "Read \"Packed.root\" files to search for anti-sexaquark reactions");

        // -- verify mode
        CLI_APP.add_subcommand("verify", "Read \"AnalysisResults.root\" to verify the existence of h-dibaryons");

        CLI_APP.require_subcommand(1);
    }

    CLI::App CLI_APP;
    std::string InputFile;
};

}  // namespace R2DS
