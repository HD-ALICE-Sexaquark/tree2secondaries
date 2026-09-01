#pragma once

#include <cstdint>
#include <string>

#include <CLI/CLI.hpp>

namespace Skimmer {

// # Command-line Interface Settings # //

struct Settings {
    std::string ConfigPath;
    std::string OutputDir{"cache"};
    std::uint64_t NEvents{0};
};

// # Parser # //

class Parser : CLI::App {
   public:
    explicit Parser(const std::string& app_name) : CLI::App{app_name} {}

    int Parse(int argc, char* argv[], Settings& cli) {
        argv = ensure_utf8(argv);

        add_option("-c,--config", cli.ConfigPath, "Path of JSON config file")  //
            ->required()                                                       //
            ->check(CLI::ExistingFile);
        add_option("-o,--output", cli.OutputDir, "Output directory, or a path ending in .root")  //
            ->expected(1);
        add_option("-n,--nevents", cli.NEvents, "Limit to N events per sample (0 = all)")  //
            ->expected(1);

        try {
            parse(argc, argv);
        } catch (const CLI::ParseError& e) {
            ExitCode = e.get_exit_code();
            HelpOrError = (e.get_name() == "CallForHelp") || (ExitCode != 0);
            return exit(e);
        }

        return 0;
    }

    int ExitCode{0};
    bool HelpOrError{false};
};

}  // namespace Skimmer
