#pragma once

#include <string>

#include <CLI/CLI.hpp>

namespace Skimmer {

// # Command-line Interface Settings # //

struct Settings {
    std::string ConfigPath;
};

// # Parser # //

class Parser : CLI::App {
   public:
    explicit Parser(const std::string& app_name) : CLI::App{app_name} {}

    int Parse(int argc, char* argv[], Settings& cli) {
        argv = ensure_utf8(argv);

        add_option("config", cli.ConfigPath, "Path of JSON config file")  //
            ->required()                                                  //
            ->check(CLI::ExistingFile);

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
