#pragma once

#include <array>
#include <string>
#include <vector>

#include <CLI/CLI.hpp>

#include "T2DS/Settings.hxx"

namespace T2DS {

class Parser {
   public:
    Parser(const Parser&) = delete;
    Parser(Parser&&) = delete;
    Parser& operator=(const Parser&) = delete;
    Parser& operator=(Parser&&) = delete;

    explicit Parser(const std::string& app_name) : CLI_APP{app_name} { AddOptions(); }
    ~Parser() = default;

    int Parse(int argc, char* argv[]);
    void Assign(Settings& settings);

    int ExitCode{0};
    bool HelpOrError{false};

   protected:
    void AddOptions();

    static constexpr std::array<double, 5> AllowedMasses{1.73, 1.8, 1.87, 1.94, 2.01};

    CLI::App CLI_APP;
    std::vector<std::string> InputFiles;
};

}  // namespace T2DS
