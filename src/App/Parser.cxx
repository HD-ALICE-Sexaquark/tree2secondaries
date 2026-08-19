#include <format>
#include <string>
#include <vector>

#include <CLI/CLI.hpp>

#include "App/Parser.hxx"
#include "App/Settings.hxx"

namespace T2DS {

int Parser::Parse(int argc, char* argv[]) {
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

void Parser::Assign(Settings& settings) {
    // -- mode
    CLI::App* mode_cmd = nullptr;
    if (CLI_APP.got_subcommand("find")) {
        settings.Mode = EProgramMode::FINDER;
        mode_cmd = CLI_APP.get_subcommand("find");
    } else if (CLI_APP.got_subcommand("verify")) {
        settings.Mode = EProgramMode::VERIFIER;
        mode_cmd = CLI_APP.get_subcommand("verify");
    }

    // -- data kind
    settings.IsMC = mode_cmd->got_subcommand("mc");
    auto* data_kind_cmd = settings.IsMC ? mode_cmd->get_subcommand("mc") : mode_cmd->get_subcommand("data");

    // -- input paths
    settings.PathInputFiles = InputFiles;

    // -- n events limit
    auto* opt_n = CLI_APP.get_option("-n");
    if (opt_n->count() > 0) {
        settings.LimitToNEvents = opt_n->as<long long>();
    }

    // -- injected mass
    if (settings.Mode == EProgramMode::FINDER && settings.IsMC) {
        settings.SexaquarkMass = data_kind_cmd->get_option("-m")->as<double>();
    }

    // -- output path
    settings.PathOutputFile = CLI_APP.get_option("-o")->as<std::string>();
    if (settings.PathOutputFile.empty()) {
        if (settings.Mode == EProgramMode::FINDER) {
            if (settings.IsMC) {
                settings.PathOutputFile = std::format("FoundRNT_{:.2f}.root", settings.SexaquarkMass);
            } else {
                settings.PathOutputFile = "FoundRNT.root";
            }
        }
        if (settings.Mode == EProgramMode::VERIFIER) {
            settings.PathOutputFile = "VerifiedRNT.root";
        }
    }
}

void Parser::AddOptions() {

    CLI_APP.add_option("-i,--input", InputFiles, "Path(s) of input file(s)")->required();
    CLI_APP.add_option("-o,--output", "Path of output file")->expected(1);
    CLI_APP.add_option("-n,--nevents", "Limit to N events")->expected(1)->check(CLI::PositiveNumber);

    auto add_mass_opt = [](CLI::App* subcmd) {
        subcmd
            ->add_option("-m,--mass", "Assign injected sexaquark mass")  //
            ->expected(1)
            ->check(CLI::IsMember(AllowedMasses))
            ->required();
    };

    // find mode
    auto* find_cmd = CLI_APP.add_subcommand("find", "Read \"AnalysisResults.root\" to search for anti-sexaquark reactions");
    // -- mc
    auto* find_mc_cmd = find_cmd->add_subcommand("mc", "Process MC");
    add_mass_opt(find_mc_cmd);
    // -- rd
    find_cmd->add_subcommand("data", "Process data");
    find_cmd->require_subcommand(1);

    // verify mode
    auto* verify_cmd = CLI_APP.add_subcommand("verify", "Read \"AnalysisResults.root\" to verify the existence of h-dibaryons");
    verify_cmd->add_subcommand("mc", "Process MC");
    verify_cmd->add_subcommand("data", "Process data");
    verify_cmd->require_subcommand(1);

    CLI_APP.require_subcommand(1);
}

}  // namespace T2DS
