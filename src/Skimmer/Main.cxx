#include <exception>

#include "Skimmer/Config.hxx"
#include "Skimmer/Parser.hxx"
#include "Utils/Logger.hxx"

#include "Skimmer/Skimmer.hxx"

int main(int argc, char* argv[]) {

    Skimmer::Settings cli;
    Skimmer::Parser parser("Skimmer");
    parser.Parse(argc, argv, cli);
    if (parser.HelpOrError) return parser.ExitCode;

    try {
        const Skimmer::Config config = Skimmer::Load(cli.ConfigPath);
        config.Print();
        Skimmer::Run(config, cli.OutputDir, cli.NEvents);
    } catch (const std::exception& exc) {
        Logger::Error("Main", "{}", exc.what());
        return 1;
    }

    return 0;
}
