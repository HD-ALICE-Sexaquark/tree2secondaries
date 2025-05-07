#include "App/Parser.hxx"
#include "App/Logger.hxx"

#include <cstdlib>

namespace Tree2Secondaries {

Parser::Parser(int t_argc, char* t_argv[]) : exitCode(0), argc(t_argc), argv(t_argv) {}

void Parser::PrintHelp() {
    INFO("Usage: ./App [OPTIONS]");
    INFO("Options:");
    INFO("  -h,--help                   Print this help message and exit");
    INFO("  -i,--input REQUIRED         Path of input files");
    INFO("  -o,--output REQUIRED        Path of output file");
    INFO("  -m,--mc                     Flag to process MC");
    INFO("  -s,--signal NEEDS: --mc     Flag to process Signal MC");
    INFO("  -n,--nevents                Limit to N events");
}

// return true if parsing was signal is help or error
bool Parser::Parse(Settings& Settings) {

    if (argc < 2) {
        ERROR("No arguments provided");
        PrintHelp();
        exitCode = 1;
        return false;
    }

    for (int i{1}; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            PrintHelp();
            exitCode = 0;
            return false;
        } else if (arg == "-i" || arg == "--input") {
            if (i + 1 < argc) {
                options["input"] = argv[i + 1];
            }
        } else if (arg == "-o" || arg == "--output") {
            if (i + 1 < argc) {
                options["output"] = argv[i + 1];
            }
        } else if (arg == "-m" || arg == "--mc") {
            options["mc"] = "true";
        } else if (arg == "-s" || arg == "--signal") {
            options["signal"] = "true";
        } else if (arg == "-n" || arg == "--nevents") {
            if (i + 1 < argc) {
                options["nevents"] = argv[i + 1];
            }
        }
    }

    return Validate(Settings);
}

bool Parser::Validate(Settings& Settings) {

    if (options.count("signal") && !options.count("mc")) {
        ERROR("--signal requires --mc");
        exitCode = 1;
        return false;
    }
    if (!options.count("input")) {
        ERROR("Missing value for --input");
        exitCode = 1;
        return false;
    }
    if (!options.count("output")) {
        ERROR("Missing value for --output");
        exitCode = 1;
        return false;
    }

    Settings.isMC = options["mc"] == "true";
    Settings.isSignalMC = options["signal"] == "true";
    Settings.pathInputFile = options["input"];
    Settings.pathOutputFile = options["output"];
    Settings.limitToNEvents = options["nevents"].empty() ? 0 : std::atoi(options["nevents"].c_str());

    return true;
}

}  // namespace Tree2Secondaries
