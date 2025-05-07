#ifndef PARSER_HXX
#define PARSER_HXX

#include <string>
#include <unordered_map>

#include "App/Settings.hxx"

namespace Tree2Secondaries {

class Parser {
   public:
    Parser(int t_argc, char* t_argv[]);
    bool Parse(Settings& options);
    static void PrintHelp();
    int exitCode;

   private:
    int argc;
    char** argv;
    std::unordered_map<std::string, std::string> options;

    bool Validate(Settings& Settings);
};

}  // namespace Tree2Secondaries

#endif  // PARSER_HXX
