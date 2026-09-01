#include "T2DS/Finder.hxx"
#include "T2DS/Parser.hxx"
#include "T2DS/RunOverInputs.hxx"
#include "T2DS/Settings.hxx"
#include "T2DS/Verifier.hxx"

int main(int argc, char *argv[]) {

    T2DS::Settings settings;
    T2DS::Parser parser("Tree2DoubleStrangeness");
    parser.Parse(argc, argv);
    if (parser.HelpOrError) return parser.ExitCode;
    parser.Assign(settings);
    settings.Print();

    switch (settings.Mode) {
        case (T2DS::EProgramMode::FINDER): {
            T2DS::Finder fndr(settings);
            T2DS::RunOverInputs(fndr, settings, [&] {
                fndr.ProcessEvent();
                if (settings.IsMC) fndr.ProcessInjected();
                fndr.ProcessTracks();
                fndr.FindV0s();
                fndr.FindSexaquarks();
                fndr.EndOfEvent();
            });
            if (!fndr.EndOfAnalysis()) return 1;
            break;
        }
        case (T2DS::EProgramMode::VERIFIER): {
            T2DS::Verifier vrfr(settings);
            T2DS::RunOverInputs(vrfr, settings, [&] {
                vrfr.ProcessEvent();
                if (settings.IsMC) vrfr.ProcessInjected();
                vrfr.ProcessPreFoundLambda();
                vrfr.Verify();
                vrfr.EndOfEvent();
            });
            if (!vrfr.EndOfAnalysis()) return 1;
            break;
        }
    }

    return 0;
}
