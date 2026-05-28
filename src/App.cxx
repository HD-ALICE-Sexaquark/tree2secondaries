#include <ROOT/RNTupleReader.hxx>

#include "App/Parser.hxx"
#include "App/Settings.hxx"
#include "Finder/Finder.hxx"
#include "Packager/Packager.hxx"
#include "Verifier/Verifier.hxx"

int main(int argc, char *argv[]) {

    R2DS::Settings settings;
    R2DS::Parser parser("RNTuple2DoubleStrangeness");
    parser.Parse(argc, argv);
    if (parser.HelpOrError) return parser.ExitCode;
    parser.Assign(settings);
    settings.Print();

    if (settings.Mode == R2DS::EProgramMode::PACKAGER) {

        R2DS::Packager pkgr(settings);
        if (!pkgr.Initialize()) return 1;

        for (auto id_event : *pkgr.fReader) {
            pkgr.fReader->LoadEntry(id_event);
            pkgr.ProcessEvent();
            if (settings.IsMC) pkgr.ProcessInjected();
            pkgr.ProcessTracks();
            pkgr.Pack();
            pkgr.EndOfEvent();
        }
        pkgr.EndOfAnalysis();
    } else if (settings.Mode == R2DS::EProgramMode::FINDER) {

        R2DS::Finder finder(settings);
        if (!finder.Initialize()) return 1;

        for (auto id_event : *finder.fReader) {
            finder.fReader->LoadEntry(id_event);
            finder.ProcessEvent();
            if (settings.IsMC) finder.ProcessInjected();
            finder.Find();
        }
        finder.EndOfAnalysis();

    } else if (settings.Mode == R2DS::EProgramMode::VERIFIER) {

        R2DS::Verifier verifier(settings);
        if (!verifier.Initialize()) return 1;

        for (auto id_event : *verifier.fReader) {
            verifier.fReader->LoadEntry(id_event);
            verifier.ProcessEvent();
            verifier.Verify();
        }
        verifier.EndOfAnalysis();
    } else {
        Logger::Error("main", "Invalid mode.");
        return 1;
    }

    return 0;
}
