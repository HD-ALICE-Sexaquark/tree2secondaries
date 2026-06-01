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
        if (!pkgr.CheckArguments()) return parser.ExitCode;
        for (unsigned long id_event = 0; id_event < pkgr.NumberEventsToRead(); ++id_event) {
            pkgr.Load(id_event);
            pkgr.ProcessEvent();
            pkgr.ProcessInjected();
            pkgr.ProcessTracks();
            pkgr.Pack();
            pkgr.EndOfEvent();
        }
        pkgr.EndOfAnalysis();

    } else if (settings.Mode == R2DS::EProgramMode::FINDER) {

        R2DS::Finder fndr(settings);
        for (unsigned long id_event = 0; id_event < fndr.NumberEventsToRead(); ++id_event) {
            fndr.Load(id_event);
            fndr.ProcessEvent();
            fndr.ProcessInjected();
            fndr.Find();
        }
        fndr.EndOfAnalysis();

    } else if (settings.Mode == R2DS::EProgramMode::VERIFIER) {

        R2DS::Verifier vrfr(settings);
        for (unsigned long id_event = 0; id_event < vrfr.NumberEventsToRead(); ++id_event) {
            vrfr.Load(id_event);
            vrfr.ProcessEvent();
            vrfr.Verify();
        }
        vrfr.EndOfAnalysis();

    } else {

        Logger::Error("main", "Invalid mode.");
        return 1;
    }

    return 0;
}
