#include "App/Parser.hxx"
#include "App/Settings.hxx"
#include "Finder/Finder.hxx"
#include "Packager/Packager.hxx"
#include "Verifier/Verifier.hxx"

int main(int argc, char *argv[]) {

    T2DS::Settings settings;
    T2DS::Parser parser("Tree2DoubleStrangeness");
    parser.Parse(argc, argv);
    if (parser.HelpOrError) return parser.ExitCode;
    parser.Assign(settings);
    settings.Print();

    if (settings.Mode == T2DS::EProgramMode::PACKAGER) {

        T2DS::Packager pkgr(settings);
        const auto n_events = pkgr.NumberEventsToRead();
        for (long long id_event = 0; id_event < n_events; ++id_event) {
            pkgr.Load(id_event);
            pkgr.ProcessEvent();
            if (settings.IsMC) pkgr.ProcessInjectedSexa();
            pkgr.ProcessTracks();
            pkgr.Pack();
            pkgr.EndOfEvent();
        }
        pkgr.EndOfAnalysis();

    } else if (settings.Mode == T2DS::EProgramMode::FINDER) {

        T2DS::Finder fndr(settings);
        const auto n_events = fndr.NumberEventsToRead();
        for (unsigned long id_event = 0; id_event < n_events; ++id_event) {
            fndr.Load(id_event);
            fndr.ProcessEvent();
            if (settings.IsMC) fndr.ProcessInjected();
            fndr.Find();
            fndr.EndOfEvent();
        }
        fndr.EndOfAnalysis();

    } else if (settings.Mode == T2DS::EProgramMode::VERIFIER) {

        T2DS::Verifier vrfr(settings);
        const auto n_events = vrfr.NumberEventsToRead();
        for (long long id_event = 0; id_event < n_events; ++id_event) {
            vrfr.Load(id_event);
            vrfr.ProcessEvent();
            if (settings.IsMC) vrfr.ProcessInjected();
            vrfr.ProcessPreFoundLambda();
            vrfr.Verify();
            vrfr.EndOfEvent();
        }
        vrfr.EndOfAnalysis();

    } else {

        Logger::Error("main", "Invalid mode.");
        return 1;
    }

    return 0;
}
