#include "App/Parser.hxx"
#include "App/Settings.hxx"
#include "Finder/Finder.hxx"
#include "Packager/Packager.hxx"

namespace T2S = Tree2Secondaries;

int main(int argc, char *argv[]) {

    T2S::Settings settings;
    T2S::Parser parser("Tree2Secondaries");
    parser.Parse(argc, argv);
    if (parser.HelpOrError) return parser.ExitCode;
    parser.Assign(settings);
    settings.Print();

    if (!settings.DoTheSearch) {

        T2S::Packager pkgr(settings);
        if (!pkgr.Initialize()) return 1;

        for (long long i_event = 0; i_event < pkgr.NumberEventsToRead(); ++i_event) {
            pkgr.GetEvent(i_event);
            pkgr.ProcessEvent();
            if (settings.IsMC) pkgr.ProcessInjected();
            pkgr.ProcessTracks();
            pkgr.Pack();
            pkgr.EndOfEvent();
        }
        pkgr.EndOfAnalysis();

    } else {

        T2S::Finder finder(settings);
        if (!finder.Initialize()) return 1;

        for (long long i_event = 0; i_event < finder.NumberEventsToRead(); ++i_event) {
            finder.GetEvent(i_event);
            finder.ProcessEvent();
            if (settings.IsMC) finder.ProcessInjected();
            finder.Find();
        }
        finder.EndOfAnalysis();
    }

    return 0;
}
