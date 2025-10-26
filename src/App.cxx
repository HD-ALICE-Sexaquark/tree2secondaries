#include "App/Parser.hxx"
#include "App/Settings.hxx"
#include "Finder/Finder.hxx"
#include "Math/Constants.hxx"
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

        for (int i_event{0}; i_event < pkgr.NumberEventsToRead(); ++i_event) {
            pkgr.GetEvent(i_event);

            pkgr.ProcessEvent();
            if (pkgr.IsMC()) pkgr.ProcessInjected();
            pkgr.ProcessTracks();

            switch (pkgr.GetReactionChannel()) {
                // standard channels //
                case T2S::EReactionChannel::A:
                    pkgr.FindV0s(T2S::EParticle::AntiLambda);
                    pkgr.FindV0s(T2S::EParticle::Lambda);
                    pkgr.FindV0s(T2S::EParticle::KaonZeroShort);
                    break;
                case T2S::EReactionChannel::D:
                    pkgr.FindV0s(T2S::EParticle::AntiLambda);
                    pkgr.FindV0s(T2S::EParticle::Lambda);
                    pkgr.PackTracks(T2S::EParticle::NegKaon);
                    pkgr.PackTracks(T2S::EParticle::PosKaon);
                    break;
                case T2S::EReactionChannel::E:
                    pkgr.FindV0s(T2S::EParticle::AntiLambda);
                    pkgr.FindV0s(T2S::EParticle::Lambda);
                    pkgr.PackTracks(T2S::EParticle::NegKaon);
                    pkgr.PackTracks(T2S::EParticle::PosKaon);
                    pkgr.PackTracks(T2S::EParticle::PiMinus);
                    pkgr.PackTracks(T2S::EParticle::PiPlus);
                    break;
                case T2S::EReactionChannel::H:
                    pkgr.PackTracks(T2S::EParticle::NegKaon);
                    pkgr.PackTracks(T2S::EParticle::PosKaon);
                    break;
            }
            pkgr.EndOfEvent();
        }
        pkgr.EndOfAnalysis();

    } else {

        T2S::Finder finder(settings);
        if (!finder.Initialize()) return 1;

        for (int i_event{0}; i_event < finder.NumberEventsToRead(); ++i_event) {
            finder.GetEvent(i_event);
            if (finder.IsMC()) finder.Injected_FlattenAndStore();
            finder.Find(finder.GetReactionChannel());
        }
        finder.EndOfAnalysis();
    }

    return 0;
}
