#include <cstdlib>

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

        for (size_t i_event = 0; i_event < pkgr.NumberEventsToRead(); ++i_event) {
            pkgr.GetEvent(i_event);

            pkgr.ProcessEvent();
            if (pkgr.IsMC()) pkgr.ProcessInjected();
            pkgr.ProcessTracks();

            switch (pkgr.GetReactionChannel()) {
                // standard channels //
                case T2S::EReactionChannel::A:
                    pkgr.FindV0s(T2S::PID_V0::AntiLambda);
                    pkgr.FindV0s(T2S::PID_V0::Lambda);
                    pkgr.FindV0s(T2S::PID_V0::KaonZeroShort);
                    break;
                case T2S::EReactionChannel::D:
                    pkgr.FindV0s(T2S::PID_V0::AntiLambda);
                    pkgr.FindV0s(T2S::PID_V0::Lambda);
                    pkgr.PackTracks(T2S::PID_StableParticle::NegKaon);
                    pkgr.PackTracks(T2S::PID_StableParticle::PosKaon);
                    break;
                case T2S::EReactionChannel::E:
                    pkgr.FindV0s(T2S::PID_V0::AntiLambda);
                    pkgr.FindV0s(T2S::PID_V0::Lambda);
                    pkgr.PackTracks(T2S::PID_StableParticle::NegKaon);
                    pkgr.PackTracks(T2S::PID_StableParticle::PosKaon);
                    pkgr.PackTracks(T2S::PID_StableParticle::PiMinus);
                    pkgr.PackTracks(T2S::PID_StableParticle::PiPlus);
                    break;
                case T2S::EReactionChannel::H:
                    pkgr.PackTracks(T2S::PID_StableParticle::NegKaon);
                    pkgr.PackTracks(T2S::PID_StableParticle::PosKaon);
                    break;
                default:
                    break;
            }
            pkgr.EndOfEvent();
        }
        pkgr.EndOfAnalysis();

    } else {

        T2S::Finder finder(settings);
        if (!finder.Initialize()) return 1;

        for (size_t i_event = 0; i_event < finder.NumberEventsToRead(); ++i_event) {
            finder.GetEvent(i_event);

            finder.ProcessEvent();
            if (finder.IsMC()) finder.ProcessInjected();
            finder.Find(finder.GetReactionChannel());
        }
        finder.EndOfAnalysis();
    }

    return 0;
}
