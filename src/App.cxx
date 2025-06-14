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

        for (long long i_event{0}; i_event < pkgr.NumberEventsToRead(); ++i_event) {
            pkgr.GetEvent(i_event);

            pkgr.ProcessEvent();
            if (pkgr.IsSignalMC()) pkgr.ProcessInjected();
            pkgr.ProcessTracks();

            switch (pkgr.GetReactionChannel()) {
                // standard channels //
                case T2S::ReactionChannel::A:
                    pkgr.FindV0s(T2S::PdgCode::AntiLambda);
                    pkgr.FindV0s(T2S::PdgCode::KaonZeroShort);
                    break;
                case T2S::ReactionChannel::D:
                    pkgr.FindV0s(T2S::PdgCode::AntiLambda);
                    pkgr.PackTracks(T2S::PdgCode::PosKaon);
                    break;
                case T2S::ReactionChannel::E:
                    pkgr.FindV0s(T2S::PdgCode::AntiLambda);
                    pkgr.PackTracks(T2S::PdgCode::PosKaon);
                    pkgr.PackTracks(T2S::PdgCode::PiMinus);
                    pkgr.PackTracks(T2S::PdgCode::PiPlus);
                    break;
                case T2S::ReactionChannel::H:
                    pkgr.PackTracks(T2S::PdgCode::PosKaon);
                    break;
                // anti-channels //
                case T2S::ReactionChannel::AntiA:
                    pkgr.FindV0s(T2S::PdgCode::Lambda);
                    pkgr.FindV0s(T2S::PdgCode::KaonZeroShort);
                    break;
                case T2S::ReactionChannel::AntiD:
                    pkgr.FindV0s(T2S::PdgCode::Lambda);
                    pkgr.PackTracks(T2S::PdgCode::NegKaon);
                    break;
                case T2S::ReactionChannel::AntiE:
                    pkgr.FindV0s(T2S::PdgCode::Lambda);
                    pkgr.PackTracks(T2S::PdgCode::NegKaon);
                    pkgr.PackTracks(T2S::PdgCode::PiMinus);
                    pkgr.PackTracks(T2S::PdgCode::PiPlus);
                    break;
                case T2S::ReactionChannel::AntiH:
                    pkgr.PackTracks(T2S::PdgCode::NegKaon);
                    break;
                // for data //
                case T2S::ReactionChannel::All:
                    pkgr.FindV0s(T2S::PdgCode::AntiLambda);
                    pkgr.FindV0s(T2S::PdgCode::Lambda);
                    pkgr.FindV0s(T2S::PdgCode::KaonZeroShort);
                    pkgr.PackTracks(T2S::PdgCode::NegKaon);
                    pkgr.PackTracks(T2S::PdgCode::PosKaon);
                    pkgr.PackTracks(T2S::PdgCode::PiMinus);
                    pkgr.PackTracks(T2S::PdgCode::PiPlus);
                    break;
            }

            pkgr.EndOfEvent();
        }
        pkgr.EndOfAnalysis();

    } else {

        T2S::Finder finder(settings);
        if (!finder.Initialize()) return 1;

        for (long long i_event{0}; i_event < finder.NumberEventsToRead(); ++i_event) {
            finder.GetEvent(i_event);

            switch (finder.GetReactionChannel()) {
                // standard channels //
                case T2S::ReactionChannel::A:
                    finder.UnpackV0s(T2S::PdgCode::AntiLambda);
                    finder.UnpackV0s(T2S::PdgCode::KaonZeroShort);
                    break;
                case T2S::ReactionChannel::D:
                    finder.UnpackV0s(T2S::PdgCode::AntiLambda);
                    finder.UnpackTracks(T2S::PdgCode::PosKaon);
                    break;
                case T2S::ReactionChannel::E:
                    finder.UnpackV0s(T2S::PdgCode::AntiLambda);
                    finder.UnpackTracks(T2S::PdgCode::PosKaon);
                    finder.UnpackTracks(T2S::PdgCode::PiMinus);
                    finder.UnpackTracks(T2S::PdgCode::PiPlus);
                    break;
                case T2S::ReactionChannel::H:
                    finder.UnpackTracks(T2S::PdgCode::PosKaon);
                    break;
                // anti-channels //
                case T2S::ReactionChannel::AntiA:
                    finder.UnpackV0s(T2S::PdgCode::Lambda);
                    finder.UnpackV0s(T2S::PdgCode::KaonZeroShort);
                    break;
                case T2S::ReactionChannel::AntiD:
                    finder.UnpackV0s(T2S::PdgCode::Lambda);
                    finder.UnpackTracks(T2S::PdgCode::NegKaon);
                    break;
                case T2S::ReactionChannel::AntiE:
                    finder.UnpackV0s(T2S::PdgCode::Lambda);
                    finder.UnpackTracks(T2S::PdgCode::NegKaon);
                    finder.UnpackTracks(T2S::PdgCode::PiMinus);
                    finder.UnpackTracks(T2S::PdgCode::PiPlus);
                    break;
                case T2S::ReactionChannel::AntiH:
                    finder.UnpackTracks(T2S::PdgCode::NegKaon);
                    break;
                default:
                    break;
            }

            finder.Find(finder.GetReactionChannel());

            finder.EndOfEvent();
        }
        finder.EndOfAnalysis();
    }

    return 0;
}
