#include "Analysis/Manager.hxx"
#include "Analysis/Settings.hxx"
#include "Math/Constants.hxx"
#include "Utilities/Parser.hxx"

namespace T2S = Tree2Secondaries;

int main(int argc, char *argv[]) {

    T2S::Settings settings;
    T2S::Parser parser(argc, argv, settings, "Tree2Secondaries");
    if (parser.HelpOrError) return parser.ExitCode;
    settings.Print();

    T2S::Analysis::Manager mgr(settings);
    if (!mgr.Initialize()) return 1;

    for (long long i_event{0}; i_event < mgr.NumberEventsToRead(); ++i_event) {
        mgr.GetEvent(i_event);

        mgr.ProcessEvent();
        if (mgr.IsSignalMC()) mgr.ProcessInjected();
        mgr.ProcessTracks();

        switch (mgr.GetReactionChannel()) {
            // standard channels //
            case T2S::ReactionChannel::A:
                mgr.FindV0s(T2S::PdgCode::AntiLambda);
                mgr.FindV0s(T2S::PdgCode::KaonZeroShort);
                break;
            case T2S::ReactionChannel::D:
                mgr.FindV0s(T2S::PdgCode::AntiLambda);
                mgr.StoreTracks(T2S::PdgCode::PosKaon);
                break;
            case T2S::ReactionChannel::E:
                mgr.FindV0s(T2S::PdgCode::AntiLambda);
                mgr.StoreTracks(T2S::PdgCode::PosKaon);
                mgr.StoreTracks(T2S::PdgCode::PiMinus);
                mgr.StoreTracks(T2S::PdgCode::PiPlus);
                break;
            case T2S::ReactionChannel::H:
                mgr.StoreTracks(T2S::PdgCode::PosKaon);
                break;
            // anti-channels //
            case T2S::ReactionChannel::AntiA:
                mgr.FindV0s(T2S::PdgCode::Lambda);
                mgr.FindV0s(T2S::PdgCode::KaonZeroShort);
                break;
            case T2S::ReactionChannel::AntiD:
                mgr.FindV0s(T2S::PdgCode::Lambda);
                mgr.StoreTracks(T2S::PdgCode::NegKaon);
                break;
            case T2S::ReactionChannel::AntiE:
                mgr.FindV0s(T2S::PdgCode::Lambda);
                mgr.StoreTracks(T2S::PdgCode::NegKaon);
                mgr.StoreTracks(T2S::PdgCode::PiMinus);
                mgr.StoreTracks(T2S::PdgCode::PiPlus);
                break;
            case T2S::ReactionChannel::AntiH:
                mgr.StoreTracks(T2S::PdgCode::NegKaon);
                break;
            // for data //
            case T2S::ReactionChannel::All:
                mgr.FindV0s(T2S::PdgCode::AntiLambda);
                mgr.FindV0s(T2S::PdgCode::Lambda);
                mgr.FindV0s(T2S::PdgCode::KaonZeroShort);
                mgr.StoreTracks(T2S::PdgCode::NegKaon);
                mgr.StoreTracks(T2S::PdgCode::PosKaon);
                mgr.StoreTracks(T2S::PdgCode::PiMinus);
                mgr.StoreTracks(T2S::PdgCode::PiPlus);
                break;
        }

        mgr.EndOfEvent();
    }
    mgr.EndOfAnalysis();

    return 0;
}
