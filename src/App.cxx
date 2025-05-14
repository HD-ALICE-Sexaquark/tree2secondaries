#include "Analysis/Manager.hxx"
#include "Analysis/Settings.hxx"
#include "Math/Constants.hxx"
#include "Utilities/Parser.hxx"

using namespace Tree2Secondaries;

int main(int argc, char *argv[]) {

    Settings settings;
    Parser parser(argc, argv, settings, "II. Tree2Secondaries");
    if (parser.HelpOrError) return parser.ExitCode;
    settings.Print();

    Analysis::Manager mgr(settings);
    if (!mgr.Initialize()) return 1;

    for (long long i_event{0}; i_event < mgr.NumberEventsToRead(); i_event++) {
        mgr.GetEvent(i_event);

        mgr.ProcessEvent();
        if (mgr.IsSignalMC()) mgr.ProcessInjected();
        mgr.ProcessTracks();

        switch (mgr.GetReactionChannel()) {
            // standard channels //
            case ReactionChannel::A:
                mgr.FindV0s(PdgCode::AntiLambda, PdgCode::AntiProton, PdgCode::PiPlus);
                mgr.FindV0s(PdgCode::KaonZeroShort, PdgCode::PiMinus, PdgCode::PiPlus);
                break;
            case ReactionChannel::D:
                mgr.FindV0s(PdgCode::AntiLambda, PdgCode::AntiProton, PdgCode::PiPlus);
                mgr.StoreTracks(PdgCode::PosKaon);
                break;
            case ReactionChannel::E:
                mgr.FindV0s(PdgCode::AntiLambda, PdgCode::AntiProton, PdgCode::PiPlus);
                mgr.StoreTracks(PdgCode::PosKaon);
                mgr.StoreTracks(PdgCode::PiMinus);
                mgr.StoreTracks(PdgCode::PiPlus);
                break;
            case ReactionChannel::H:
                // SUPER PENDING
                break;
            // anti-channels //
            case ReactionChannel::AntiA:
                mgr.FindV0s(PdgCode::Lambda, PdgCode::PiMinus, PdgCode::Proton);
                mgr.FindV0s(PdgCode::KaonZeroShort, PdgCode::PiMinus, PdgCode::PiPlus);
                break;
            case ReactionChannel::AntiD:
                mgr.FindV0s(PdgCode::Lambda, PdgCode::PiMinus, PdgCode::Proton);
                mgr.StoreTracks(PdgCode::NegKaon);
                break;
            case ReactionChannel::AntiE:
                mgr.FindV0s(PdgCode::Lambda, PdgCode::PiMinus, PdgCode::Proton);
                mgr.StoreTracks(PdgCode::NegKaon);
                mgr.StoreTracks(PdgCode::PiMinus);
                mgr.StoreTracks(PdgCode::PiPlus);
                break;
            case ReactionChannel::AntiH:
                // SUPER PENDING
                break;
        }  // end of switch statement

        mgr.EndOfEvent();
    }
    mgr.EndOfAnalysis();

    return 0;
}
