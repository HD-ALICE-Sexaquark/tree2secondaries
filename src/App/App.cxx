#include "Analysis/Manager.hxx"
#include "App/Parser.hxx"
#include "App/Settings.hxx"
#include "Math/Constants.hxx"

using namespace Tree2Secondaries;

int main(int argc, char *argv[]) {

    Settings settings{"AnalysisResults.root",   //
                      "SecondaryResults.root",  //
                      false,                    //
                      false,                    //
                      0};
    Parser parser(argc, argv);
    if (!parser.Parse(settings)) return parser.exitCode;
    settings.Print();

    Analysis::Manager mgr(settings);
    if (!mgr.Initialize()) return 1;

    for (long long i_event{0}; i_event < mgr.NumberEventsToRead(); i_event++) {
        mgr.GetEvent(i_event);
        mgr.ProcessEvent();
        // if (mgr.IsMC()) {
        // mgr.ProcessMC();
        // if (mgr.IsSignalMC()) mgr.ProcessInjected(); // PENDING
        // }
        mgr.ProcessTracks();
        mgr.FindV0s(PdgCode::Lambda, PdgCode::PiMinus, PdgCode::Proton);
        mgr.FindV0s(PdgCode::AntiLambda, PdgCode::AntiProton, PdgCode::PiPlus);
        mgr.FindV0s(PdgCode::KaonZeroShort, PdgCode::PiMinus, PdgCode::PiPlus);

        // mgr.FindSexaquarks(Channel::A);
        // mgr.FindSexaquarks(Channel::D);

        mgr.EndOfEvent();
    }
    // mgr.WriteOutputFile();
    mgr.EndOfAnalysis();

    return 0;
}
