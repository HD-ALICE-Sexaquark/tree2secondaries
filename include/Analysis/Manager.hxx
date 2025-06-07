#ifndef T2S_ANALYSIS_MANAGER_HXX
#define T2S_ANALYSIS_MANAGER_HXX

#include <memory>
#include <utility>

#include "TChain.h"
#include "TFile.h"

#include "Analysis/InputFormat.hxx"
#include "Analysis/OutputFormat.hxx"
#include "Analysis/Settings.hxx"
#include "Analysis/TruthHandler.hxx"
#include "Math/Constants.hxx"
#include "Math/KFWrapper.hxx"

namespace Tree2Secondaries::Analysis {

class Manager {
   public:
    Manager(const Manager &) = delete;
    Manager(Manager &&) = delete;
    Manager &operator=(const Manager &) = delete;
    Manager &operator=(Manager &&) = delete;
    ~Manager() = default;

    explicit Manager(Settings settings) : fSettings{std::move(settings)} {}

    ReactionChannel GetReactionChannel() const { return fSettings.Channel; }

    bool Initialize();
    void ConnectInputBranches();

    void ConnectBranchesEvents();
    void ConnectBranchesInjected();
    void ConnectBranchesMC();
    void ConnectBranchesTracks();

    bool PrepareOutputFile();
    bool PrepareOutputTree();
    void CreateOutputBranches();

    void CreateOutputBranchesEvents();
    void CreateOutputBranchesInjected();
    void CreateOutputBranchesV0s(std::string_view v0_sv, OutputSOA::V0s &out_branches);
    void CreateOutputBranchesTracks(std::string_view charged_sv, OutputSOA::Tracks &out_branches);

    long long NumberEventsToRead() const { return fSettings.LimitToNEvents ? fSettings.LimitToNEvents : fEventsTree->GetEntries(); }
    bool IsMC() const { return fSettings.IsMC; }
    bool IsSignalMC() const { return fSettings.IsSignalMC; }
    void GetEvent(long long i_event) { fEventsTree->GetEntry(i_event); }

    size_t NumberInjected() const { return fInput_Injected.ReactionID->size(); }
    size_t NumberTracks() const { return fInput_Tracks.Px->size(); }

    void ProcessEvent();
    void ProcessInjected();
    void ProcessMC();
    void ProcessTracks();

    void StoreTracks(PdgCode pdg_code);

    void FindV0s(PdgCode pdg_code_v0);
    bool PassesV0Cuts(const KF::V0 &v0, PdgCode pdg_code_v0) const {
        switch (pdg_code_v0) {
            case PdgCode::AntiLambda:
            case PdgCode::Lambda:
                return PassesLambdaCuts(v0);
            case PdgCode::KaonZeroShort:
                return PassesKaonZeroCuts(v0);
            default:
                return false;
        }
    }

    void EndOfEvent();
    void EndOfAnalysis();

   private:
    std::array<double, 6> PackParams_KF(size_t esd_idx) { return KF::PackParams(fInput_Tracks, esd_idx); }
    std::array<float, 5> PackParams_ALICE(size_t esd_idx) { return ALICE::PackParams(fInput_Tracks, esd_idx); }
    std::array<float, 15> PackCovMatrix_ALICE(size_t esd_idx) { return ALICE::PackCovMatrix(fInput_Tracks, esd_idx); }

    void Store(const KF::Track &track, OutputSOA::Tracks &out_branches);

    bool PassesLambdaCuts(const KF::V0 &v0) const;
    bool PassesKaonZeroCuts(const KF::V0 &v0) const;
    void Store(const KF::V0 &v0, OutputSOA::V0s &out_branches);

    Settings fSettings;
    std::unique_ptr<TChain> fEventsTree;

    std::unique_ptr<TFile> fOutputFile;
    std::unique_ptr<TTree> fOutputTree;

    // input structs //

    Struct::Event fInput_Event;
    Struct::InjectedSOA fInput_Injected;

    InputSOA::Tracks fInput_Tracks;

    // helpers //

    Helper::TruthHandler fTruthHandler;

    // indices //

    std::vector<size_t> fVec_AntiProtons;
    std::vector<size_t> fVec_Protons;
    std::vector<size_t> fVec_NegKaons;
    std::vector<size_t> fVec_PosKaons;
    std::vector<size_t> fVec_PiMinus;
    std::vector<size_t> fVec_PiPlus;

    // output structs //

    Struct::Event fOutput_Event;
    Struct::InjectedSOA fOutput_Injected;

    OutputSOA::V0s fOutput_AntiLambdas;
    OutputSOA::V0s fOutput_Lambdas;
    OutputSOA::V0s fOutput_KaonsZeroShort;

    OutputSOA::Tracks fOutput_NegKaons;
    OutputSOA::Tracks fOutput_PosKaons;
    OutputSOA::Tracks fOutput_PiMinus;
    OutputSOA::Tracks fOutput_PiPlus;
};

}  // namespace Tree2Secondaries::Analysis

#endif  // T2S_ANALYSIS_MANAGER_HXX
