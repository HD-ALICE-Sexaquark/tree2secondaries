#ifndef T2S_ANALYSIS_MANAGER_HXX
#define T2S_ANALYSIS_MANAGER_HXX

#include <memory>
#include <utility>

#include "Math/Point3D.h"
#include "TFile.h"
#include "TTree.h"

#include "Analysis/InputFormat.hxx"
#include "Analysis/OutputFormat.hxx"
#include "Analysis/Settings.hxx"
#include "Math/Constants.hxx"
#include "Secondary/Charged.hxx"
#include "Secondary/Neutral.hxx"
// #include "Secondary/True.hxx"

namespace Tree2Secondaries::Analysis {

class Manager {
   public:
    Manager(const Manager &) = delete;
    Manager(Manager &&) = delete;
    Manager &operator=(const Manager &) = delete;
    Manager &operator=(Manager &&) = delete;
    ~Manager() = default;

    explicit Manager(Settings settings) : fSettings{std::move(settings)} {}

    bool Initialize();
    bool OpenInputFile();
    bool LoadInputTree();
    void ConnectInputBranches();

    void ConnectBranchesEvents();
    void ConnectBranchesInjected();
    void ConnectBranchesMC();
    void ConnectBranchesTracks();

    bool PrepareOutputFile();
    bool PrepareOutputTree();
    void PrepareOutputBranches();

    void PrepareBranchesEvents();
    void PrepareBranchesV0s();

    [[nodiscard]] long long NumberEventsToRead() const { return fSettings.LimitToNEvents ? fSettings.LimitToNEvents : fEventsTree->GetEntries(); }
    [[nodiscard]] bool IsMC() const { return fSettings.IsMC; }
    [[nodiscard]] bool IsSignalMC() const { return fSettings.IsSignalMC; }
    void GetEvent(long long i_event) { fEventsTree->GetEntry(i_event); }

    [[nodiscard]] int NumberMC() const { return static_cast<int>(fInput_MC.PdgCode->size()); }
    [[nodiscard]] int NumberInjected() const { return static_cast<int>(fInput_Injected.ReactionID->size()); }
    [[nodiscard]] int NumberTracks() const { return static_cast<int>(fInput_Tracks.Px->size()); }

    void ProcessEvent();
    void ProcessInjected();
    void ProcessMC();
    void ProcessTracks();

    void FindV0s(int pdg_code_v0, int pdg_code_neg, int pdg_code_pos);
    [[nodiscard]] bool PassesLambdaCuts(const std::shared_ptr<Neutral> &this_v0) const;
    [[nodiscard]] bool PassesKaonZeroCuts(const std::shared_ptr<Neutral> &this_v0) const;
    [[nodiscard]] bool PassesV0Cuts(const std::shared_ptr<Neutral> &this_v0, int pdg_code_v0) const {
        if (std::abs(pdg_code_v0) == PdgCode::Lambda) return PassesLambdaCuts(this_v0);
        return PassesKaonZeroCuts(this_v0);
    }
    void Store(const std::shared_ptr<Neutral> &v0, int pdg_code_v0);

    void EndOfEvent();
    void EndOfAnalysis();

   private:
    Settings fSettings;
    std::unique_ptr<TFile> fInputFile;
    std::unique_ptr<TTree> fEventsTree;

    std::unique_ptr<TFile> fOutputFile;
    std::unique_ptr<TTree> fOutputTree;

    // input structs //
    Struct::Event fInput_Event;
    Input::Injected fInput_Injected;
    Input::MC fInput_MC;
    Input::Tracks fInput_Tracks;

    // helpers //
    Helper::Propagator fPropagator{0.};

    // transitory containers //
    // std::vector<std::shared_ptr<True>> fMCParticles;

    std::vector<std::shared_ptr<Charged>> fAntiProtons;
    std::vector<std::shared_ptr<Charged>> fProtons;
    std::vector<std::shared_ptr<Charged>> fNegKaons;
    std::vector<std::shared_ptr<Charged>> fPosKaons;
    std::vector<std::shared_ptr<Charged>> fPiMinus;
    std::vector<std::shared_ptr<Charged>> fPiPlus;

    std::vector<std::shared_ptr<Neutral>> fAntiLambdas;
    std::vector<std::shared_ptr<Neutral>> fLambdas;
    std::vector<std::shared_ptr<Neutral>> fNeutralKaons;  // KaonZeroShort

    // output structs //
    Struct::Event fOutput_Event;
    Output::V0s fOutput_V0s;
};

}  // namespace Tree2Secondaries::Analysis

#endif  // T2S_ANALYSIS_MANAGER_HXX
