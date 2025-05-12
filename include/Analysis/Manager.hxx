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

    [[nodiscard]] ReactionChannel GetReactionChannel() const { return fSettings.Channel; }

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
    void CreateOutputBranches();

    void CreateOutputBranchesEvents();
    void CreateOutputBranchesV0s(std::string_view v0_sv, Output::V0s &out_branches);
    void CreateOutputBranchesTracks(std::string_view charged_sv, Output::Tracks &out_branches);

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

    void StoreTracks(int pdg_code) {
        if (pdg_code == PdgCode::AntiProton) StoreTracks(fVec_AntiProtons, fOutput_AntiProtons);
        if (pdg_code == PdgCode::Proton) StoreTracks(fVec_Protons, fOutput_Protons);
        if (pdg_code == PdgCode::NegKaon) StoreTracks(fVec_NegKaons, fOutput_NegKaons);
        if (pdg_code == PdgCode::PosKaon) StoreTracks(fVec_PosKaons, fOutput_PosKaons);
        if (pdg_code == PdgCode::PiMinus) StoreTracks(fVec_PiMinus, fOutput_PiMinus);
        if (pdg_code == PdgCode::PiPlus) StoreTracks(fVec_PiPlus, fOutput_PiPlus);
    }

    void FindV0s(int pdg_code_v0, int pdg_code_neg, int pdg_code_pos);

    [[nodiscard]] bool PassesV0Cuts(const std::shared_ptr<Neutral> &this_v0, int pdg_code_v0) const {
        if (std::abs(pdg_code_v0) == PdgCode::Lambda) return PassesLambdaCuts(this_v0);
        return PassesKaonZeroCuts(this_v0);
    }

    void Store(const std::shared_ptr<Neutral> &v0, int pdg_code_v0) {
        if (pdg_code_v0 == PdgCode::AntiLambda)
            Store(v0, fOutput_AntiLambdas);
        else if (pdg_code_v0 == PdgCode::KaonZeroShort)
            Store(v0, fOutput_NeutralKaons);
        else
            Store(v0, fOutput_Lambdas);
    }

    void EndOfEvent();
    void EndOfAnalysis();

   private:
    static void StoreTracks(const std::vector<std::shared_ptr<Charged>> &charged_vec, const Output::Tracks &out_branch);
    [[nodiscard]] bool PassesLambdaCuts(const std::shared_ptr<Neutral> &this_v0) const;
    [[nodiscard]] bool PassesKaonZeroCuts(const std::shared_ptr<Neutral> &this_v0) const;

    static void Store(const std::shared_ptr<Neutral> &v0, const Output::V0s &out_branches);

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

    std::vector<std::shared_ptr<Charged>> fVec_AntiProtons;
    std::vector<std::shared_ptr<Charged>> fVec_Protons;
    std::vector<std::shared_ptr<Charged>> fVec_NegKaons;
    std::vector<std::shared_ptr<Charged>> fVec_PosKaons;
    std::vector<std::shared_ptr<Charged>> fVec_PiMinus;
    std::vector<std::shared_ptr<Charged>> fVec_PiPlus;

    // output structs //
    Struct::Event fOutput_Event;
    Output::V0s fOutput_AntiLambdas;
    Output::V0s fOutput_Lambdas;
    Output::V0s fOutput_NeutralKaons;
    Output::Tracks fOutput_AntiProtons;
    Output::Tracks fOutput_Protons;
    Output::Tracks fOutput_NegKaons;
    Output::Tracks fOutput_PosKaons;
    Output::Tracks fOutput_PiMinus;
    Output::Tracks fOutput_PiPlus;
};

}  // namespace Tree2Secondaries::Analysis

#endif  // T2S_ANALYSIS_MANAGER_HXX
