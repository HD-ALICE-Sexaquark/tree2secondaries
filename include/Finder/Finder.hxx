#pragma once

#include <cstddef>
#include <memory>

#include <TChain.h>
#include <TFile.h>
#include <TH1.h>

#include "App/Settings.hxx"
#include "KalmanFitter/KalmanFitterChannelA.hxx"
#include "KalmanFitter/KalmanFitterChannelD.hxx"
#include "MonteCarlo/MonteCarloChannelA.hxx"
#include "MonteCarlo/MonteCarloChannelD.hxx"
#include "Schema/SchemaFlatChannelA.hxx"
#include "Schema/SchemaFlatChannelD.hxx"
#include "Schema/SchemaFlatInjected.hxx"
#include "Schema/SchemaFlatTrack.hxx"
#include "Schema/SchemaFlatV0.hxx"
#include "Schema/SchemaVectorInjected.hxx"
#include "Schema/SchemaVectorMcTracks.hxx"
#include "Schema/SchemaVectorMcV0s.hxx"
#include "Schema/SchemaVectorTracks.hxx"
#include "Schema/SchemaVectorV0s.hxx"
#include "View/ViewVectorMcTracks.hxx"
#include "View/ViewVectorMcV0s.hxx"
#include "View/ViewVectorTracks.hxx"
#include "View/ViewVectorV0s.hxx"

class TTree;

namespace Tree2Secondaries {

// Read Vector events and find anti-sexaquarks reactions.
class Finder {
   public:
    Finder(const Finder &) = delete;
    Finder(Finder &&) = delete;
    Finder &operator=(const Finder &) = delete;
    Finder &operator=(Finder &&) = delete;
    ~Finder() = default;

    explicit Finder(const Settings &settings) : fSettings{settings} {}

    bool Initialize();
    void ReadInputBranches();

    bool PrepareOutputFile();
    void PrepareOutputHistograms();
    bool PrepareOutputTree();
    void CreateOutputBranches();

    void ProcessEvent();

    bool Injected_PrepareOutputTree();
    void ProcessInjected();

    [[nodiscard]] long long NumberEventsToRead() const {
        long long total_entries = fInputChain_PackedEvents->GetEntries();
        return 0 < fSettings.LimitToNEvents && fSettings.LimitToNEvents < total_entries ? fSettings.LimitToNEvents : total_entries;
    }
    void GetEvent(long long i_event) { fInputChain_PackedEvents->GetEntry(i_event); }

    [[nodiscard]] std::size_t NumberInjected() const { return fInput_Injected.reaction_id->size(); }

    void FindSexaquarks_ChannelA(bool anti_channel);
    void FindSexaquarks_ChannelD(bool anti_channel);
    void Find() {
        switch (fSettings.ReactionChannel.name) {
            case 'A': {
                FindSexaquarks_ChannelA(false);
                FindSexaquarks_ChannelA(true);
                break;
            }
            case 'D': {
                FindSexaquarks_ChannelD(false);
                FindSexaquarks_ChannelD(true);
                break;
            }
            default:
                return;
        }
    }

    [[nodiscard]] bool FastCuts_ChannelA(const Seeder::Seed &seed_v0a, const Seeder::Seed &seed_v0b, TH1D *cut_flow_hist) const;
    [[nodiscard]] bool FastCuts_ChannelD(const Seeder::Seed &seed_ka, const Seeder::Seed &seed_v0, TH1D *cut_flow_hist) const;

    [[nodiscard]] bool SlowCuts(const KalmanFitter::ChannelA &sexa, TH1D *cut_flow_hist) const;
    [[nodiscard]] bool SlowCuts(const KalmanFitter::ChannelD &sexa, TH1D *cut_flow_hist) const;

    void EndOfAnalysis();

   private:
    void Assign(const View::VecTracks &v, Schema::FlatTrack &out);
    void Assign(const View::VecV0s &v, Schema::FlatV0 &out);
    void Assign(const View::VecMcTracks &v, Schema::FlatMcTrack &out, bool ascendants_info);
    void Assign(const View::VecMcV0s &v, Schema::FlatMcV0 &out);

    void Assign(const KalmanFitter::ChannelA &sexa, bool anti_channel);
    void Assign(const std::optional<MonteCarlo::ChannelA> &mc_sexa);

    void Assign(const KalmanFitter::ChannelD &sexa, bool anti_channel);
    void Assign(const std::optional<MonteCarlo::ChannelD> &mc_sexa);

    const Settings &fSettings;
    std::unique_ptr<TChain> fInputChain_PackedEvents;

    std::unique_ptr<TFile> fOutputFile;
    std::unique_ptr<TTree> fOutputTree;
    std::unique_ptr<TTree> fOutputTree_Injected;

    std::unique_ptr<TH1D> fHist_EventCounter;

    std::unique_ptr<TH1D> fHist_CutFlow;
    std::unique_ptr<TH1D> fHist_CutFlow_AntiChannel;

    // input //

    Schema::FlatEvent fInput_Event;
    Schema::VecInjected fInput_Injected;

    Schema::VecV0s fInput_AntiLambdas;
    Schema::VecV0s fInput_Lambdas;
    Schema::VecV0s fInput_KaonsZeroShort;

    Schema::VecTracks fInput_NegKaons;
    Schema::VecTracks fInput_PosKaons;

    Schema::VecMcV0s fInput_MC_AntiLambdas;
    Schema::VecMcTracks fInput_MC_AntiLambdas_Neg;
    Schema::VecMcTracks fInput_MC_AntiLambdas_Pos;
    Schema::VecMcV0s fInput_MC_Lambdas;
    Schema::VecMcTracks fInput_MC_Lambdas_Neg;
    Schema::VecMcTracks fInput_MC_Lambdas_Pos;
    Schema::VecMcV0s fInput_MC_KaonsZeroShort;
    Schema::VecMcTracks fInput_MC_KaonsZeroShort_Neg;
    Schema::VecMcTracks fInput_MC_KaonsZeroShort_Pos;

    Schema::VecMcTracks fInput_MC_NegKaons;
    Schema::VecMcTracks fInput_MC_PosKaons;
    Schema::VecMcTracks fInput_MC_PiMinus;
    Schema::VecMcTracks fInput_MC_PiPlus;

    // output //

    Schema::FlatInjected fOutput_Injected;

    Schema::FlatChannelA fOutput_ChannelA;
    Schema::FlatChannelD fOutput_ChannelD;

    Schema::FlatMcChannelA fOutput_MC_ChannelA;
    Schema::FlatMcV0 fOutput_MC_ChannelA_V0A;
    Schema::FlatMcTrack fOutput_MC_ChannelA_V0A_Neg;
    Schema::FlatMcTrack fOutput_MC_ChannelA_V0A_Pos;
    Schema::FlatMcV0 fOutput_MC_ChannelA_V0B;
    Schema::FlatMcTrack fOutput_MC_ChannelA_V0B_Neg;
    Schema::FlatMcTrack fOutput_MC_ChannelA_V0B_Pos;

    Schema::FlatMcChannelD fOutput_MC_ChannelD;
    Schema::FlatMcV0 fOutput_MC_ChannelD_V0;
    Schema::FlatMcTrack fOutput_MC_ChannelD_V0_Neg;
    Schema::FlatMcTrack fOutput_MC_ChannelD_V0_Pos;
    Schema::FlatMcTrack fOutput_MC_ChannelD_Kaon;
};

}  // namespace Tree2Secondaries
