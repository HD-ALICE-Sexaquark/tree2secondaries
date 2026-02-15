#pragma once

#include <memory>
#include <utility>

#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

#include "App/Settings.hxx"
#include "KalmanFitter/KalmanFitterChannelA.hxx"
#include "KalmanFitter/KalmanFitterChannelD.hxx"
#include "Math/Constants.hxx"
#include "Storage/Flat/FlatChannelA.hxx"
#include "Storage/Flat/FlatChannelD.hxx"
#include "Storage/Flat/FlatEvent.hxx"
#include "Storage/Flat/FlatInjected.hxx"
#include "Storage/Vector/VectorInjected.hxx"
#include "Storage/Vector/VectorTracks.hxx"
#include "Storage/Vector/VectorV0s.hxx"
#include "View/MC/ViewMcInjected.hxx"

namespace Tree2Secondaries {

// Read Vector events and find anti-sexaquarks reactions.
class Finder {
   public:
    Finder(const Finder &) = delete;
    Finder(Finder &&) = delete;
    Finder &operator=(const Finder &) = delete;
    Finder &operator=(Finder &&) = delete;
    ~Finder() = default;

    explicit Finder(Settings settings) : fSettings{std::move(settings)} {}

    [[nodiscard]] EReactionChannel GetReactionChannel() const { return fSettings.ReactionChannel; }

    bool Initialize();
    void ReadInputBranches();

    bool PrepareOutputFile();
    void PrepareOutputHistograms();
    bool PrepareOutputTree();
    void CreateOutputBranches();

    void ProcessEvent();

    bool Injected_PrepareOutputTree();
    void ProcessInjected();

    [[nodiscard]] size_t NumberEventsToRead() const {
        return fSettings.LimitToNEvents ? fSettings.LimitToNEvents : static_cast<size_t>(fInputChain_PackedEvents->GetEntries());
    }
    [[nodiscard]] bool IsMC() const { return fSettings.IsMC; }
    void GetEvent(size_t i_event) { fInputChain_PackedEvents->GetEntry(static_cast<long long>(i_event)); }

    [[nodiscard]] size_t NumberInjected() const { return fInput_Injected.ReactionID->size(); }

    void FindSexaquarks_ChannelA(bool anti_channel);
    void FindSexaquarks_ChannelD(bool anti_channel);
    void Find(EReactionChannel reaction_channel) {
        switch (reaction_channel) {
            case EReactionChannel::A:
                FindSexaquarks_ChannelA(false);
                FindSexaquarks_ChannelA(true);
                break;
            case EReactionChannel::D:
                FindSexaquarks_ChannelD(false);
                FindSexaquarks_ChannelD(true);
                break;
            default:
                return;
        }
    }

    [[nodiscard]] bool FastCuts_ChannelA(const Seeder::Seed &seed_v0a, const Seeder::Seed &seed_v0b, TH1D *cut_flow_hist) const;
    [[nodiscard]] bool FastCuts_ChannelD(const Seeder::Seed &seed_ka, const Seeder::Seed &seed_v0, TH1D *cut_flow_hist) const;

    [[nodiscard]] bool SlowCuts(const KalmanFitter::ChannelA &sexa, TH1D *cut_flow_hist) const;
    [[nodiscard]] bool SlowCuts(const KalmanFitter::ChannelD &sexa, TH1D *cut_flow_hist) const;

    void Store(const KalmanFitter::ChannelA &sexa, bool anti_channel);
    void StoreMC(const View::MC::ChannelA &sexa);
    void StoreDummyMC_ChannelA();

    void Store(const KalmanFitter::ChannelD &sexa, bool anti_channel);
    void StoreMC(const View::MC::ChannelD &sexa);
    void StoreDummyMC_ChannelD();

    void EndOfAnalysis();

   private:
    Settings fSettings;
    std::unique_ptr<TChain> fInputChain_PackedEvents;

    std::unique_ptr<TFile> fOutputFile;
    std::unique_ptr<TTree> fOutputTree;
    std::unique_ptr<TTree> fOutputTree_Injected;

    std::unique_ptr<TH1D> fHist_EventCounter;

    std::unique_ptr<TH1D> fHist_CutFlow;
    std::unique_ptr<TH1D> fHist_CutFlow_AntiChannel;

    // input //

    Storage::Flat::Event fInput_Event;
    Storage::Flat::Coordinates fInput_MC_PV;
    Storage::Vector::Injected fInput_Injected;

    Storage::Vector::V0s fInput_AntiLambdas;
    Storage::Vector::V0s fInput_Lambdas;
    Storage::Vector::V0s fInput_KaonsZeroShort;

    Storage::Vector::Tracks fInput_NegKaons;
    Storage::Vector::Tracks fInput_PosKaons;
    Storage::Vector::Tracks fInput_PiMinus;
    Storage::Vector::Tracks fInput_PiPlus;

    Storage::Vector::MC_V0s fInput_MC_AntiLambdas;
    Storage::Vector::MC_V0s fInput_MC_Lambdas;
    Storage::Vector::MC_V0s fInput_MC_KaonsZeroShort;

    Storage::Vector::MC_Tracks fInput_MC_NegKaons;
    Storage::Vector::MC_Tracks fInput_MC_PosKaons;
    Storage::Vector::MC_Tracks fInput_MC_PiMinus;
    Storage::Vector::MC_Tracks fInput_MC_PiPlus;

    // output //

    Storage::Flat::Injected fOutput_Injected;

    Storage::Flat::ChannelA fOutput_ChannelA;
    Storage::Flat::ChannelD fOutput_ChannelD;

    Storage::Flat::MC_ChannelA fOutput_MC_ChannelA;
    Storage::Flat::MC_ChannelD fOutput_MC_ChannelD;
};

}  // namespace Tree2Secondaries
