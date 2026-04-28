#pragma once

#include <cstddef>
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
#include "Storage/Schema/SchemaFlat.hxx"
#include "Storage/Schema/SchemaVector.hxx"
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
    void StoreInjected(const View::MC::Injected &sexa);

    [[nodiscard]] std::size_t NumberEventsToRead() const {
        return fSettings.LimitToNEvents ? fSettings.LimitToNEvents : static_cast<std::size_t>(fInputChain_PackedEvents->GetEntries());
    }
    [[nodiscard]] bool IsMC() const { return fSettings.IsMC; }
    void GetEvent(std::size_t i_event) { fInputChain_PackedEvents->GetEntry(static_cast<long long>(i_event)); }

    [[nodiscard]] std::size_t NumberInjected() const { return fInput_Injected.reaction_id->size(); }

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

    void Store(const KalmanFitter::ChannelD &sexa, bool anti_channel);
    void StoreMC(const View::MC::ChannelD &sexa);

    void EndOfAnalysis();

   private:
    void Assign(Schema::Flat::Track &out, const View::Rec::Track &t);
    void Assign(Schema::Flat::V0 &out, const View::Rec::V0 &v);
    void Assign(Schema::Flat::MC_Track &out, const View::MC::PackedTrack &t, bool ascendants_info);
    void Assign(Schema::Flat::MC_V0 &out, const View::MC::PackedV0 &t);

    Settings fSettings;
    std::unique_ptr<TChain> fInputChain_PackedEvents;

    std::unique_ptr<TFile> fOutputFile;
    std::unique_ptr<TTree> fOutputTree;
    std::unique_ptr<TTree> fOutputTree_Injected;

    std::unique_ptr<TH1D> fHist_EventCounter;

    std::unique_ptr<TH1D> fHist_CutFlow;
    std::unique_ptr<TH1D> fHist_CutFlow_AntiChannel;

    // input //

    Schema::Flat::Event fInput_Event;
    Schema::Vector::Injected fInput_Injected;

    Schema::Vector::V0s fInput_AntiLambdas;
    Schema::Vector::V0s fInput_Lambdas;
    Schema::Vector::V0s fInput_KaonsZeroShort;

    Schema::Vector::Tracks fInput_NegKaons;
    Schema::Vector::Tracks fInput_PosKaons;
    Schema::Vector::Tracks fInput_PiMinus;
    Schema::Vector::Tracks fInput_PiPlus;

    Schema::Vector::MC_V0s fInput_MC_AntiLambdas;
    Schema::Vector::MC_V0s fInput_MC_Lambdas;
    Schema::Vector::MC_V0s fInput_MC_KaonsZeroShort;

    Schema::Vector::MC_Tracks fInput_MC_NegKaons;
    Schema::Vector::MC_Tracks fInput_MC_PosKaons;
    Schema::Vector::MC_Tracks fInput_MC_PiMinus;
    Schema::Vector::MC_Tracks fInput_MC_PiPlus;

    // output //

    Schema::Flat::Injected fOutput_Injected;

    Schema::Flat::ChannelA fOutput_ChannelA;
    Schema::Flat::ChannelD fOutput_ChannelD;

    Schema::Flat::MC_ChannelA fOutput_MC_ChannelA;
    Schema::Flat::MC_ChannelD fOutput_MC_ChannelD;
};

}  // namespace Tree2Secondaries
