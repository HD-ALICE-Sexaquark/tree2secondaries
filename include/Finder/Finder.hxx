#pragma once

#include <memory>
#include <utility>

#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

#include "App/Settings.hxx"
#include "DataFormats/ChannelA.hxx"
#include "DataFormats/ChannelD.hxx"
#include "DataFormats/Events.hxx"
#include "DataFormats/Injected.hxx"
#include "DataFormats/PackedEvents.hxx"
#include "Fit/ChannelA.hxx"
#include "Fit/ChannelD.hxx"
#include "Math/Constants.hxx"
#include "References/Injected.hxx"

namespace Tree2Secondaries {

// Read packed events and find anti-sexaquarks reactions.
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

    void ReadBranches_Events() {  //
        fInput_Event.ReadBranches_Event(fInputChain_PackedEvents.get(), IsMC());
    }
    void ReadBranches_Injected() {  //
        fInput_Injected.ReadBranches_SOV_Injected(fInputChain_PackedEvents.get(), true);
    }
    void ReadBranches_V0s(EParticle pid, DF::Packed::V0s &df) {
        df.ReadBranches_PackedV0s(fInputChain_PackedEvents.get(), Const::Particle_Acronym[pid]);
    }
    void ReadBranches_Tracks(EParticle pid, DF::Packed::Tracks &df) {
        df.ReadBranches_PackedTracks(fInputChain_PackedEvents.get(), Const::Particle_Acronym[pid]);
    }
    void ReadBranches_LinkedV0s(EParticle pid, DF::Packed::LinkedV0s &df) {
        df.ReadBranches_LinkedV0s(fInputChain_PackedEvents.get(), Const::Particle_Acronym[pid]);
    }
    void ReadBranches_LinkedTracks(EParticle pid, DF::Packed::LinkedTracks &df) {
        df.ReadBranches_LinkedTracks(fInputChain_PackedEvents.get(), Const::Particle_Acronym[pid]);
    }

    bool PrepareOutputFile();
    bool PrepareOutputTree();
    void PrepareOutputHistograms();

    void CreateOutputBranches(DF::Found::ChannelA &df) { df.CreateBranches_ChannelA(fOutputTree.get(), IsMC()); }
    void CreateOutputBranches(DF::Found::ChannelD &df) { df.CreateBranches_ChannelD(fOutputTree.get(), IsMC()); }
    void CreateOutputBranches(DF::Found::MC_ChannelA &df) { df.CreateBranches_MC_ChannelA(fOutputTree.get()); }
    void CreateOutputBranches(DF::Found::MC_ChannelD &df) { df.CreateBranches_MC_ChannelD(fOutputTree.get()); }
    void CreateOutputBranches() {
        switch (GetReactionChannel()) {
            // standard channels //
            case EReactionChannel::A:
                CreateOutputBranches(fOutput_ChannelA);
                if (IsMC()) CreateOutputBranches(fOutput_MC_ChannelA);
                break;
            case EReactionChannel::D:
                CreateOutputBranches(fOutput_ChannelD);
                if (IsMC()) CreateOutputBranches(fOutput_MC_ChannelD);
                break;
            default:
                break;
        }  // end of switch statement
    }

    void ProcessEvent() { fHist_EventCounter->Fill(0.); }

    bool Injected_PrepareOutputTree();
    void Injected_CreateOutputBranches() { fOutput_Injected.CreateBranches_Flat_Injected(fOutputTree_Injected.get()); };
    void Injected_FlattenAndStore();

    [[nodiscard]] int NumberEventsToRead() const {
        return fSettings.LimitToNEvents ? fSettings.LimitToNEvents : static_cast<int>(fInputChain_PackedEvents->GetEntries());
    }
    [[nodiscard]] bool IsMC() const { return fSettings.IsMC; }
    void GetEvent(int i_event) { fInputChain_PackedEvents->GetEntry(i_event); }

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

    [[nodiscard]] bool PassesCuts(const Fit::ChannelA &sexa, TH1D *cut_flow_hist) const;
    [[nodiscard]] bool PassesCuts(const Fit::ChannelD &sexa, TH1D *cut_flow_hist) const;

    void Store(const Fit::ChannelA &sexa, bool anti_channel);
    void Store(const Fit::ChannelD &sexa, bool anti_channel);
    void StoreMC(const Ref::ChannelA &sexa);
    void StoreMC(const Ref::ChannelD &sexa);

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

    DF::Events::Event fInput_Event;
    DF::SOV::Injected fInput_Injected;

    DF::Packed::V0s fInput_AntiLambdas;
    DF::Packed::V0s fInput_Lambdas;
    DF::Packed::V0s fInput_KaonsZeroShort;

    DF::Packed::Tracks fInput_NegKaons;
    DF::Packed::Tracks fInput_PosKaons;
    DF::Packed::Tracks fInput_PiMinus;
    DF::Packed::Tracks fInput_PiPlus;

    DF::Packed::LinkedV0s fInput_Linked_AntiLambdas;
    DF::Packed::LinkedV0s fInput_Linked_Lambdas;
    DF::Packed::LinkedV0s fInput_Linked_KaonsZeroShort;

    DF::Packed::LinkedTracks fInput_Linked_NegKaons;
    DF::Packed::LinkedTracks fInput_Linked_PosKaons;
    DF::Packed::LinkedTracks fInput_Linked_PiMinus;
    DF::Packed::LinkedTracks fInput_Linked_PiPlus;

    // output //

    DF::Found::ChannelA fOutput_ChannelA;
    DF::Found::ChannelD fOutput_ChannelD;

    DF::Flat::Injected fOutput_Injected;

    DF::Found::MC_ChannelA fOutput_MC_ChannelA;
    DF::Found::MC_ChannelD fOutput_MC_ChannelD;
};

}  // namespace Tree2Secondaries
