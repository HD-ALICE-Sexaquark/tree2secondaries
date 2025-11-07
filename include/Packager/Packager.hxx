#pragma once

#include <memory>

#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

#include "App/Settings.hxx"
#include "DataFormats/Events.hxx"
#include "DataFormats/Injected.hxx"
#include "DataFormats/PackedEvents.hxx"
#include "Fit/Track.hxx"
#include "Fit/V0.hxx"
#include "Math/Constants.hxx"
#include "References/Events.hxx"

namespace Tree2Secondaries {

// Pack secondary V0s and tracks.
class Packager {
   public:
    Packager(const Packager &) = delete;
    Packager(Packager &&) = delete;
    Packager &operator=(const Packager &) = delete;
    Packager &operator=(Packager &&) = delete;
    ~Packager() = default;

    explicit Packager(Settings settings) : fSettings{std::move(settings)} {}

    [[nodiscard]] EReactionChannel GetReactionChannel() const { return fSettings.ReactionChannel; }

    bool Initialize();
    void ReadInputBranches();

    void ReadBranches_Events() { fInput_Event.ReadBranches_Event(fInputChain_Events.get(), IsMC()); }
    void ReadBranches_Injected() { fInput_Injected.ReadBranches_SOV_Injected(fInputChain_Events.get(), false); }
    void ReadBranches_MC() { fInput_MC.ReadBranches_MCParticles(fInputChain_Events.get()); }
    void ReadBranches_Tracks() { fInput_Tracks.ReadBranches_MCParticles(fInputChain_Events.get(), IsMC()); }

    bool PrepareOutputFile();
    bool PrepareOutputTree();
    void CreateOutputBranches();
    void PrepareOutputHistograms();

    void CreateOutputBranches_Events() {  //
        fOutput_Event.CreateBranches_Event(fOutputTree.get(), IsMC());
    }
    void CreateOutputBranches_Injected() {  //
        fOutput_Injected.CreateBranches_SOV_Injected(fOutputTree.get(), true);
    }
    void CreateOutputBranches_V0s(EParticle pid, DF::Packed::V0s &df) {
        df.CreateBranches_PackedV0s(fOutputTree.get(), Const::Particle_Acronym[pid]);
    }
    void CreateOutputBranches_Tracks(EParticle pid, DF::Packed::Tracks &df) {
        df.CreateBranches_PackedTracks(fOutputTree.get(), Const::Particle_Acronym[pid]);
    }
    void CreateOutputBranches_LinkedV0s(EParticle pid, DF::Packed::LinkedV0s &df) {
        df.CreateBranches_LinkedV0s(fOutputTree.get(), Const::Particle_Acronym[pid]);
    }
    void CreateOutputBranches_LinkedTracks(EParticle pid, DF::Packed::LinkedTracks &df) {
        df.CreateBranches_LinkedTracks(fOutputTree.get(), Const::Particle_Acronym[pid]);
    }

    [[nodiscard]] int NumberEventsToRead() const {
        return fSettings.LimitToNEvents ? fSettings.LimitToNEvents : static_cast<int>(fInputChain_Events->GetEntries());
    }
    [[nodiscard]] bool IsMC() const { return fSettings.IsMC; }
    void GetEvent(int i_event) { fInputChain_Events->GetEntry(i_event); }

    [[nodiscard]] int NumberMC() const { return static_cast<int>(fInput_MC.Px->size()); }
    [[nodiscard]] int NumberInjected() const { return static_cast<int>(fInput_Injected.ReactionID->size()); }
    [[nodiscard]] int NumberTracks() const { return static_cast<int>(fInput_Tracks.Px->size()); }

    void ProcessEvent() {
        fOutput_Event = fInput_Event;
        fHist_EventCounter->Fill(0.);
    }

    void Injected_GetSecondaryVertex();
    void Injected_Store();
    void ProcessInjected() {
        Injected_GetSecondaryVertex();
        Injected_Store();
    }

    void ProcessTracks();

    void PackTracks(EParticle pid);

    void FindV0s(EParticle pid);
    [[nodiscard]] bool PassesCuts(const Fit::V0 &v0, EParticle pid) const {
        switch (pid) {
            case EParticle::AntiLambda:
                return PassesCuts_Lambda(v0, fHist_CutFlow_AntiLambda.get());
            case EParticle::Lambda:
                return PassesCuts_Lambda(v0, fHist_CutFlow_Lambda.get());
            case EParticle::KaonZeroShort:
                return PassesCuts_KaonZeroShort(v0, fHist_CutFlow_KaonZeroShort.get());
            default:
                return false;
        }
    }

    void EndOfEvent();
    void EndOfAnalysis();

   private:
    bool PassesCuts_Proton(const Ref::Track &track, TH1D *cut_flow_hist) const;
    bool PassesCuts_Kaon(const Ref::Track &track, TH1D *cut_flow_hist) const;
    bool PassesCuts_Pion(const Ref::Track &track, TH1D *cut_flow_hist) const;
    void Store(const Fit::Track &track, DF::Packed::Tracks &df);
    void StoreMC(const Ref::MC_Track &mc, DF::Packed::LinkedTracks &df);

    bool PassesCuts_Lambda(const Fit::V0 &v0, TH1D *cut_flow_hist) const;
    bool PassesCuts_KaonZeroShort(const Fit::V0 &v0, TH1D *cut_flow_hist) const;
    void Store(const Fit::V0 &v0, DF::Packed::V0s &df);
    void StoreMC(const Ref::MC_V0 &mc, DF::Packed::LinkedV0s &df);

    Settings fSettings;
    std::unique_ptr<TChain> fInputChain_Events;

    std::unique_ptr<TFile> fOutputFile;
    std::unique_ptr<TTree> fOutputTree;

    std::unique_ptr<TH1D> fHist_EventCounter;

    std::unique_ptr<TH1D> fHist_CutFlow_AntiProton;
    std::unique_ptr<TH1D> fHist_CutFlow_Proton;
    std::unique_ptr<TH1D> fHist_CutFlow_NegKaon;
    std::unique_ptr<TH1D> fHist_CutFlow_PosKaon;
    std::unique_ptr<TH1D> fHist_CutFlow_PiMinus;
    std::unique_ptr<TH1D> fHist_CutFlow_PiPlus;

    std::unique_ptr<TH1D> fHist_CutFlow_AntiLambda;
    std::unique_ptr<TH1D> fHist_CutFlow_Lambda;
    std::unique_ptr<TH1D> fHist_CutFlow_KaonZeroShort;

    // input branches //

    DF::Events::Event fInput_Event;
    DF::SOV::Injected fInput_Injected;

    DF::Events::MCParticles fInput_MC;
    DF::Events::Tracks fInput_Tracks;

    // temporary containers, cleaned after event loop //

    std::vector<float> fVec_SV_X;
    std::vector<float> fVec_SV_Y;
    std::vector<float> fVec_SV_Z;

    std::vector<int> fVec_AntiProtons;
    std::vector<int> fVec_Protons;
    std::vector<int> fVec_NegKaons;
    std::vector<int> fVec_PosKaons;
    std::vector<int> fVec_PiMinus;
    std::vector<int> fVec_PiPlus;

    // output branches //

    DF::Events::Event fOutput_Event;
    DF::SOV::Injected fOutput_Injected;

    DF::Packed::LinkedV0s fOutput_Linked_AntiLambdas;
    DF::Packed::LinkedV0s fOutput_Linked_Lambdas;
    DF::Packed::LinkedV0s fOutput_Linked_KaonsZeroShort;

    DF::Packed::LinkedTracks fOutput_Linked_NegKaons;
    DF::Packed::LinkedTracks fOutput_Linked_PosKaons;
    DF::Packed::LinkedTracks fOutput_Linked_PiMinus;
    DF::Packed::LinkedTracks fOutput_Linked_PiPlus;

    DF::Packed::V0s fOutput_AntiLambdas;
    DF::Packed::V0s fOutput_Lambdas;
    DF::Packed::V0s fOutput_KaonsZeroShort;

    DF::Packed::Tracks fOutput_NegKaons;
    DF::Packed::Tracks fOutput_PosKaons;
    DF::Packed::Tracks fOutput_PiMinus;
    DF::Packed::Tracks fOutput_PiPlus;
};

}  // namespace Tree2Secondaries
