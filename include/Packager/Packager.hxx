#pragma once

#include <memory>

#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

#include "App/Settings.hxx"
#include "Math/Constants.hxx"
#include "Storage/Flat/FlatEvent.hxx"
#include "Storage/Vector/VectorInjected.hxx"
#include "Storage/Vector/VectorMcParticles.hxx"
#include "Storage/Vector/VectorTracks.hxx"
#include "Storage/Vector/VectorV0s.hxx"
#include "View/MC/ViewMcParticle.hxx"
#ifdef T2S_LEGACY_KF
#include "Fit/Legacy/FitV0_Legacy.hxx"
#include "View/Reconstructed/Legacy/ViewTrack_Legacy.hxx"
#else
#include "Fit/FitV0.hxx"
#include "View/Reconstructed/ViewTrack.hxx"
#endif

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

    bool PrepareOutputFile();
    bool PrepareOutputTree();
    void CreateOutputBranches();
    void PrepareOutputHistograms();

    [[nodiscard]] int NumberEventsToRead() const {
        return fSettings.LimitToNEvents ? fSettings.LimitToNEvents : static_cast<int>(fInputChain_Events->GetEntries());
    }
    [[nodiscard]] bool IsMC() const { return fSettings.IsMC; }
    void GetEvent(int i_event) { fInputChain_Events->GetEntry(i_event); }

    [[nodiscard]] int NumberMC() const { return static_cast<int>(fInput_MC.Px->size()); }
    [[nodiscard]] int NumberInjected() const { return static_cast<int>(fInput_Injected.ReactionID->size()); }
    [[nodiscard]] int NumberTracks() const { return static_cast<int>(fInput_Tracks.Px->size()); }

    void ProcessEvent();

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
    void Store(const View::Rec::Track &view, Storage::Vector::Tracks &df);
    void StoreMC(const View::MC::Particle &view, Storage::Vector::MC_Tracks &df, EParticle pid);
    bool PassesCuts_Proton(const View::Rec::Track &track, TH1D *cut_flow_hist) const;
    bool PassesCuts_Kaon(const View::Rec::Track &track, TH1D *cut_flow_hist) const;
    bool PassesCuts_Pion(const View::Rec::Track &track, TH1D *cut_flow_hist) const;

    bool PassesCuts_Lambda(const Fit::V0 &v0, TH1D *cut_flow_hist) const;
    bool PassesCuts_KaonZeroShort(const Fit::V0 &v0, TH1D *cut_flow_hist) const;
    void Store(const Fit::V0 &v0, Storage::Vector::V0s &df);
    void StoreMC(const View::MC::V0 &v0_view, Storage::Vector::MC_V0s &df, EParticle v0_pid);

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

    // temporary containers, cleaned after event loop //

    std::vector<float> fVec_SV_X;
    std::vector<float> fVec_SV_Y;
    std::vector<float> fVec_SV_Z;

    std::vector<Fit::Track> fVec_AntiProtons;
    std::vector<Fit::Track> fVec_Protons;
    std::vector<Fit::Track> fVec_NegKaons;
    std::vector<Fit::Track> fVec_PosKaons;
    std::vector<Fit::Track> fVec_PiMinus;
    std::vector<Fit::Track> fVec_PiPlus;

    // input //

    Storage::Flat::Event fInput_Event;
    Storage::Flat::Coordinates fInput_MC_PV;
    Storage::Vector::Injected fInput_Injected;

    Storage::Vector::MCParticles fInput_MC;
    Storage::Vector::Tracks fInput_Tracks;

    // output //

    Storage::Flat::Event fOutput_Event;
    Storage::Flat::Coordinates fOutput_MC_PV;
    Storage::Vector::Injected fOutput_Injected;

    Storage::Vector::V0s fOutput_AntiLambdas;
    Storage::Vector::V0s fOutput_Lambdas;
    Storage::Vector::V0s fOutput_KaonsZeroShort;

    Storage::Vector::Tracks fOutput_NegKaons;
    Storage::Vector::Tracks fOutput_PosKaons;
    Storage::Vector::Tracks fOutput_PiMinus;
    Storage::Vector::Tracks fOutput_PiPlus;

    Storage::Vector::MC_V0s fOutput_MC_AntiLambdas;
    Storage::Vector::MC_V0s fOutput_MC_Lambdas;
    Storage::Vector::MC_V0s fOutput_MC_KaonsZeroShort;

    Storage::Vector::MC_Tracks fOutput_MC_NegKaons;
    Storage::Vector::MC_Tracks fOutput_MC_PosKaons;
    Storage::Vector::MC_Tracks fOutput_MC_PiMinus;
    Storage::Vector::MC_Tracks fOutput_MC_PiPlus;
};

}  // namespace Tree2Secondaries
