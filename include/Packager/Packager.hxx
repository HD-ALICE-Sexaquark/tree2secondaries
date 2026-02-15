#pragma once

#include <memory>

#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

#include "App/Settings.hxx"
#include "KalmanFitter/KalmanFitterV0.hxx"
#include "Math/Constants.hxx"
#include "Seeder/BaseSeeder.hxx"
#include "Storage/Flat/FlatEvent.hxx"
#include "Storage/Vector/VectorInjected.hxx"
#include "Storage/Vector/VectorMcParticles.hxx"
#include "Storage/Vector/VectorTracks.hxx"
#include "Storage/Vector/VectorV0s.hxx"
#include "View/MC/ViewMcParticle.hxx"
#include "View/Reconstructed/ViewRecTrack.hxx"

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

    [[nodiscard]] size_t NumberEventsToRead() const {
        return fSettings.LimitToNEvents ? fSettings.LimitToNEvents : static_cast<size_t>(fInputChain_Events->GetEntries());
    }
    [[nodiscard]] bool IsMC() const { return fSettings.IsMC; }
    void GetEvent(size_t i_event) { fInputChain_Events->GetEntry(static_cast<long long>(i_event)); }

    [[nodiscard]] size_t NumberMC() const { return fInput_MC.Px->size(); }
    [[nodiscard]] size_t NumberInjected() const { return fInput_Injected.ReactionID->size(); }
    [[nodiscard]] size_t NumberTracks() const { return fInput_Tracks.Px->size(); }

    void ProcessEvent();

    void Injected_GetSecondaryVertex();
    void Injected_Store();
    void ProcessInjected() {
        Injected_GetSecondaryVertex();
        Injected_Store();
    }

    void ProcessTracks();

    void PackTracks(PID_StableParticle pid);

    void FindV0s(PID_V0 pid);

    [[nodiscard]] bool FastCuts(const Seeder::Seed &pca_neg, const Seeder::Seed &pca_pos, PID_V0 pid) const {
        switch (pid) {
            case PID_V0::AntiLambda:
                return FastCuts_Lambda(pca_neg, pca_pos, fHist_CutFlow_AntiLambda.get());
            case PID_V0::Lambda:
                return FastCuts_Lambda(pca_neg, pca_pos, fHist_CutFlow_Lambda.get());
            case PID_V0::KaonZeroShort:
                return FastCuts_KaonZeroShort(pca_neg, pca_pos, fHist_CutFlow_KaonZeroShort.get());
            default:
                return false;
        }
    }

    [[nodiscard]] bool SlowCuts(const KalmanFitter::V0 &v0, PID_V0 pid) const {
        switch (pid) {
            case PID_V0::AntiLambda:
                return SlowCuts_Lambda(v0, fHist_CutFlow_AntiLambda.get());
            case PID_V0::Lambda:
                return SlowCuts_Lambda(v0, fHist_CutFlow_Lambda.get());
            case PID_V0::KaonZeroShort:
                return SlowCuts_KaonZeroShort(v0, fHist_CutFlow_KaonZeroShort.get());
            default:
                return false;
        }
    }

    void EndOfEvent();
    void EndOfAnalysis();

   private:
    void Store(const View::Rec::Track &track, Storage::Vector::Tracks &df);
    void StoreMC(const View::MC::Particle &mc, Storage::Vector::MC_Tracks &df, PID_StableParticle pid);
    void StoreDummyMC(Storage::Vector::MC_Tracks &df);

    bool Cuts_Proton(const View::Rec::Track &track, TH1D *cut_flow_hist) const;
    bool Cuts_Kaon(const View::Rec::Track &track, TH1D *cut_flow_hist) const;
    bool Cuts_Pion(const View::Rec::Track &track, TH1D *cut_flow_hist) const;

    bool FastCuts_Lambda(const Seeder::Seed &seed_neg, const Seeder::Seed &seed_pos, TH1D *cut_flow_hist) const;
    bool FastCuts_KaonZeroShort(const Seeder::Seed &seed_neg, const Seeder::Seed &seed_pos, TH1D *cut_flow_hist) const;

    bool SlowCuts_Lambda(const KalmanFitter::V0 &v0, TH1D *cut_flow_hist) const;
    bool SlowCuts_KaonZeroShort(const KalmanFitter::V0 &v0, TH1D *cut_flow_hist) const;

    void Store(const KalmanFitter::V0 &v0, Storage::Vector::V0s &df);
    void StoreMC(const View::MC::V0 &v0, Storage::Vector::MC_V0s &df, PID_V0 pid);
    void StoreDummyMC(Storage::Vector::MC_V0s &df);

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

    std::vector<size_t> fIndices_AntiProtons;
    std::vector<size_t> fIndices_Protons;
    std::vector<size_t> fIndices_NegKaons;
    std::vector<size_t> fIndices_PosKaons;
    std::vector<size_t> fIndices_PiMinus;
    std::vector<size_t> fIndices_PiPlus;

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
