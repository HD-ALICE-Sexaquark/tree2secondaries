#pragma once

#include <cstddef>
#include <memory>

#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

#include "App/Settings.hxx"
#include "KalmanFitter/KalmanFitterV0.hxx"
#include "Math/Constants.hxx"
#include "Seeder/BaseSeeder.hxx"
#include "Storage/Schema/SchemaFlat.hxx"
#include "Storage/Schema/SchemaVector.hxx"
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

    [[nodiscard]] std::size_t NumberMC() const { return fInput_MC.lv.px->size(); }
    [[nodiscard]] std::size_t NumberInjected() const { return fInput_Injected.reaction_id->size(); }
    [[nodiscard]] std::size_t NumberTracks() const { return fInput_Tracks.state.px->size(); }

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
    void Store(const View::Rec::Track &track, Schema::Vector::Tracks &df);
    void StoreMC(const View::MC::Particle &mc, Schema::Vector::MC_Tracks &df, PID_StableParticle pid);
    void StoreDummyMC(Schema::Vector::MC_Tracks &df);

    bool Cuts_Proton(const View::Rec::Track &track, TH1D *cut_flow_hist) const;
    bool Cuts_Kaon(const View::Rec::Track &track, TH1D *cut_flow_hist) const;
    bool Cuts_Pion(const View::Rec::Track &track, TH1D *cut_flow_hist) const;

    bool FastCuts_Lambda(const Seeder::Seed &seed_neg, const Seeder::Seed &seed_pos, TH1D *cut_flow_hist) const;
    bool FastCuts_KaonZeroShort(const Seeder::Seed &seed_neg, const Seeder::Seed &seed_pos, TH1D *cut_flow_hist) const;

    bool SlowCuts_Lambda(const KalmanFitter::V0 &v0, TH1D *cut_flow_hist) const;
    bool SlowCuts_KaonZeroShort(const KalmanFitter::V0 &v0, TH1D *cut_flow_hist) const;

    void Store(const KalmanFitter::V0 &v0, Schema::Vector::V0s &df);
    void StoreMC(const View::MC::V0 &mc_v0, Schema::Vector::MC_V0s &df, PID_V0 pid);
    void StoreDummyMC(Schema::Vector::MC_V0s &df);

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

    std::vector<unsigned int> fIndices_AntiProtons;
    std::vector<unsigned int> fIndices_Protons;
    std::vector<unsigned int> fIndices_NegKaons;
    std::vector<unsigned int> fIndices_PosKaons;
    std::vector<unsigned int> fIndices_PiMinus;
    std::vector<unsigned int> fIndices_PiPlus;

    // input //

    Schema::Flat::Event fInput_Event;
    Schema::Vector::Injected fInput_Injected;
    Schema::Vector::MCParticles fInput_MC;
    Schema::Vector::Tracks fInput_Tracks;

    // output //

    Schema::Flat::Event fOutput_Event;
    Schema::Vector::Injected fOutput_Injected;

    Schema::Vector::V0s fOutput_AntiLambdas;
    Schema::Vector::V0s fOutput_Lambdas;
    Schema::Vector::V0s fOutput_KaonsZeroShort;

    Schema::Vector::Tracks fOutput_NegKaons;
    Schema::Vector::Tracks fOutput_PosKaons;
    Schema::Vector::Tracks fOutput_PiMinus;
    Schema::Vector::Tracks fOutput_PiPlus;

    Schema::Vector::MC_V0s fOutput_MC_AntiLambdas;
    Schema::Vector::MC_V0s fOutput_MC_Lambdas;
    Schema::Vector::MC_V0s fOutput_MC_KaonsZeroShort;

    Schema::Vector::MC_Tracks fOutput_MC_NegKaons;
    Schema::Vector::MC_Tracks fOutput_MC_PosKaons;
    Schema::Vector::MC_Tracks fOutput_MC_PiMinus;
    Schema::Vector::MC_Tracks fOutput_MC_PiPlus;
};

}  // namespace Tree2Secondaries
