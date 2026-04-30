#pragma once

#include <cstddef>
#include <memory>
#include <optional>
#include <vector>

#include <TChain.h>
#include <TFile.h>
#include <TH1.h>

#include "App/DB_Particles.hxx"
#include "App/Settings.hxx"
#include "KalmanFitter/KalmanFitterV0.hxx"
#include "MonteCarlo/MonteCarloParticle.hxx"
#include "MonteCarlo/MonteCarloV0.hxx"
#include "Schema/SchemaFlatEvent.hxx"
#include "Schema/SchemaVectorInjected.hxx"
#include "Schema/SchemaVectorMcParticles.hxx"
#include "Schema/SchemaVectorMcTracks.hxx"
#include "Schema/SchemaVectorMcV0s.hxx"
#include "Schema/SchemaVectorTracks.hxx"
#include "Schema/SchemaVectorV0s.hxx"
#include "Seeder/BaseSeeder.hxx"
#include "View/ViewVectorTracks.hxx"

class TTree;

namespace Tree2Secondaries {

// Pack secondary V0s and tracks.
class Packager {
   public:
    Packager(const Packager &) = delete;
    Packager(Packager &&) = delete;
    Packager &operator=(const Packager &) = delete;
    Packager &operator=(Packager &&) = delete;
    ~Packager() = default;

    explicit Packager(const Settings &settings) : fSettings{settings} {}

    bool Initialize();
    void ReadInputBranches();

    bool PrepareOutputFile();
    bool PrepareOutputTree();
    void CreateOutputBranches();
    void PrepareOutputHistograms();

    [[nodiscard]] long long NumberEventsToRead() const {
        long long total_entries = fInputChain_Events->GetEntries();
        return 0 < fSettings.LimitToNEvents && fSettings.LimitToNEvents < total_entries ? fSettings.LimitToNEvents : total_entries;
    }
    void GetEvent(long long i_event) { fInputChain_Events->GetEntry(i_event); }

    [[nodiscard]] std::size_t NumberMC() const { return fInput_MC.lv.px->size(); }
    [[nodiscard]] std::size_t NumberInjected() const { return fInput_Injected.reaction_id->size(); }
    [[nodiscard]] std::size_t NumberTracks() const { return fInput_Tracks.state.px->size(); }

    void ProcessEvent();

    void ProcessInjected();

    void ProcessTracks();
    void PackTracks(const Particles::Definition &pid);
    void FindV0s(const Particles::Definition &pid);

    void Pack() {
        FindV0s(Particles::Particle("AntiLambda"));
        FindV0s(Particles::Particle("Lambda"));
        FindV0s(Particles::Particle("KaonZeroShort"));
        PackTracks(Particles::Particle("NegKaon"));
        PackTracks(Particles::Particle("PosKaon"));
    }

    [[nodiscard]] bool FastCuts(const Seeder::Seed &pca_neg, const Seeder::Seed &pca_pos, const Particles::Definition &pid) const {
        switch (pid.pdg_code) {
            case Particles::Particle("AntiLambda").pdg_code:
                return FastCuts_SecondaryLambdas(pca_neg, pca_pos, fHist_CutFlow_AntiLambda.get());
            case Particles::Particle("Lambda").pdg_code:
                return FastCuts_SecondaryLambdas(pca_neg, pca_pos, fHist_CutFlow_Lambda.get());
            case Particles::Particle("KaonZeroShort").pdg_code:
                return FastCuts_KaonZeroShort(pca_neg, pca_pos, fHist_CutFlow_KaonZeroShort.get());
            default:
                return false;
        }
    }

    [[nodiscard]] bool SlowCuts(const KalmanFitter::V0 &v0, const Particles::Definition &pid) const {
        switch (pid.pdg_code) {
            case Particles::Particle("AntiLambda").pdg_code:
                return SlowCuts_SecondaryLambdas(v0, fHist_CutFlow_AntiLambda.get());
            case Particles::Particle("Lambda").pdg_code:
                return SlowCuts_SecondaryLambdas(v0, fHist_CutFlow_Lambda.get());
            case Particles::Particle("KaonZeroShort").pdg_code:
                return SlowCuts_KaonZeroShort(v0, fHist_CutFlow_KaonZeroShort.get());
            default:
                return false;
        }
    }

    void EndOfEvent();
    void EndOfAnalysis();

   private:
    void Store(const View::VecTracks &v, Schema::VecTracks &df);
    void Store(const std::optional<MonteCarlo::Particle> &mc, Schema::VecMcTracks &df, bool ascendants_info);

    bool PassesProtonCuts(const View::VecTracks &v, TH1D *cut_flow_hist) const;
    bool PassesKaonCuts(const View::VecTracks &v, TH1D *cut_flow_hist) const;
    bool PassesPionCuts(const View::VecTracks &v, TH1D *cut_flow_hist) const;

    bool FastCuts_SecondaryLambdas(const Seeder::Seed &seed_neg, const Seeder::Seed &seed_pos, TH1D *cut_flow_hist) const;
    bool FastCuts_KaonZeroShort(const Seeder::Seed &seed_neg, const Seeder::Seed &seed_pos, TH1D *cut_flow_hist) const;

    bool SlowCuts_SecondaryLambdas(const KalmanFitter::V0 &v0, TH1D *cut_flow_hist) const;
    bool SlowCuts_KaonZeroShort(const KalmanFitter::V0 &v0, TH1D *cut_flow_hist) const;

    void Store(const KalmanFitter::V0 &v0, Schema::VecV0s &df);
    void Store(const std::optional<MonteCarlo::V0> &mc_v0, Schema::VecMcV0s &df);

    const Settings &fSettings;
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

    // temporary indices containers, cleaned after event loop //

    std::vector<std::size_t> fIndices_AntiProtons;
    std::vector<std::size_t> fIndices_Protons;
    std::vector<std::size_t> fIndices_NegKaons;
    std::vector<std::size_t> fIndices_PosKaons;
    std::vector<std::size_t> fIndices_PiMinus;
    std::vector<std::size_t> fIndices_PiPlus;

    // input //

    Schema::FlatEvent fInput_Event;
    Schema::VecInjected fInput_Injected;
    Schema::VecMcParticles fInput_MC;
    Schema::VecTracks fInput_Tracks;

    // output //

    Schema::FlatEvent fOutput_Event;
    Schema::VecInjected fOutput_Injected;

    Schema::VecV0s fOutput_AntiLambdas;
    Schema::VecV0s fOutput_Lambdas;
    Schema::VecV0s fOutput_KaonsZeroShort;

    Schema::VecTracks fOutput_NegKaons;
    Schema::VecTracks fOutput_PosKaons;

    Schema::VecMcV0s fOutput_MC_AntiLambdas;
    Schema::VecMcTracks fOutput_MC_AntiLambdas_Neg;
    Schema::VecMcTracks fOutput_MC_AntiLambdas_Pos;
    Schema::VecMcV0s fOutput_MC_Lambdas;
    Schema::VecMcTracks fOutput_MC_Lambdas_Neg;
    Schema::VecMcTracks fOutput_MC_Lambdas_Pos;
    Schema::VecMcV0s fOutput_MC_KaonsZeroShort;
    Schema::VecMcTracks fOutput_MC_KaonsZeroShort_Neg;
    Schema::VecMcTracks fOutput_MC_KaonsZeroShort_Pos;

    Schema::VecMcTracks fOutput_MC_NegKaons;
    Schema::VecMcTracks fOutput_MC_PosKaons;
};

}  // namespace Tree2Secondaries
