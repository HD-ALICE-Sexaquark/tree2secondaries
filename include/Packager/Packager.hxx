#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include <TFile.h>
#include <TH1.h>

#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleWriter.hxx>

#include "common/DB_Particles.hpp"
#include "common/FL_Event.hpp"
#include "common/VC_InjectedSexa.hpp"
#include "common/VC_McParticle.hpp"
#include "common/VC_McParticleView.hpp"
#include "common/VC_Track.hpp"
#include "common/VC_TrackView.hpp"
#include "common/VC_V0.hpp"

#include "App/Settings.hxx"
#include "KalmanFitter/KalmanFitterV0.hxx"
#include "Seeder/BaseSeeder.hxx"

namespace R2DS {

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

    bool PrepareOutputFile();
    void PrepareOutputHistograms();

    [[nodiscard]] unsigned long NumberEventsToRead() const {
        unsigned long total_entries = fReader->GetNEntries();
        return 0 < fSettings.LimitToNEvents && fSettings.LimitToNEvents < total_entries ? fSettings.LimitToNEvents : total_entries;
    }

    void ProcessEvent();
    void ProcessInjected();
    void ProcessTracks();

    void Pack() {
        FindV0s(DB::Particles::Particle("AntiLambda"));
        FindV0s(DB::Particles::Particle("Lambda"));
        FindV0s(DB::Particles::Particle("KaonZeroShort"));
        PackTracks(DB::Particles::Particle("NegKaon"));
        PackTracks(DB::Particles::Particle("PosKaon"));
    }

    [[nodiscard]] bool FastCuts(const Seeder::Seed &pca_neg, const Seeder::Seed &pca_pos, const DB::Particles::Definition &pid) const {
        switch (pid.pdg_code) {
            case DB::Particles::Particle("AntiLambda").pdg_code: {
                return FastCuts_Lambda(pca_neg, pca_pos, fHist_CutFlow_AntiLambda.get());
            }
            case DB::Particles::Particle("Lambda").pdg_code: {
                return FastCuts_Lambda(pca_neg, pca_pos, fHist_CutFlow_Lambda.get());
            }
            case DB::Particles::Particle("KaonZeroShort").pdg_code: {
                return FastCuts_KaonZeroShort(pca_neg, pca_pos, fHist_CutFlow_KaonZeroShort.get());
            }
            default:
                return false;
        }
    }

    [[nodiscard]] bool SlowCuts(const KalmanFitter::V0 &v0, const DB::Particles::Definition &pid) const {
        switch (pid.pdg_code) {
            case DB::Particles::Particle("AntiLambda").pdg_code: {
                return SlowCuts_Lambda(v0, fHist_CutFlow_AntiLambda.get());
            }
            case DB::Particles::Particle("Lambda").pdg_code: {
                return SlowCuts_Lambda(v0, fHist_CutFlow_Lambda.get());
            }
            case DB::Particles::Particle("KaonZeroShort").pdg_code: {
                return SlowCuts_KaonZeroShort(v0, fHist_CutFlow_KaonZeroShort.get());
            }
            default:
                return false;
        }
    }

    void EndOfEvent();
    void EndOfAnalysis();

    std::unique_ptr<ROOT::RNTupleModel> fInput_Model;
    std::unique_ptr<ROOT::RNTupleReader> fReader;

    std::unique_ptr<ROOT::RNTupleModel> fOutput_Model;
    std::unique_ptr<ROOT::RNTupleWriter> fWriter;

   private:
    void PackTracks(const DB::Particles::Definition &pid);
    void FindV0s(const DB::Particles::Definition &pid);

    void Store(const Vector::TrackView *track, const Vector::McParticleView *linked_mc, Vector::Track &df);

    bool PassesProtonCuts(const Vector::TrackView &track, TH1D *cut_flow_hist) const;
    bool PassesKaonCuts(const Vector::TrackView &track, TH1D *cut_flow_hist) const;
    bool PassesPionCuts(const Vector::TrackView &track, TH1D *cut_flow_hist) const;

    bool FastCuts_Lambda(const Seeder::Seed &seed_neg, const Seeder::Seed &seed_pos, TH1D *cut_flow_hist) const;
    bool FastCuts_KaonZeroShort(const Seeder::Seed &seed_neg, const Seeder::Seed &seed_pos, TH1D *cut_flow_hist) const;

    bool SlowCuts_Lambda(const KalmanFitter::V0 &v0, TH1D *cut_flow_hist) const;
    bool SlowCuts_KaonZeroShort(const KalmanFitter::V0 &v0, TH1D *cut_flow_hist) const;

    void Store(const KalmanFitter::V0 *v0, const Vector::McParticleView *mc_neg, const Vector::McParticleView *mc_pos,
               const Vector::McParticleView *mc_v0, Vector::V0 &df);

    const Settings &fSettings;

    std::unique_ptr<TFile> fOutput_File;

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

    // temporary entries containers, cleaned after event loop //

    std::vector<std::size_t> fEntries_AntiProton;
    std::vector<std::size_t> fEntries_Proton;
    std::vector<std::size_t> fEntries_NegKaon;
    std::vector<std::size_t> fEntries_PosKaon;
    std::vector<std::size_t> fEntries_PiMinus;
    std::vector<std::size_t> fEntries_PiPlus;

    // input //

    Flat::Event fInput_Event;
    ROOT::Math::XYZPoint fPrimaryVertex;

    Vector::InjectedSexa fInput_InjectedSexa;
    Vector::McParticle fInput_McParticle;
    Vector::Track fInput_Track;

    // output //

    Flat::Event fOutput_Event;

    Vector::InjectedSexa fOutput_InjectedSexa;

    Vector::V0 fOutput_AntiLambda;
    Vector::V0 fOutput_Lambda;
    Vector::V0 fOutput_KaonZeroShort;

    Vector::Track fOutput_NegKaon;
    Vector::Track fOutput_PosKaon;
};

}  // namespace R2DS
