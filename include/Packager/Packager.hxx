#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include <TFile.h>
#include <TH1.h>

#include <Math/Point3D.h>

#include "App/Logger.hxx"
#include "common/Constants.hpp"
#include "common/DB_Particles.hpp"
#include "common/DB_ReactionChannels.hpp"
#include "common/Framework.hpp"
#include "common/Framework_TeeTree.hpp"
#include "common/Schema_Events.hpp"
#include "common/Schema_PackedEvents.hpp"

#include "App/Settings.hxx"

// forward declarations //
// clang-format off
namespace POD { struct Track; struct V0; }
namespace Cached { struct V0; }

namespace T2DS {

namespace Seeder{ struct PCA; }
namespace KF{ struct Particle; }
// clang-format on

// Pack secondary V0s and tracks.
class Packager {
   public:
    Packager(const Packager &) = delete;
    Packager(Packager &&) = delete;
    Packager &operator=(const Packager &) = delete;
    Packager &operator=(Packager &&) = delete;
    ~Packager() = default;

    explicit Packager(const Settings &settings)
        : fSettings{settings},
          fReactionChannel{DB::ReactionChannels::FindReactionChannel(fSettings.ReactionChannel)},
          // input
          fInput_File{std::make_unique<TFile>(fSettings.PathInputFile.c_str(), "READ")},
          fInput{},
          fReader{nullptr},
          // output
          fOutput_File{std::make_unique<TFile>(fSettings.PathOutputFile.c_str(), "RECREATE")},
          fOutput{},
          fWriter{nullptr} {

        fReader = std::make_unique<Framework::TeeTree::Reader>(fInput.CreateModel_TeeTree(fSettings.IsMC, fSettings.IsMC), E2T::Name_OutputTree,
                                                               *fInput_File);
        fWriter = std::make_unique<Framework::Writer>(fOutput.CreateModel(fSettings.IsMC), T2DS::Name_PackedRNT, *fOutput_File);

        PrepareOutputHistograms();

        Logger::Info(__FUNCTION__, "Packager initialized successfully.");
    }

    void PrepareOutputHistograms();

    void Load(long long entry_idx) { fReader->Load(entry_idx); }
    [[nodiscard]] long long NumberEventsToRead() {
        auto total = fReader->GetEntries();
        return fSettings.LimitToNEvents.has_value() ? std::min(static_cast<long long>(fSettings.LimitToNEvents.value()), total) : total;
    }

    void ProcessEvent();
    void ProcessInjectedSexa();
    void ProcessTracks();

    void Pack() {
        FindV0s(DB::Particles::Particle("AntiLambda"));
        FindV0s(DB::Particles::Particle("Lambda"));
        FindV0s(DB::Particles::Particle("KaonZeroShort"));
        PackTracks(DB::Particles::Particle("NegKaon"));
        PackTracks(DB::Particles::Particle("PosKaon"));
    }

    [[nodiscard]] bool FastCuts(const Seeder::PCA &pca_neg, const Seeder::PCA &pca_pos, const DB::Particles::Definition &pid) const {
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

    [[nodiscard]] bool SlowCuts(const Cached::V0 &c_v0, const DB::Particles::Definition &pid) const {
        switch (pid.pdg_code) {
            case DB::Particles::Particle("AntiLambda").pdg_code: {
                return SlowCuts_Lambda(c_v0, fHist_CutFlow_AntiLambda.get());
            }
            case DB::Particles::Particle("Lambda").pdg_code: {
                return SlowCuts_Lambda(c_v0, fHist_CutFlow_Lambda.get());
            }
            case DB::Particles::Particle("KaonZeroShort").pdg_code: {
                return SlowCuts_KaonZeroShort(c_v0, fHist_CutFlow_KaonZeroShort.get());
            }
            default:
                return false;
        }
    }

    void EndOfEvent();
    void EndOfAnalysis();

   private:
    // tracks //

    void PackTracks(const DB::Particles::Definition &pid);

    bool PassesProtonCuts(const POD::Track &track, TH1D *hist_cut_flow) const;
    bool PassesKaonCuts(const POD::Track &track, TH1D *hist_cut_flow) const;
    bool PassesPionCuts(const POD::Track &track, TH1D *hist_cut_flow) const;

    POD::Extended::McParticle BuildMcTrack(unsigned int track_mc_entry, int pdg_code_hypothesis, bool include_gm);

    // V0s //

    void FindV0s(const DB::Particles::Definition &pid);

    bool FastCuts_Lambda(const Seeder::PCA &pca_neg, const Seeder::PCA &pca_pos, TH1D *hist_cut_flow) const;
    bool FastCuts_KaonZeroShort(const Seeder::PCA &pca_neg, const Seeder::PCA &pca_pos, TH1D *hist_cut_flow) const;

    bool SlowCuts_Lambda(const Cached::V0 &v0, TH1D *hist_cut_flow) const;
    bool SlowCuts_KaonZeroShort(const Cached::V0 &v0, TH1D *hist_cut_flow) const;

    POD::Extended::McParticle BuildMcV0(const POD::Extended::McParticle &mc_neg, const POD::Extended::McParticle &mc_pos, int pdg_code_hypothesis);
    POD::V0 CreateV0(const KF::Particle &fit, const Seeder::PCA &neg_pca_wrt_v0, const Seeder::PCA &pos_pca_wrt_v0);

    // member variables //

    const Settings &fSettings;
    DB::ReactionChannels::Definition fReactionChannel;

    // temporary entries containers, cleaned after event loop //

    std::vector<std::size_t> fEntries_AntiProton;
    std::vector<std::size_t> fEntries_Proton;
    std::vector<std::size_t> fEntries_NegKaon;
    std::vector<std::size_t> fEntries_PosKaon;
    std::vector<std::size_t> fEntries_PiMinus;
    std::vector<std::size_t> fEntries_PiPlus;

    // input //

    std::unique_ptr<TFile> fInput_File;
    Schema::Events fInput;
    std::unique_ptr<Framework::TeeTree::Reader> fReader;
    // -- cached
    ROOT::Math::XYZPoint fPrimaryVertex;

    // output //

    std::unique_ptr<TFile> fOutput_File;
    Schema::PackedEvents fOutput;
    std::unique_ptr<Framework::Writer> fWriter;

    // histograms
    // -- event counter
    std::unique_ptr<TH1D> fHist_EventCounter;
    // -- cut flow for tracks
    std::unique_ptr<TH1D> fHist_CutFlow_AntiProton;
    std::unique_ptr<TH1D> fHist_CutFlow_Proton;
    std::unique_ptr<TH1D> fHist_CutFlow_NegKaon;
    std::unique_ptr<TH1D> fHist_CutFlow_PosKaon;
    std::unique_ptr<TH1D> fHist_CutFlow_PiMinus;
    std::unique_ptr<TH1D> fHist_CutFlow_PiPlus;
    // -- cut flow for v0s
    std::unique_ptr<TH1D> fHist_CutFlow_AntiLambda;
    std::unique_ptr<TH1D> fHist_CutFlow_Lambda;
    std::unique_ptr<TH1D> fHist_CutFlow_KaonZeroShort;
};

}  // namespace T2DS
