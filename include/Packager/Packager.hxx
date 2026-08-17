#pragma once

#include <cstddef>
#include <exception>
#include <memory>
#include <string>
#include <string_view>
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
namespace KF{ struct Particle; struct FitResult; }
// clang-format on

// Pack secondary V0s and tracks.
class Packager {
    // PENDING for author: missing fit constraints + cuts enums+structs

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
          fInput{},
          fReader{nullptr},
          // output
          fOutput_File{std::make_unique<TFile>(fSettings.PathOutputFile.c_str(), "RECREATE")},
          fOutput{},
          fWriter{nullptr} {

        fWriter = std::make_unique<Framework::Writer>(fOutput.CreateModel(fSettings.IsMC), T2DS::Name_PackedRNT, *fOutput_File);

        PrepareOutputHistograms();

        Logger::Info(__FUNCTION__, "Packager initialized successfully.");
    }

    [[nodiscard]] bool OpenInput(std::string_view path) {
        fReader.reset();  // raw `TTree*` must die first
        try {
            fReader =
                std::make_unique<Framework::TeeTree::Reader>(fInput.CreateModel_TeeTree(fSettings.IsMC, fSettings.IsMC), E2T::Name_OutputTree, path);
        } catch (const std::exception &exc) {
            Logger::Error(__FUNCTION__, "Couldn't read {} ({}) -- skipping it.", path, exc.what());
            return false;
        }
        return true;
    }

    void PrepareOutputHistograms();

    void Load(long long entry_idx) { fReader->Load(entry_idx); }
    [[nodiscard]] unsigned long NumberEventsToRead() { return static_cast<unsigned long>(fReader->GetEntries()); }

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
    bool EndOfAnalysis();

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
    POD::V0 CreateV0(const KF::FitResult &fit, const Seeder::PCA &neg_pca_wrt_v0, const Seeder::PCA &pos_pca_wrt_v0);

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

    Schema::Events fInput;
    std::unique_ptr<Framework::TeeTree::Reader> fReader;
    // -- cached
    ROOT::Math::XYZPoint fPrimaryVertex;

    // output //

    std::unique_ptr<TFile> fOutput_File;  // single file, kept alive across every input file, if multiple
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
