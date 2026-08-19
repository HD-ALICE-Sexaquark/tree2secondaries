#pragma once

#include <exception>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include <TFile.h>
#include <TH1.h>

#include <Math/Point3D.h>

#include "common/Constants.hpp"
#include "common/DB_Particles.hpp"
#include "common/DB_ReactionChannels.hpp"
#include "common/Framework.hpp"
#include "common/Framework_TeeTree.hpp"
#include "common/Schema_Events.hpp"
#include "common/Schema_FoundSexaquark.hpp"

#include "App/Logger.hxx"
#include "App/Settings.hxx"

// forward declarations //
// clang-format off
namespace Cached { struct Sexaquark; struct V0; }

namespace T2DS {

namespace Seeder{ struct PCA; }
namespace KF{ struct Particle; struct FitResult; }
// clang-format on

// Collect secondary charged kaons and reconstruct V0s, then combine them into anti-sexaquark candidates.
// The secondaries are built once per event and shared by every requested reaction channel.
class Finder {
    // PENDING for author: missing fit constraints + cuts enums+structs + duplications treatments
    // Cut-related //

    enum class EProton : int {
        kAllPossibleProtons,
        // PENDING
        kPasses_NSigmasProtons,
        // --
        kNProtonCuts,
    };
    enum class EKaon : int {
        kAllPossibleKaons,
        // PENDING
        kPasses_NSigmasKaons,
        // --
        kNKaonCuts,
    };
    enum class EPion : int {
        kAllPossiblePions,
        // PENDING
        kPasses_NSigmasPions,
        // --
        kNPionCuts,
    };
    enum class ELambda : int {
        kAllCombinations,
        // PENDING
        kPasses_DcaBtwDaughters,
        // --
        kNLambdaCuts,
    };
    enum class EKaonZeroShort : int {
        kAllCombinations,
        // PENDING
        kPasses_DcaBtwDaughters,
        // --
        kNKaonZeroShortCuts,
    };
    enum class EChannelA : int {
        kAllCombinations,
        // PENDING
        kNChannelACuts,
    };
    enum class EChannelD : int {
        kAllCombinations,
        // PENDING
        kNChannelDCuts,
    };
    enum class EChannelH : int {
        kAllCombinations,
        // PENDING
        kNChannelHCuts,
    };
    template <typename E>
    static void FillHist(TH1D *hist, const E &bin_n) {
        hist->Fill(static_cast<double>(bin_n));
    }

    // Injected-related //

    struct InjectedReaction {
        double after_px{0.};
        double after_py{0.};
        double after_pz{0.};
        double after_e{0.};
        float sv_x{Common::DummyFloat};
        float sv_y{Common::DummyFloat};
        float sv_z{Common::DummyFloat};
        bool found{false};
    };

   public:
    Finder() = delete;
    Finder(const Finder &) = delete;
    Finder(Finder &&) = delete;
    Finder &operator=(const Finder &) = delete;
    Finder &operator=(Finder &&) = delete;
    ~Finder() = default;

    explicit Finder(const Settings &settings)
        : fSettings{settings},
          // input
          fInput{},
          fReader{nullptr},
          // output
          fOutput_File{std::make_unique<TFile>(fSettings.PathOutputFile.c_str(), "RECREATE")},
          fOutput{},
          fWriter{std::make_unique<Framework::Writer>(fOutput.CreateModel(fSettings.IsMC), T2DS::Name_FoundSexaquarkRNT, *fOutput_File)} {

        PrepareOutputHistograms();

        Logger::Info(__FUNCTION__, "Finder initialized successfully.");
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

    void Load(long long entry_idx) { fReader->Load(entry_idx); }
    [[nodiscard]] unsigned long NumberEventsToRead() { return static_cast<unsigned long>(fReader->GetEntries()); }

    void ProcessEvent();
    void ProcessInjected();
    void ProcessTracks();

    void FindV0s() {
        FindV0s(DB::Particles::Particle("AntiLambda"));
        FindV0s(DB::Particles::Particle("Lambda"));
        FindV0s(DB::Particles::Particle("KaonZeroShort"));
    }

    void FindSexaquarks() {
        FindSexaquarks_ChannelA(false);
        FindSexaquarks_ChannelA(true);
        FindSexaquarks_ChannelD(false);
        FindSexaquarks_ChannelD(true);
        FindSexaquarks_ChannelH(false);
        FindSexaquarks_ChannelH(true);
    }

    void EndOfEvent();
    bool EndOfAnalysis();

   private:
    void PrepareOutputHistograms();

    [[nodiscard]] bool PostSeedCuts(const Seeder::PCA &pca_neg, const Seeder::PCA &pca_pos, const DB::Particles::Definition &pid) {
        switch (pid.pdg_code) {
            case DB::Particles::Particle("AntiLambda").pdg_code: {
                return PostSeedCuts_Lambda(pca_neg, pca_pos, fHist_CutFlow_AntiLambda.get());
            }
            case DB::Particles::Particle("Lambda").pdg_code: {
                return PostSeedCuts_Lambda(pca_neg, pca_pos, fHist_CutFlow_Lambda.get());
            }
            case DB::Particles::Particle("KaonZeroShort").pdg_code: {
                return PostSeedCuts_KaonZeroShort(pca_neg, pca_pos, fHist_CutFlow_KaonZeroShort.get());
            }
            default:
                return false;
        }
    }

    [[nodiscard]] bool PostFitCuts(const Cached::V0 &c_v0, const DB::Particles::Definition &pid) {
        switch (pid.pdg_code) {
            case DB::Particles::Particle("AntiLambda").pdg_code: {
                return PostFitCuts_Lambda(c_v0, fHist_CutFlow_AntiLambda.get());
            }
            case DB::Particles::Particle("Lambda").pdg_code: {
                return PostFitCuts_Lambda(c_v0, fHist_CutFlow_Lambda.get());
            }
            case DB::Particles::Particle("KaonZeroShort").pdg_code: {
                return PostFitCuts_KaonZeroShort(c_v0, fHist_CutFlow_KaonZeroShort.get());
            }
            default:
                return false;
        }
    }

    // tracks //

    bool PassesCuts_Proton(const POD::Track &track, TH1D *hist_cut_flow) const;
    bool PassesCuts_Kaon(const POD::Track &track, TH1D *hist_cut_flow) const;
    bool PassesCuts_Pion(const POD::Track &track, TH1D *hist_cut_flow) const;

    POD::Extended::McParticle BuildMcTrack(unsigned int track_mc_entry, int pdg_code_hypothesis, bool include_gm);

    // V0s //

    void FindV0s(const DB::Particles::Definition &pid);

    bool PreSeedCuts_Lambda() const;         // PENDING
    bool PreSeedCuts_KaonZeroShort() const;  // PENDING

    bool PostSeedCuts_Lambda(const Seeder::PCA &pca_neg, const Seeder::PCA &pca_pos, TH1D *hist_cut_flow) const;
    bool PostSeedCuts_KaonZeroShort(const Seeder::PCA &pca_neg, const Seeder::PCA &pca_pos, TH1D *hist_cut_flow) const;

    bool PostFitCuts_Lambda(const Cached::V0 &v0, TH1D *hist_cut_flow) const;
    bool PostFitCuts_KaonZeroShort(const Cached::V0 &v0, TH1D *hist_cut_flow) const;

    POD::Extended::McParticle BuildMcV0(const POD::Extended::McParticle &mc_neg, const POD::Extended::McParticle &mc_pos, int pdg_code_hypothesis);
    POD::V0 CreateV0(const KF::FitResult &fit, const Seeder::PCA &neg_pca_wrt_v0, const Seeder::PCA &pos_pca_wrt_v0);

    // mc sexaquark //
    POD::Linked::InjectedSexa BuildMcSexaquark(const POD::Extended::McParticle &mc_dau1, const POD::Extended::McParticle &mc_dau2);

    // channel A //
    void FindSexaquarks_ChannelA(bool is_bkg_channel);
    bool PreSeedCuts_ChannelA() const;  // PENDING
    [[nodiscard]] bool PostSeedCuts_ChannelA(const Seeder::PCA &pca_v0a, const Seeder::PCA &pca_v0b, TH1D *hist_cut_flow) const;
    [[nodiscard]] bool PostFitCuts_ChannelA(const Cached::Sexaquark &c_sexa, TH1D *hist_cut_flow) const;
    POD::Sexaquark Create_ChannelA(const KF::FitResult &fit, const Seeder::PCA &pca_v0a, const Seeder::PCA &pca_v0b, bool is_bkg_channel);

    // channel D //
    void FindSexaquarks_ChannelD(bool is_bkg_channel);
    bool PreSeedCuts_ChannelD() const;  // PENDING
    [[nodiscard]] bool PostSeedCuts_ChannelD(const Seeder::PCA &pca_v0, const Seeder::PCA &pca_ka, TH1D *hist_cut_flow) const;
    [[nodiscard]] bool PostFitCuts_ChannelD(const Cached::Sexaquark &c_sexa, TH1D *hist_cut_flow) const;
    POD::Sexaquark Create_ChannelD(const KF::FitResult &fit, const Seeder::PCA &pca_v0, const Seeder::PCA &pca_ka, bool is_bkg_channel);

    // channel H //
    void FindSexaquarks_ChannelH(bool is_bkg_channel);
    bool PreSeedCuts_ChannelH() const;  // PENDING
    [[nodiscard]] bool PostSeedCuts_ChannelH(const Seeder::PCA &pca_kaon1, const Seeder::PCA &pca_kaon2, TH1D *hist_cut_flow) const;
    [[nodiscard]] bool PostFitCuts_ChannelH(const Cached::Sexaquark &c_sexa, TH1D *hist_cut_flow) const;
    POD::Sexaquark Create_ChannelH(const KF::FitResult &fit, const Seeder::PCA &pca_kaon1, const Seeder::PCA &pca_kaon2, bool is_bkg_channel);

    // member variables //

    const Settings &fSettings;

    // temporary in-memory data //

    std::vector<POD::Track> fTemp_AntiProton;
    std::vector<POD::Track> fTemp_Proton;
    std::vector<POD::Track> fTemp_NegKaon;
    std::vector<POD::Track> fTemp_PosKaon;
    std::vector<POD::Track> fTemp_PiMinus;
    std::vector<POD::Track> fTemp_PiPlus;

    std::vector<POD::Extended::McParticle> fTemp_MC_AntiProton;
    std::vector<POD::Extended::McParticle> fTemp_MC_Proton;
    std::vector<POD::Extended::McParticle> fTemp_MC_NegKaon;
    std::vector<POD::Extended::McParticle> fTemp_MC_PosKaon;
    std::vector<POD::Extended::McParticle> fTemp_MC_PiMinus;
    std::vector<POD::Extended::McParticle> fTemp_MC_PiPlus;

    std::vector<POD::V0> fTemp_AntiLambda;
    std::vector<POD::Track> fTemp_AntiLambda_Neg;
    std::vector<POD::Track> fTemp_AntiLambda_Pos;
    std::vector<POD::V0> fTemp_Lambda;
    std::vector<POD::Track> fTemp_Lambda_Neg;
    std::vector<POD::Track> fTemp_Lambda_Pos;
    std::vector<POD::V0> fTemp_KaonZeroShort;
    std::vector<POD::Track> fTemp_KaonZeroShort_Neg;
    std::vector<POD::Track> fTemp_KaonZeroShort_Pos;

    std::vector<POD::Extended::McParticle> fTemp_MC_AntiLambda;
    std::vector<POD::Extended::McParticle> fTemp_MC_AntiLambda_Neg;
    std::vector<POD::Extended::McParticle> fTemp_MC_AntiLambda_Pos;
    std::vector<POD::Extended::McParticle> fTemp_MC_Lambda;
    std::vector<POD::Extended::McParticle> fTemp_MC_Lambda_Neg;
    std::vector<POD::Extended::McParticle> fTemp_MC_Lambda_Pos;
    std::vector<POD::Extended::McParticle> fTemp_MC_KaonZeroShort;
    std::vector<POD::Extended::McParticle> fTemp_MC_KaonZeroShort_Neg;
    std::vector<POD::Extended::McParticle> fTemp_MC_KaonZeroShort_Pos;

    // input //

    Schema::Events fInput;
    std::unique_ptr<Framework::TeeTree::Reader> fReader;
    // -- cached
    ROOT::Math::XYZPoint fPrimaryVertex;
    double fMagneticField{0.};
    // -- injected reaction channel of dedicated sexa mc production, identified on the first event
    //    without value until detected injected channel; value '0' if there are no injected particles
    std::optional<DB::ReactionChannels::Definition> fMcSignalChannel;

    // output //

    std::unique_ptr<TFile> fOutput_File;  // single file, kept alive across every input file, if multiple
    Schema::FoundSexaquark fOutput;
    std::unique_ptr<Framework::Writer> fWriter;

    // histograms //

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
    // -- cut flow for anti-sexaquarks + bkg. sexaquarks
    std::unique_ptr<TH1D> fHist_CutFlow_ChannelA;
    std::unique_ptr<TH1D> fHist_CutFlow_ChannelA_Bkg;
    std::unique_ptr<TH1D> fHist_CutFlow_ChannelD;
    std::unique_ptr<TH1D> fHist_CutFlow_ChannelD_Bkg;
    std::unique_ptr<TH1D> fHist_CutFlow_ChannelH;
    std::unique_ptr<TH1D> fHist_CutFlow_ChannelH_Bkg;
};

}  // namespace T2DS
