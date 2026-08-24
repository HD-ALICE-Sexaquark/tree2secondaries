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
#include "KalmanFitter/BaseKalmanFitter.hxx"

// forward declarations //
// clang-format off
namespace Cached { struct ChannelA; struct ChannelD; struct ChannelH; struct V0; }

namespace T2DS {

namespace Seeder{ struct PCA; }
// clang-format on

// Collect secondary charged kaons and reconstruct V0s, then combine them into antisexaquark candidates.
// The secondaries are built once per event and shared by every requested reaction channel.
class Finder {
    // Cut-related //

    enum class EProton : int {
        kAllPossibleProtons,
        kPasses_NSigmasProtons,  // PENDING
        kNProtonCuts,
    };
    enum class EKaon : int {
        kAllPossibleKaons,
        kPasses_NSigmasKaons,  // PENDING
        kNKaonCuts,
    };
    enum class EPion : int {
        kAllPossiblePions,
        kPasses_NSigmasPions,  // PENDING
        kNPionCuts,
    };
    enum class ELambda : int {
        kAllCombinations,
        // pre-seed cuts
        kPasses_DiffDaughters_Logical,
        kPasses_DiffDaughters_Physical,
        // post-seed cuts
        // PENDING
        kPasses_DcaBtwDaughters,
        // post-fit cuts
        // PENDING
        // --
        kNLambdaCuts,
    };
    enum class EKaonZeroShort : int {
        kAllCombinations,
        // pre-seed cuts
        kPasses_DiffDaughters_Logical,
        kPasses_DiffDaughters_Physical,
        // post-seed cuts
        // PENDING
        kPasses_DcaBtwDaughters,
        // post-fit cuts
        // PENDING
        // --
        kNKaonZeroShortCuts,
    };
    enum class EChannelA : int {
        kAllCombinations,
        // pre-seed cuts
        kPasses_DiffTracks_Logical,
        kPasses_DiffTracks_Physical,
        kPasses_DiffV0s_Physical,
        // post-seed cuts
        // PENDING
        // post-fit cuts
        // PENDING
        // --
        kNChannelACuts,
    };
    enum class EChannelD : int {
        kAllCombinations,
        // pre-seed cuts
        kPasses_DiffTracks_Logical,
        kPasses_DiffTracks_Physical,
        // post-seed cuts
        // PENDING
        // post-fit cuts
        // PENDING
        // --
        kNChannelDCuts,
    };
    enum class EChannelH : int {
        kAllCombinations,
        // pre-seed cuts
        kPasses_DiffTracks_Logical,
        kPasses_DiffTracks_Physical,
        // post-seed cuts
        // PENDING
        // post-fit cuts
        // PENDING
        // --
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

    // tracks //

    bool PassesCuts_Proton(const POD::Track &track, TH1D *hist_cut_flow) const;
    bool PassesCuts_Kaon(const POD::Track &track, TH1D *hist_cut_flow) const;
    bool PassesCuts_Pion(const POD::Track &track, TH1D *hist_cut_flow) const;

    POD::Extended::McParticle BuildMcTrack(unsigned int track_mc_entry, int pdg_code_hypothesis, bool include_gm);

    // V0s //
    void FindV0s(const DB::Particles::Definition &pid);

    template <typename E>
    bool PreSeedCuts(const POD::Track &track_neg, const POD::Track &track_pos, TH1D *hist_cut_flow) const;
    bool PostSeedCuts_Lambda(const Seeder::PCA &pca_neg, const Seeder::PCA &pca_pos, TH1D *hist_cut_flow) const;
    bool PostSeedCuts_KaonZeroShort(const Seeder::PCA &pca_neg, const Seeder::PCA &pca_pos, TH1D *hist_cut_flow) const;
    bool PostFitCuts_Lambda(const Cached::V0 &v0, TH1D *hist_cut_flow) const;
    bool PostFitCuts_KaonZeroShort(const Cached::V0 &v0, TH1D *hist_cut_flow) const;

    POD::Extended::McParticle BuildMcV0(const POD::Extended::McParticle &mc_neg, const POD::Extended::McParticle &mc_pos, int pdg_code_hypothesis);
    POD::V0 Create_V0(const KF::FitResult &fit, const Seeder::PCA &neg_pca_wrt_v0, const Seeder::PCA &pos_pca_wrt_v0);

    // mc sexaquark //
    POD::Linked::InjectedSexa BuildMcSexaquark(const POD::Extended::McParticle &mc_dau1, const POD::Extended::McParticle &mc_dau2);

    // channel A //
    void FindSexaquarks_ChannelA(bool is_bkg_channel);
    [[nodiscard]] bool PreSeedCuts_ChannelA(const POD::V0 &lambda, const POD::Track &lambda_neg, const POD::Track &lambda_pos, const POD::V0 &k0s,
                                            const POD::Track &k0s_neg, const POD::Track &k0s_pos, TH1D *hist_cut_flow) const;
    [[nodiscard]] bool PostSeedCuts_ChannelA(const Seeder::PCA &pca_v0a, const Seeder::PCA &pca_v0b, TH1D *hist_cut_flow) const;
    [[nodiscard]] bool PostFitCuts_ChannelA(const Cached::ChannelA &c_sexa, TH1D *hist_cut_flow) const;
    POD::Sexaquark Create_ChannelA(const KF::FitResult &fit, const Seeder::PCA &pca_v0a, const Seeder::PCA &pca_v0b, bool is_bkg_channel);

    // channel D //
    void FindSexaquarks_ChannelD(bool is_bkg_channel);
    [[nodiscard]] bool PreSeedCuts_ChannelD(const POD::Track &lambda_neg, const POD::Track &lambda_pos, const POD::Track &kaon,
                                            TH1D *hist_cut_flow) const;
    [[nodiscard]] bool PostSeedCuts_ChannelD(const Seeder::PCA &pca_v0, const Seeder::PCA &pca_ka, TH1D *hist_cut_flow) const;
    [[nodiscard]] bool PostFitCuts_ChannelD(const Cached::ChannelD &c_sexa, TH1D *hist_cut_flow) const;
    POD::Sexaquark Create_ChannelD(const KF::FitResult &fit, const Seeder::PCA &pca_v0, const Seeder::PCA &pca_ka, bool is_bkg_channel);

    // channel H //
    void FindSexaquarks_ChannelH(bool is_bkg_channel);
    [[nodiscard]] bool PreSeedCuts_ChannelH(const POD::Track &kaon1, const POD::Track &kaon2, TH1D *hist_cut_flow) const;
    [[nodiscard]] bool PostSeedCuts_ChannelH(const Seeder::PCA &pca_kaon1, const Seeder::PCA &pca_kaon2, TH1D *hist_cut_flow) const;
    [[nodiscard]] bool PostFitCuts_ChannelH(const Cached::ChannelH &c_sexa, TH1D *hist_cut_flow) const;
    POD::Sexaquark Create_ChannelH(const KF::FitResult &fit, const Seeder::PCA &pca_kaon1, const Seeder::PCA &pca_kaon2, bool is_bkg_channel);

    // fit configuration //

    static constexpr bool kMassConstraints = true;  // HARDCODED; if enabled, allow all mass constraints, except for the antisexaquark mass
    static constexpr KF::FitPolicy GetPolicy_V0s(double mass) {
        if (kMassConstraints) {
            return {
                .pin_daughters = true,
                .daughters_already_pinned = false,
                .mother_mass = mass,
                .prod_vertex = std::nullopt,
            };
        }
        return {};
    }
    static constexpr KF::FitPolicy GetPolicy_SV() {
        if (kMassConstraints) {
            return {
                .pin_daughters = true,
                .daughters_already_pinned = true,
                .mother_mass = std::nullopt,
                .prod_vertex = std::nullopt,
            };
        }
        return {};
    }

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
    // -- cut flow for antisexaquarks + bkg. sexaquarks
    std::unique_ptr<TH1D> fHist_CutFlow_ChannelA;
    std::unique_ptr<TH1D> fHist_CutFlow_ChannelA_Bkg;
    std::unique_ptr<TH1D> fHist_CutFlow_ChannelD;
    std::unique_ptr<TH1D> fHist_CutFlow_ChannelD_Bkg;
    std::unique_ptr<TH1D> fHist_CutFlow_ChannelH;
    std::unique_ptr<TH1D> fHist_CutFlow_ChannelH_Bkg;
};

}  // namespace T2DS
