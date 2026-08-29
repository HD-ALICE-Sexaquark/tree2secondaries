#pragma once

#include <exception>
#include <memory>
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
#include "KalmanFitter/KalmanFitterFitTypes.hxx"

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

    // V0s //

    void FindV0s(const DB::Particles::Definition &pid);

    template <typename E>
    bool PreSeedCuts(const POD::Track &track_neg, const POD::Track &track_pos, TH1D *hist_cut_flow) const;
    bool PostSeedCuts_Lambda(const Seeder::PCA &pca_neg, const Seeder::PCA &pca_pos, TH1D *hist_cut_flow) const;
    bool PostSeedCuts_KaonZeroShort(const Seeder::PCA &pca_neg, const Seeder::PCA &pca_pos, TH1D *hist_cut_flow) const;
    bool PostFitCuts_Lambda(const Cached::V0 &v0, TH1D *hist_cut_flow) const;
    bool PostFitCuts_KaonZeroShort(const Cached::V0 &v0, TH1D *hist_cut_flow) const;

    POD::V0 Create_V0(const KF::FitResult &fit, const Seeder::PCA &neg_pca_wrt_v0, const Seeder::PCA &pos_pca_wrt_v0);

    // mc sexaquark //
    POD::Linked::InjectedSexa Link_InjectedSexa(const POD::Extended::McParticle &mc_dau1, const POD::Extended::McParticle &mc_dau2);

    // channel A //
    void FindSexaquarks_ChannelA(bool wrong_sign_channel);
    [[nodiscard]] bool PreSeedCuts_ChannelA(const POD::V0 &lambda, const POD::Track &lambda_neg, const POD::Track &lambda_pos, const POD::V0 &k0s,
                                            const POD::Track &k0s_neg, const POD::Track &k0s_pos, TH1D *hist_cut_flow) const;
    [[nodiscard]] bool PostSeedCuts_ChannelA(const Seeder::PCA &pca_v0a, const Seeder::PCA &pca_v0b, TH1D *hist_cut_flow) const;
    [[nodiscard]] bool PostFitCuts_ChannelA(const Cached::ChannelA &c_sexa, TH1D *hist_cut_flow) const;
    POD::Sexaquark Create_ChannelA(const KF::FitResult &fit, const Seeder::PCA &pca_v0a, const Seeder::PCA &pca_v0b, bool wrong_sign_channel);

    // channel D //
    void FindSexaquarks_ChannelD(bool wrong_sign_channel);
    [[nodiscard]] bool PreSeedCuts_ChannelD(const POD::Track &lambda_neg, const POD::Track &lambda_pos, const POD::Track &kaon,
                                            TH1D *hist_cut_flow) const;
    [[nodiscard]] bool PostSeedCuts_ChannelD(const Seeder::PCA &pca_v0, const Seeder::PCA &pca_ka, TH1D *hist_cut_flow) const;
    [[nodiscard]] bool PostFitCuts_ChannelD(const Cached::ChannelD &c_sexa, TH1D *hist_cut_flow) const;
    POD::Sexaquark Create_ChannelD(const KF::FitResult &fit, const Seeder::PCA &pca_v0, const Seeder::PCA &pca_ka, bool wrong_sign_channel);

    // channel H //
    void FindSexaquarks_ChannelH(bool wrong_sign_channel);
    [[nodiscard]] bool PreSeedCuts_ChannelH(const POD::Track &kaon1, const POD::Track &kaon2, TH1D *hist_cut_flow) const;
    [[nodiscard]] bool PostSeedCuts_ChannelH(const Seeder::PCA &pca_kaon1, const Seeder::PCA &pca_kaon2, TH1D *hist_cut_flow) const;
    [[nodiscard]] bool PostFitCuts_ChannelH(const Cached::ChannelH &c_sexa, TH1D *hist_cut_flow) const;
    POD::Sexaquark Create_ChannelH(const KF::FitResult &fit, const Seeder::PCA &pca_kaon1, const Seeder::PCA &pca_kaon2, bool wrong_sign_channel);

    // fit configuration //

    static constexpr bool kMassConstraints = true;  // HARDCODED; if enabled, allow all mass constraints, except for the antisexaquark mass
    static constexpr KF::FitPolicy GetPolicy_V0s(double mass) {
        if (kMassConstraints) {
            return {
                .pin_daughters = true,
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
                .mother_mass = std::nullopt,
                .prod_vertex = std::nullopt,
            };
        }
        return {};
    }
    static constexpr bool kV0sArePinned = kMassConstraints;  // re-assert that V0s were pinned to their nominal mass

    // member variables //

    const Settings &fSettings;

    // temporary in-memory data //

    std::vector<POD::Track> fTemp_AntiProton;
    std::vector<POD::Track> fTemp_Proton;
    std::vector<POD::Track> fTemp_NegKaon;
    std::vector<POD::Track> fTemp_PosKaon;
    std::vector<POD::Track> fTemp_PiMinus;
    std::vector<POD::Track> fTemp_PiPlus;

    std::vector<unsigned int> fTemp_McEntry_AntiProton;
    std::vector<unsigned int> fTemp_McEntry_Proton;
    std::vector<unsigned int> fTemp_McEntry_NegKaon;
    std::vector<unsigned int> fTemp_McEntry_PosKaon;
    std::vector<unsigned int> fTemp_McEntry_PiMinus;
    std::vector<unsigned int> fTemp_McEntry_PiPlus;

    std::vector<POD::V0> fTemp_AntiLambda;
    std::vector<POD::Track> fTemp_AntiLambda_Neg;
    std::vector<POD::Track> fTemp_AntiLambda_Pos;
    std::vector<POD::V0> fTemp_Lambda;
    std::vector<POD::Track> fTemp_Lambda_Neg;
    std::vector<POD::Track> fTemp_Lambda_Pos;
    std::vector<POD::V0> fTemp_KaonZeroShort;
    std::vector<POD::Track> fTemp_KaonZeroShort_Neg;
    std::vector<POD::Track> fTemp_KaonZeroShort_Pos;

    std::vector<unsigned int> fTemp_McEntry_AntiLambda_Neg;
    std::vector<unsigned int> fTemp_McEntry_AntiLambda_Pos;
    std::vector<unsigned int> fTemp_McEntry_Lambda_Neg;
    std::vector<unsigned int> fTemp_McEntry_Lambda_Pos;
    std::vector<unsigned int> fTemp_McEntry_KaonZeroShort_Neg;
    std::vector<unsigned int> fTemp_McEntry_KaonZeroShort_Pos;

    // input //

    Schema::Events fInput;
    std::unique_ptr<Framework::TeeTree::Reader> fReader;
    // -- cached
    ROOT::Math::XYZPoint fPrimaryVertex;
    double fMagneticField{0.};
    DB::ReactionChannels::Definition fMcSignalChannel{DB::ReactionChannels::ReactionChannel('0')};  // identified from the first event of the mc stack
    bool fMcSignalChannelChecked{false};

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
