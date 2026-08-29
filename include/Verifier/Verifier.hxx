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
#include "common/Framework.hpp"
#include "common/Framework_TeeTree.hpp"
#include "common/Schema_Events.hpp"
#include "common/Schema_FoundHdibaryon.hpp"

#include "App/Logger.hxx"
#include "App/Settings.hxx"
#include "KalmanFitter/KalmanFitterVertex.hxx"

// forward declarations //
// clang-format off
namespace POD { struct InjectedHdib; struct PreFoundLambda; struct LambdaPair; }
namespace Cached { struct PreFoundLambda; struct Hdibaryon; }

namespace T2DS {

namespace Seeder { struct PCA; }
namespace KF { struct FitResult; }
// clang-format on

class Verifier {

    // Fit-Related //

    enum class EMassConstraints : int {
        kNone,                                           // baseline, don't use any kind of mass constraint during/after fits
        kPinDaughtersOnFirstFit_DontPinMothers,          // protons+pions on shell; mass cuts available
        kPinDaughtersOnBothFits_PinMothersOnlyFirstFit,  // costs +1 NDF and a mass-pull chi2 term; lambdas' masses become deltas
    };
    struct FitSetup {
        bool pin_lambda_daughters{false};  // 1st fit: (anti)protons and pions onto their own shells
        bool pin_lambda_mass{false};       // 1st fit: the (anti)lambda onto m(lambda); costs +1 NDF and a chi2 term
        bool lambdas_on_shell{false};      // 2nd fit: the (anti)lambdas arrive carrying an exact mass hypothesis ...
        bool pin_lambdas{false};           // 2nd fit: ... and get re-pinned to it, after the update, before the sum
        bool pin_hdib_to_pv{false};        // 2nd fit, final step: pin particle to known production vertex
    };

    // Cut-related //

    enum class EPreFoundLambda : int {
        kAllPreFoundLambdas,
        // pre-seed cuts
        kPasses_DiffDaughters_Logical,
        kPasses_DiffDaughters_Physical,
        // post-seed cuts
        kPasses_Max_DCAbtwDaughters,
        // post-fit cuts
        kPasses_AbsMax_Pz,
        kPasses_Max_Pt,
        kPasses_Min_Pt,
        kPasses_AbsMax_Rapidity,
        kPasses_Min_CPAwrtPV,
        kPasses_AbsMax_ArmRadiusDev,
        kPasses_Max_DCAwrtPV,
        kPasses_Max_Chi2NDF,
        // -- depend on (anti)protons
        kPasses_AbsMax_Pz_Proton,
        kPasses_Max_Pt_Proton,
        kPasses_Min_Pt_Proton,
        // -- depend on pi(minus/plus)
        kPasses_AbsMax_Pz_Pion,
        kPasses_Max_Pt_Pion,
        kPasses_Min_Pt_Pion,
        // --
        kNPreFoundLambdaCuts,
    };
    enum class ELambdaPair : int {
        kAllCombinations,
        // pre-seed cuts
        kPasses_DiffLambdas_Logical,
        kPasses_DiffTracks_Logical,
        kPasses_DiffTracks_Physical,
        kPasses_DiffLambdas_Physical,
        // post-seed cuts
        kPasses_Max_DCAbtwDau,
        // post-fit cuts
        kPasses_AbsMax_Pz,
        kPasses_Max_Pt,
        kPasses_Min_Pt,
        kPasses_Min_Mass,
        kPasses_Max_Mass,
        kPasses_AbsMax_Rapidity,
        kPasses_Max_DecayLength,
        kPasses_Min_CPAwrtPV,
        kPasses_Max_Chi2NDF,
        kPasses_Max_Chi2CV,
        // further cuts on (anti)lambdas, as they depend on (anti)h-dibaryon decay vertex
        kPasses_Max_L1_DecayLength,
        kPasses_Min_L1_DecayLength,
        kPasses_Min_L1_CPAwrtDV,
        kPasses_Max_L2_DecayLength,
        kPasses_Min_L2_DecayLength,
        kPasses_Min_L2_CPAwrtDV,
        // --
        kNLambdaPairCuts,
    };
    template <typename E>
    void FillHist(TH1D *hist, const E &bin_n) {
        hist->Fill(static_cast<double>(bin_n));
    }

   public:
    Verifier(const Verifier &) = delete;
    Verifier(Verifier &&) = delete;
    Verifier &operator=(const Verifier &) = delete;
    Verifier &operator=(Verifier &&) = delete;
    ~Verifier() = default;

    explicit Verifier(const Settings &settings)
        : fSettings{settings},
          // input
          fInput{},
          fReader{nullptr},
          // output
          fOutput_File{std::make_unique<TFile>(fSettings.PathOutputFile.c_str(), "RECREATE")},
          fOutput{},
          fWriter{nullptr} {

        fWriter = std::make_unique<Framework::Writer>(fOutput.CreateModel(fSettings.IsMC), T2DS::Name_FoundHdibaryonRNT, *fOutput_File);

        PrepareOutputHistograms();

        Logger::Info(__FUNCTION__, "Verifier initialized successfully.");
    }

    [[nodiscard]] bool OpenInput(std::string_view path) {
        fReader.reset();  // raw `TTree*` must die first
        try {
            fReader = std::make_unique<Framework::TeeTree::Reader>(fInput.CreateModel_TeeTree(fSettings.IsMC, false), E2T::Name_OutputTree, path);
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
    void ProcessInjected();
    void ProcessPreFoundLambda();

    void Verify() {
        VerifyLambdaPair(true, true);    // anti-lambda + anti-lambda
        VerifyLambdaPair(false, false);  // lambda + lambda
        VerifyLambdaPair(false, true);   // mixed: lambda + anti-lambda
    }

    void EndOfEvent();
    bool EndOfAnalysis();

   private:
    // injected //
    POD::InjectedHdib Create_InjectedHdibaryon(const POD::McParticle &mc);

    // mc charged track //
    [[nodiscard]] POD::Track ExtractTrack(const POD::PreFoundLambda &pod_lambda, short charge) const;

    // pre-found on-the-fly (anti)lambdas //
    [[nodiscard]] bool PreSeedCuts_Lambda(const POD::PreFoundLambda &lambda);
    [[nodiscard]] bool PostSeedCuts_Lambda(const Seeder::PCA &pca_neg, const Seeder::PCA &pca_pos);
    bool PostFitCuts_Lambda(const Cached::PreFoundLambda &c_lambda);

    POD::Extended::PreFoundLambda Extend_PreFoundLambda(const POD::PreFoundLambda &old_lambda, const KF::FitResult &fit, const Seeder::PCA &pca_neg,
                                                        const Seeder::PCA &pca_pos, bool anti_lambda);

    // h-dibaryon //
    void VerifyLambdaPair(bool anti_channel_l1, bool anti_channel_l2);

    [[nodiscard]] bool PreSeedCuts_Hdibaryon(const POD::Extended::PreFoundLambda &lambda1, const POD::Extended::PreFoundLambda &lambda2,
                                             TH1D *hist_cut_flow);
    [[nodiscard]] bool PostSeedCuts_Hdibaryon(const Seeder::PCA &pca_lambda1, const Seeder::PCA &pca_lambda2, TH1D *hist_cut_flow);
    [[nodiscard]] bool PostFitCuts_Hdibaryon(const Cached::Hdibaryon &c_hdib, TH1D *hist_cut_flow);

    POD::LambdaPair Create_LambdaPair(const KF::FitResult &fit, const Seeder::PCA &pca_lambda1, const Seeder::PCA &pca_lambda2);

    // fit configuration //
    // -- (anti)h-dibaryon mass is never pinned in any configuration, as it's the property under study
    [[nodiscard]] static constexpr FitSetup GetFitSetup() {
        switch (kMassConstraints) {
            // -- nothing is constrained anywhere
            case EMassConstraints::kNone:
                return {.pin_hdib_to_pv = kProdVertexConstraint};
            // -- every daughter is set on shell, no mother is ever pinned.
            //    The (anti)lambda mass stays a measurement, so its cuts still select.
            case EMassConstraints::kPinDaughtersOnFirstFit_DontPinMothers:
                return {.pin_lambda_daughters = true,
                        .pin_lambda_mass = false,
                        .lambdas_on_shell = false,
                        .pin_lambdas = false,
                        .pin_hdib_to_pv = kProdVertexConstraint};
            // -- the same, plus the (anti)lambda is pinned to m(Lambda) in the first fit, which does pay +1 NDF
            //    and a mass-pull chi2 term, and turns its mass into a delta
            case EMassConstraints::kPinDaughtersOnBothFits_PinMothersOnlyFirstFit:
                return {.pin_lambda_daughters = true,
                        .pin_lambda_mass = true,
                        .lambdas_on_shell = true,
                        .pin_lambdas = true,
                        .pin_hdib_to_pv = kProdVertexConstraint};
        }
        return {.pin_hdib_to_pv = kProdVertexConstraint};  // unreachable, switch above is exhaustive
    }
    static constexpr EMassConstraints kMassConstraints = EMassConstraints::kPinDaughtersOnBothFits_PinMothersOnlyFirstFit;  // HARDCODED
    static constexpr bool kProdVertexConstraint = true;                                                                     // HARDCODED

    // member variables //

    const Settings &fSettings;

    // temporary pre-found lambdas, extended for KF usage //

    std::vector<POD::Extended::PreFoundLambda> fTemp_AntiLambda;
    std::vector<POD::Extended::PreFoundLambda> fTemp_Lambda;

    std::vector<unsigned int> fTemp_McEntry_AntiLambda_Neg;
    std::vector<unsigned int> fTemp_McEntry_AntiLambda_Pos;
    std::vector<unsigned int> fTemp_McEntry_Lambda_Neg;
    std::vector<unsigned int> fTemp_McEntry_Lambda_Pos;

    // input //

    Schema::Events fInput;
    std::unique_ptr<Framework::TeeTree::Reader> fReader;
    // -- cached
    ROOT::Math::XYZPoint fPrimaryVertex;
    KF::Vertex fPrimaryVertexKF;
    double fMagneticField{0.};

    // output //

    std::unique_ptr<TFile> fOutput_File;  // single file, kept alive across every input file, if multiple
    Schema::FoundHdibaryon fOutput;
    std::unique_ptr<Framework::Writer> fWriter;

    // histograms
    // -- event counter
    std::unique_ptr<TH1D> fHist_EventCounter;
    // -- cut flow for (anti)lambdas
    std::unique_ptr<TH1D> fHist_CutFlow_AntiLambda;
    std::unique_ptr<TH1D> fHist_CutFlow_Lambda;
    // -- cut flow for (anti)h-dibaryons
    std::unique_ptr<TH1D> fHist_CutFlow_AntiHdibaryon;
    std::unique_ptr<TH1D> fHist_CutFlow_Hdibaryon;
    std::unique_ptr<TH1D> fHist_CutFlow_MixedLambdaPair;
};

}  // namespace T2DS
