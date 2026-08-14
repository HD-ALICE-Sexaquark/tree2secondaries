#pragma once

#include <memory>

#include <TFile.h>
#include <TH1.h>

#include <Math/Point3D.h>

#include "common/Constants.hpp"
#include "common/Framework.hpp"
#include "common/Framework_TeeTree.hpp"
#include "common/HD_Library.hpp"
#include "common/Schema_Events.hpp"
#include "common/Schema_FoundHdibaryon.hpp"

#include "App/Settings.hxx"
#include "KalmanFitter/KalmanFitterParticle.hxx"

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
        // fast cuts
        kPasses_Max_DCAbtwDaughters,
        // slow cuts
        kPasses_AbsMax_Pz,
        kPasses_Max_Pt,
        kPasses_Min_Pt,
        kPasses_AbsMax_Rapidity,
        kPasses_Min_Mass,
        kPasses_Max_Mass,
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
        // fast cuts
        kPasses_Max_DCAbtwDau,
        // slow cuts
        kPasses_AbsMax_Pz,
        kPasses_Max_Pt,
        kPasses_Min_Pt,
        kPasses_Min_Mass,
        kPasses_Max_Mass,
        kPasses_AbsMax_Rapidity,
        kPasses_Max_DecayLength,
        kPasses_Min_CPAwrtPV,
        kPasses_Max_Chi2NDF,
        // (anti)lambdas : depend on (anti)h-dibaryon decay vertex
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
          fInput_File{std::make_unique<TFile>(fSettings.PathInputFile.c_str(), "READ")},
          fInput{},
          fReader{nullptr},
          // output
          fOutput_File{std::make_unique<TFile>(fSettings.PathOutputFile.c_str(), "RECREATE")},
          fOutput{},
          fWriter{nullptr} {

        fReader = std::make_unique<Framework::TeeTree::Reader>(fInput.CreateModel_TeeTree(fSettings.IsMC, false), E2T::Name_OutputTree, *fInput_File);
        fWriter = std::make_unique<Framework::Writer>(fOutput.CreateModel(fSettings.IsMC), T2DS::Name_FoundHdibaryonRNT, *fOutput_File);

        PrepareOutputHistograms();

        Logger::Info(__FUNCTION__, "Verifier initialized successfully.");
    }

    void PrepareOutputHistograms();

    void Load(long long entry_idx) { fReader->Load(entry_idx); }
    [[nodiscard]] long long NumberEventsToRead() {
        auto total = fReader->GetEntries();
        return fSettings.LimitToNEvents.has_value() ? std::min(static_cast<long long>(fSettings.LimitToNEvents.value()), total) : total;
    }

    void ProcessEvent();
    void ProcessInjected();
    void ProcessPreFoundLambda();

    void Verify() {
        VerifyLambdaPair(true);
        VerifyLambdaPair(false);
    }

    void EndOfEvent();
    void EndOfAnalysis();

   private:
    // injected //
    POD::InjectedHdib BuildInjectedHdibaryon(const POD::McParticle &mc);

    // mc charged track //
    POD::Extended::McParticle BuildMcTrack(unsigned int track_mc_entry, const HD::DecayTree &decay_pid, int pdg_code_hypothesis);
    [[nodiscard]] POD::Track ExtractTrack(const POD::PreFoundLambda &pod_lambda, short charge) const;

    // pre-found on-the-fly (anti)lambdas //
    [[nodiscard]] bool FastCuts_Lambda(const Seeder::PCA &pca_neg, const Seeder::PCA &pca_pos);
    bool SlowCuts_Lambda(const Cached::PreFoundLambda &c_lambda, bool anti_channel);

    POD::Extended::McParticle BuildMcPreFoundLambda(const POD::Extended::McParticle &mc_neg, const POD::Extended::McParticle &mc_pos,
                                                    const HD::DecayTree &decay_pid);
    POD::Extended::PreFoundLambda CreateExtendedPreFoundLambda(const POD::PreFoundLambda &old_lambda, const KF::FitResult &fit,
                                                               const Seeder::PCA &pca_neg, const Seeder::PCA &pca_pos);

    // h-dibaryon //
    void VerifyLambdaPair(bool anti_channel);

    [[nodiscard]] bool FastCuts_Hdibaryon(const Seeder::PCA &pca_lambda1, const Seeder::PCA &pca_lambda2, TH1D *hist_cut_flow);
    [[nodiscard]] bool SlowCuts_Hdibaryon(const Cached::Hdibaryon &c_hdib, TH1D *hist_cut_flow);

    POD::Extended::McParticle BuildMcHdibaryon(const POD::Extended::McParticle &mc_lambda1, const POD::Extended::McParticle &mc_lambda2,
                                               const HD::DecayTree &decay_pid);
    POD::LambdaPair CreateLambdaPair(const KF::FitResult &fit, const Seeder::PCA &pca_lambda1, const Seeder::PCA &pca_lambda2, bool anti_channel);

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
    std::vector<POD::Extended::McParticle> fTemp_MC_AntiLambda;
    std::vector<POD::Extended::McParticle> fTemp_MC_AntiLambda_Neg;
    std::vector<POD::Extended::McParticle> fTemp_MC_AntiLambda_Pos;
    std::vector<POD::Extended::McParticle> fTemp_MC_Lambda;
    std::vector<POD::Extended::McParticle> fTemp_MC_Lambda_Neg;
    std::vector<POD::Extended::McParticle> fTemp_MC_Lambda_Pos;

    // input //

    std::unique_ptr<TFile> fInput_File;
    Schema::Events fInput;
    std::unique_ptr<Framework::TeeTree::Reader> fReader;
    // -- cached
    ROOT::Math::XYZPoint fPrimaryVertex;
    KF::Vertex fPrimaryVertexKF;
    double fMagneticField{0.};

    // output //

    std::unique_ptr<TFile> fOutput_File;
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
};

}  // namespace T2DS
