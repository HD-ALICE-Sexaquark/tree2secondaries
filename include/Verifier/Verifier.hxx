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

// forward declarations //
// clang-format off
namespace POD { struct InjectedHdib; struct PreFoundLambda; struct LambdaPair; }
namespace Cached { struct PreFoundLambda; struct Hdibaryon; }

namespace T2DS {

namespace Seeder { struct PCA; }
namespace KF { struct Particle; }
// clang-format on

class Verifier {

    enum EPreFoundLambda {
        kAvailablePreFoundLambdas,
        kPassesDaughtersPID,
        kPassesRapidityCut,
        kPassesMinPtProton,
        kNPreFoundLambdaCuts,
    };

    enum ELambdaPair {
        kAllCombinations,
        kNLambdaPairCuts,
    };

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
    POD::Extended::PreFoundLambda CreateExtendedPreFoundLambda(const POD::PreFoundLambda &old_lambda, const KF::Particle &fit,
                                                               const Seeder::PCA &pca_neg, const Seeder::PCA &pca_pos, double mass_neg,
                                                               double mass_pos);

    // h-dibaryon //
    void VerifyLambdaPair(bool anti_channel);

    [[nodiscard]] bool FastCuts_Hdibaryon(const Seeder::PCA &pca_lambda1, const Seeder::PCA &pca_lambda2, TH1D *hist_cut_flow);
    [[nodiscard]] bool SlowCuts_Hdibaryon(const Cached::Hdibaryon &c_hdib, TH1D *hist_cut_flow);

    POD::Extended::McParticle BuildMcHdibaryon(const POD::Extended::McParticle &mc_lambda1, const POD::Extended::McParticle &mc_lambda2,
                                               const HD::DecayTree &decay_pid);
    POD::LambdaPair CreateLambdaPair(const KF::Particle &fit, const Seeder::PCA &pca_lambda1, const Seeder::PCA &pca_lambda2, bool anti_channel);

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
