#pragma once

#include <memory>

#include <TFile.h>
#include <TH1.h>

#include <Math/Point3D.h>
#include <ROOT/RNTupleReader.hxx>

#include "common/Constants.hpp"
#include "common/Framework.hpp"
#include "common/HD_Library.hpp"
#include "common/Schema_Events.hpp"
#include "common/Schema_FoundHdibaryon.hpp"

#include "App/Settings.hxx"

// forward declarations //
// clang-format off
namespace POD {
    struct OnTheFlyLambda;
    struct LambdaPair;
}
namespace KF { struct LambdaPair; }
// clang-format on

namespace R2DS {

class Verifier {
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

        fReader = std::make_unique<Framework::Reader>(fInput.CreateModel(fSettings.IsMC, fSettings.IsMC), E2R::Name_OutputRNT, *fInput_File);
        fWriter = std::make_unique<Framework::Writer>(fOutput.CreateModel(fSettings.IsMC), R2DS::Name_FoundHdibaryonRNT, *fOutput_File);

        PrepareOutputHistograms();

        Logger::Info(__FUNCTION__, "Verifier initialized successfully.");
    }

    void PrepareOutputHistograms();

    void Load(ROOT::NTupleSize_t entry_id) { fReader->Load(entry_id); }
    [[nodiscard]] ROOT::NTupleSize_t NumberEventsToRead() {
        auto total = fReader->Iter()->GetNEntries();
        return fSettings.LimitToNEvents.has_value() ? std::min(fSettings.LimitToNEvents.value(), total) : total;
    }

    void ProcessEvent();
    void ProcessInjected();

    void Verify() {
        VerifyLambdaPair(false);
        VerifyLambdaPair(true);
    }

    void EndOfEvent();
    void EndOfAnalysis();

   private:
    // injected //
    POD::Extended::McParticle BuildInjectedHdibaryon(const POD::McParticle &mc);

    // mc charged track //
    POD::Extended::McParticle BuildMcTrack(unsigned int track_mc_entry, const HD::DecayTree &decay_pid, int pdg_code_hypothesis);

    // on-the-fly lambda //
    POD::Extended::McParticle BuildMcOnTheFlyLambda(const POD::Extended::McParticle &mc_neg, const POD::Extended::McParticle &mc_pos,
                                                    const HD::DecayTree &decay_pid);

    // h-dibaryon //
    POD::Extended::McParticle BuildMcHdibaryon(const POD::Extended::McParticle &mc_lambda1, const POD::Extended::McParticle &mc_lambda2,
                                               const HD::DecayTree &decay_pid);
    POD::LambdaPair CreateLambdaPair(const KF::LambdaPair &kf_hdib, bool anti_channel);
    void VerifyLambdaPair(bool anti_channel);
    [[nodiscard]] bool Cuts(const KF::LambdaPair &hdibaryon, TH1D *cut_flow_hist) const;

    // member variables //

    const Settings &fSettings;

    // input //

    std::unique_ptr<TFile> fInput_File;
    Schema::Events fInput;
    std::unique_ptr<Framework::Reader> fReader;
    // -- cached
    ROOT::Math::XYZPoint fPrimaryVertex;

    // output //

    std::unique_ptr<TFile> fOutput_File;
    Schema::FoundHdibaryon fOutput;
    std::unique_ptr<Framework::Writer> fWriter;

    std::unique_ptr<TH1D> fHist_EventCounter;

    std::unique_ptr<TH1D> fHist_CutFlow;
    std::unique_ptr<TH1D> fHist_CutFlow_AntiChannel;
};

}  // namespace R2DS
