#pragma once

#include <memory>

#include <TFile.h>
#include <TH1.h>

#include <Math/Point3D.h>
#include <ROOT/RNTupleReader.hxx>

#include "common/Framework.hpp"
#include "common/HD_Library.hpp"
#include "common/Schema_Events.hpp"
#include "common/Schema_FoundHdib.hpp"

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
          fReader{E2R::Name_OutputRNT, *fInput_File},
          fInput{fReader.Data()},
          // output
          fOutput_File{std::make_unique<TFile>(fSettings.PathOutputFile.c_str(), "RECREATE")},
          fWriter{R2DS::Name_FoundHdibRNT, *fOutput_File},
          fOutput{fWriter.Data()} {

        PrepareOutputHistograms();

        Logger::Info(__FUNCTION__, "Verifier initialized successfully.");
    }

    void PrepareOutputHistograms();

    void Load(ROOT::NTupleSize_t entry_id) { fReader.Load(entry_id); }
    [[nodiscard]] ROOT::NTupleSize_t NumberEventsToRead() {
        auto total = fReader.Iter()->GetNEntries();
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
    // on-the-fly lambda //

    void BuildMcInfo(POD::OnTheFlyLambda &new_lambda, const HD::DecayTree &decay_pid);

    // h-dibaryon //

    void VerifyLambdaPair(bool anti_channel);

    [[nodiscard]] bool Cuts(const KF::LambdaPair &hdibaryon, TH1D *cut_flow_hist) const;

    void BuildMcInfo(POD::LambdaPair &new_v0, const HD::DecayTree &decay_pid);
    void BuildRecInfo(POD::LambdaPair &new_v0, const KF::LambdaPair &kf_v0, bool anti_channel);

    // member variables //

    const Settings &fSettings;

    std::unique_ptr<TH1D> fHist_EventCounter;

    std::unique_ptr<TH1D> fHist_CutFlow;
    std::unique_ptr<TH1D> fHist_CutFlow_AntiChannel;

    // input //

    std::unique_ptr<TFile> fInput_File;
    Framework::Reader<Schema::Events> fReader;
    Schema::Events &fInput;
    // -- cached
    ROOT::Math::XYZPoint fPrimaryVertex;

    // output //

    std::unique_ptr<TFile> fOutput_File;
    Framework::Writer<Schema::FoundHdib> fWriter;
    Schema::FoundHdib &fOutput;
};

}  // namespace R2DS
