#pragma once

#include <memory>

#include <TFile.h>
#include <TH1.h>

#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleWriter.hxx>

#include "common/FL_Event.hpp"
#include "common/FL_InjectedHdib.hpp"
#include "common/FL_LambdaPair.hpp"
#include "common/VC_McParticle.hpp"
#include "common/VC_McParticleView.hpp"
#include "common/VC_OnTheFlyLambda.hpp"
#include "common/VC_OnTheFlyLambdaView.hpp"

#include "App/Settings.hxx"
#include "KalmanFitter/KF_LambdaPair.hxx"

namespace R2DS {

namespace MC {
struct LambdaPair {
    const Vector::McParticleView *lambda1;
    const Vector::McParticleView *lambda1_neg;
    const Vector::McParticleView *lambda1_pos;
    const Vector::McParticleView *lambda2;
    const Vector::McParticleView *lambda2_neg;
    const Vector::McParticleView *lambda2_pos;
};
}  // namespace MC

class Verifier {
   public:
    Verifier(const Verifier &) = delete;
    Verifier(Verifier &&) = delete;
    Verifier &operator=(const Verifier &) = delete;
    Verifier &operator=(Verifier &&) = delete;
    ~Verifier() = default;

    explicit Verifier(const Settings &settings) : fSettings{settings} {}

    bool Initialize();

    bool PrepareOutputFile();
    void PrepareOutputHistograms();

    void ProcessEvent();

    void ProcessInjected();

    void Verify() {
        VerifyLambdaPair(false);
        VerifyLambdaPair(true);
    }

    [[nodiscard]] unsigned long NumberEventsToRead() const {
        unsigned long total_entries = fReader->GetNEntries();
        return 0 < fSettings.LimitToNEvents && fSettings.LimitToNEvents < total_entries ? fSettings.LimitToNEvents : total_entries;
    }

    void EndOfAnalysis();

    std::unique_ptr<ROOT::RNTupleModel> fInput_Model;
    std::unique_ptr<ROOT::RNTupleReader> fReader;

    std::unique_ptr<ROOT::RNTupleModel> fOutput_Model;
    std::unique_ptr<ROOT::RNTupleWriter> fWriter;

    std::unique_ptr<ROOT::RNTupleModel> fOutput_Model_InjectedHdib;
    std::unique_ptr<ROOT::RNTupleWriter> fWriter_InjectedHdib;

   private:
    void Assign_Event();
    void Assign_Candidate(const KalmanFitter::LambdaPair &hdibaryon, bool anti_channel);
    void Assign_InjectedHdib(const Vector::McParticleView *injected, Flat::InjectedHdib &out, bool embedded_to_rec);
    void Assign(const Vector::OnTheFlyLambdaView &lambda, const Vector::McParticleView *mc_lambda, const Vector::McParticleView *mc_neg,
                const Vector::McParticleView *mc_pos, Flat::OnTheFlyLambda &out);

    void VerifyLambdaPair(bool anti_channel);
    [[nodiscard]] bool Cuts(const KalmanFitter::LambdaPair &hdibaryon, TH1D *cut_flow_hist) const;
    std::unique_ptr<Vector::McParticleView> LinkInjectedSignal(const MC::LambdaPair &hdibaryon);
    void Assign(const KalmanFitter::LambdaPair &hdibaryon, const MC::LambdaPair *mc_hdibaryon, bool anti_channel);

    const Settings &fSettings;

    std::unique_ptr<TFile> fOutput_File;

    std::unique_ptr<TH1D> fHist_EventCounter;

    std::unique_ptr<TH1D> fHist_CutFlow;
    std::unique_ptr<TH1D> fHist_CutFlow_AntiChannel;

    // input //

    Flat::Event fInput_Event;
    ROOT::Math::XYZPoint fPrimaryVertex;

    Vector::McParticle fInput_McParticle;
    Vector::OnTheFlyLambda fInput_OnTheFlyLambda;

    // output //

    Flat::InjectedHdib fOutput_InjectedHdib;

    Flat::LambdaPair fOutput_LambdaPair;
};

}  // namespace R2DS
