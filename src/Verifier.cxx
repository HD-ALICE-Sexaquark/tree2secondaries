#include <cstddef>
#include <memory>
#include <tuple>

#include "common/Constants.hpp"
#include "common/DB_Particles.hpp"
#include "common/HD_Library.hpp"
#include "common/MC_Helpers.hpp"
#include "common/Math.hpp"
#include "common/POD_LambdaPair.hpp"
#include "common/POD_McParticle.hpp"
#include "common/POD_OnTheFlyLambda.hpp"
#include "common/R2DS_Cuts.hpp"

#include "KalmanFitter/KF_LambdaPair.hxx"
#include "Seeder/SeederLineLine.hxx"

#include "Verifier/Verifier.hxx"

namespace CMath = Common::Math;

namespace R2DS {

// ## OUTPUT ZONE ## //

void Verifier::PrepareOutputHistograms() {
    // event counter //
    fHist_EventCounter = std::make_unique<TH1D>("N_Events", ";;N_Events", 1, 0, 1);
    // cut flows //
    constexpr int x_nbins = 20;
    constexpr float x_min = 0.;
    constexpr float x_max = 20.;
    constexpr const char* hist_title = ";Cut N;N Passed Cut";
    fHist_CutFlow = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", DB::Particles::Particle("AntiProton").acronym).c_str(), hist_title, x_nbins, x_min, x_max);
    fHist_CutFlow_AntiChannel = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", DB::Particles::Particle("AntiLambda").acronym).c_str(), hist_title, x_nbins, x_min, x_max);
}

// ## Event ZONE ## //

void Verifier::ProcessEvent() {
    // -- copy event info
    fOutput.Event = fInput.Event;
    // -- update event counter
    fHist_EventCounter->Fill(0.);
    // -- cache pv
    fPrimaryVertex.SetCoordinates(fOutput.Event.PV_X, fOutput.Event.PV_Y, fOutput.Event.PV_Z);
}

// ## Injected ZONE ## //

void Verifier::ProcessInjected() {
    // alias collection //
    const auto& mc_collection = fInput.McParticle;
    // loop over mc particles //
    for (const auto& mc : mc_collection) {
        // select only injected h-dibaryons //
        if (std::abs(mc.PdgCode) != DB::Particles::Particle("Hdibaryon").pdg_code) continue;
        // store //
        fOutput.Injected.emplace_back(BuildInjectedHdibaryon(mc));
    }
}

POD::Extended::McParticle Verifier::BuildInjectedHdibaryon(const POD::McParticle& mc) {
    POD::Extended::McParticle inj(mc);
    std::tie(inj.Decay_X, inj.Decay_Y, inj.Decay_Z) = MC::GetDecayVertex(mc, fInput.McParticle);
    inj.SignalID = static_cast<int>(mc.StatusCode);
    inj.Mother_PdgCode = Common::DummyInt;
    inj.GM_McEntry = Common::DummyInt;
    inj.GM_PdgCode = Common::DummyInt;
    inj.IsTrue = true;
    inj.IsGen1Signal = false;
    inj.IsGen2Signal = false;
    inj.IsTrueSignal = true;
    inj.IsSecondary = false;
    inj.IsHybrid = false;

    return inj;
}

// ## Charged Tracks ZONE ## //

POD::Extended::McParticle Verifier::BuildMcTrack(unsigned int track_mc_entry, const HD::DecayTree& decay_pid, int pdg_code_hypothesis) {
    // copy linked mc info //
    POD::Extended::McParticle new_mc(fInput.McParticle[track_mc_entry]);
    MC::Apply(new_mc, MC::HdibaryonRules::ClassifyDownstream(new_mc, fInput.McParticle, decay_pid, pdg_code_hypothesis, false));
    return new_mc;
}

// ## Single On-The-Fly Lambda ZONE ## //

POD::Extended::McParticle Verifier::BuildMcOnTheFlyLambda(const POD::Extended::McParticle& mc_neg, const POD::Extended::McParticle& mc_pos,
                                                          const HD::DecayTree& decay_pid) {
    POD::Extended::McParticle mc_lambda;
    // clang-format off
    auto mc_entry = MC::FindCommonMotherMcEntry(mc_neg, mc_pos);
    if (!mc_entry.has_value()) { return mc_lambda; }
    // clang-format on
    // fill values //
    static_cast<POD::McParticle&>(mc_lambda) = fInput.McParticle[mc_entry.value()];
    MC::Apply(mc_lambda, MC::HdibaryonRules::ClassifyDownstream(mc_lambda, fInput.McParticle, decay_pid, decay_pid.lambda.pdg_code, true));
    mc_lambda.IsHybrid = mc_neg.IsTrueSignal != mc_pos.IsTrueSignal;

    return mc_lambda;
}

// ## H-dibaryon ZONE ## //

void Verifier::VerifyLambdaPair(bool anti_channel) {

    // determine rules based on reconstruction of lambdas or anti-lambdas //
    TH1D* hist_cut_flow = anti_channel ? fHist_CutFlow_AntiChannel.get() : fHist_CutFlow.get();
    const auto decay_pid = HD::GetDecayTree(anti_channel);

    // alias input //
    const auto& input_lambdas = fInput.OnTheFlyLambda;
    const std::size_t n_lambdas = input_lambdas.size();

    // loop over all possible pairs of (anti)lambdas //
    for (std::size_t entry_lambda1 = 0; entry_lambda1 + 1 < n_lambdas; ++entry_lambda1) {
        const POD::OnTheFlyLambda& lambda1 = input_lambdas[entry_lambda1];  // cache index lookup

        /* lambda1.CacheCalculations(entry_lambda1, fPrimaryVertex, pid_neg, pid_pos); PENDING */
        auto lv_lambda1_neg =
            CMath::CreateLorentzVector(lambda1.Neg_PCAwrtV0_Px, lambda1.Neg_PCAwrtV0_Py, lambda1.Neg_PCAwrtV0_Pz, decay_pid.neg.mass);
        auto lv_lambda1_pos =
            CMath::CreateLorentzVector(lambda1.Pos_PCAwrtV0_Px, lambda1.Pos_PCAwrtV0_Py, lambda1.Pos_PCAwrtV0_Pz, decay_pid.pos.mass);
        auto lv_lambda1 = lv_lambda1_neg + lv_lambda1_pos;

        for (std::size_t entry_lambda2 = entry_lambda1 + 1; entry_lambda2 < n_lambdas; ++entry_lambda2) {
            const POD::OnTheFlyLambda& lambda2 = input_lambdas[entry_lambda2];  // cache index lookup
            // NOTE: sanity check not needed, because loops don't intersect

            /* lambda2.CacheCalculations(entry_lambda2, fPrimaryVertex, pid_neg, pid_pos); PENDING */
            auto lv_lambda2_neg =
                CMath::CreateLorentzVector(lambda2.Neg_PCAwrtV0_Px, lambda2.Neg_PCAwrtV0_Py, lambda2.Neg_PCAwrtV0_Pz, decay_pid.neg.mass);
            auto lv_lambda2_pos =
                CMath::CreateLorentzVector(lambda2.Pos_PCAwrtV0_Px, lambda2.Pos_PCAwrtV0_Py, lambda2.Pos_PCAwrtV0_Pz, decay_pid.pos.mass);
            auto lv_lambda2 = lv_lambda2_neg + lv_lambda2_pos;

            // PCAs //
            auto [seed_lambda1, seed_lambda2] = Seeder::LineLine::FullPCAs(lambda1, lambda2);

            // create composite particle //
            KF::LambdaPair kf_hdib(lambda1, lambda2, seed_lambda1.pca, seed_lambda2.pca, lv_lambda1.E(), lv_lambda2.E());

            // apply cuts //
            if (!Cuts(kf_hdib, hist_cut_flow)) continue;

            // store reconstructed //
            fOutput.Hdibaryon.emplace_back(CreateLambdaPair(kf_hdib, anti_channel));
            fOutput.Lambda1.emplace_back(lambda1, lv_lambda1.E());
            fOutput.Lambda2.emplace_back(lambda2, lv_lambda2.E());

            // store mc //
            if (fSettings.IsMC) {
                // -- lambda 1's neg
                auto mc_hdib_l1_neg = BuildMcTrack(fInput.OnTheFlyLambda_Neg_McEntry[entry_lambda1], decay_pid, decay_pid.neg.pdg_code);
                fOutput.MC_Lambda1_Neg.emplace_back(mc_hdib_l1_neg);
                // -- lambda 1's pos
                auto mc_hdib_l1_pos = BuildMcTrack(fInput.OnTheFlyLambda_Pos_McEntry[entry_lambda1], decay_pid, decay_pid.pos.pdg_code);
                fOutput.MC_Lambda1_Pos.emplace_back(mc_hdib_l1_pos);
                // -- lambda 1
                auto mc_hdib_l1 = BuildMcOnTheFlyLambda(mc_hdib_l1_neg, mc_hdib_l1_pos, decay_pid);
                fOutput.MC_Lambda1.emplace_back(mc_hdib_l1);
                // -- lambda 2's neg
                auto mc_hdib_l2_neg = BuildMcTrack(fInput.OnTheFlyLambda_Neg_McEntry[entry_lambda2], decay_pid, decay_pid.neg.pdg_code);
                fOutput.MC_Lambda2_Neg.emplace_back(mc_hdib_l2_neg);
                // -- lambda 2's pos
                auto mc_hdib_l2_pos = BuildMcTrack(fInput.OnTheFlyLambda_Pos_McEntry[entry_lambda2], decay_pid, decay_pid.pos.pdg_code);
                fOutput.MC_Lambda2_Pos.emplace_back(mc_hdib_l2_pos);
                // -- lambda 2
                auto mc_hdib_l2 = BuildMcOnTheFlyLambda(mc_hdib_l2_neg, mc_hdib_l2_pos, decay_pid);
                fOutput.MC_Lambda2.emplace_back(mc_hdib_l2);
                // -- h-dibaryon
                fOutput.MC_Hdibaryon.emplace_back(BuildMcHdibaryon(mc_hdib_l1, mc_hdib_l2, decay_pid));
            }
        }  // end of loop over pos
    }  // end of loop over neg
}

bool Verifier::Cuts(const KF::LambdaPair& kf_hdib, TH1D* cut_flow_hist) const {
    cut_flow_hist->Fill(0.);

    if (Common::Math::SquaredDistance(kf_hdib.Lambda1_PCAwrtDV->GetXYZ_AsROOT(), kf_hdib.Lambda2_PCAwrtDV->GetXYZ_AsROOT()) >
        Cuts::LambdaPair::Max_DCAbtwDau * Cuts::LambdaPair::Max_DCAbtwDau)
        return false;
    cut_flow_hist->Fill(1.);

    auto mass = kf_hdib.M();  // cached
    if (mass < Cuts::LambdaPair::Min_Mass || mass > Cuts::LambdaPair::Max_Mass) return false;
    cut_flow_hist->Fill(2.);

    return true;
}

POD::Extended::McParticle Verifier::BuildMcHdibaryon(const POD::Extended::McParticle& mc_lambda1, const POD::Extended::McParticle& mc_lambda2,
                                                     const HD::DecayTree& decay_pid) {
    POD::Extended::McParticle mc_hdib;
    // clang-format off
    auto mc_entry = MC::FindCommonMotherMcEntry(mc_lambda1, mc_lambda2);
    if (!mc_entry.has_value()) { return mc_hdib; }
    // clang-format on
    // fill values //
    static_cast<POD::McParticle&>(mc_hdib) = fInput.McParticle[mc_entry.value()];
    mc_hdib.Decay_X = mc_lambda1.Origin_X;
    mc_hdib.Decay_Y = mc_lambda1.Origin_Y;
    mc_hdib.Decay_Z = mc_lambda1.Origin_Z;
    mc_hdib.SignalID = static_cast<int>(mc_hdib.StatusCode);
    if (mc_hdib.Mother_McEntry > Common::DummyInt) {
        auto& mc_mother = fInput.McParticle[static_cast<std::size_t>(mc_hdib.Mother_McEntry)];
        mc_hdib.Mother_PdgCode = mc_mother.PdgCode;
        if (mc_mother.Mother_McEntry > Common::DummyInt) {
            mc_hdib.GM_McEntry = mc_mother.Mother_McEntry;
            mc_hdib.GM_PdgCode = fInput.McParticle[static_cast<std::size_t>(mc_hdib.GM_McEntry)].PdgCode;
        }
    }
    mc_hdib.IsTrue = mc_hdib.PdgCode == decay_pid.hdibaryon.pdg_code;
    mc_hdib.IsGen1Signal = false;
    mc_hdib.IsGen2Signal = false;
    mc_hdib.IsTrueSignal = mc_hdib.IsTrue;
    mc_hdib.IsSecondary = mc_hdib.IsSecFromMat || mc_hdib.IsSecFromWeak;
    mc_hdib.IsHybrid = mc_lambda1.IsHybrid || mc_lambda2.IsHybrid || mc_lambda1.IsTrueSignal != mc_lambda2.IsTrueSignal;

    return mc_hdib;
}

POD::LambdaPair Verifier::CreateLambdaPair(const KF::LambdaPair& kf_hdib, bool anti_channel) {
    POD::LambdaPair hdib;
    hdib.Decay_X = static_cast<float>(kf_hdib.DV.X());
    hdib.Decay_Y = static_cast<float>(kf_hdib.DV.Y());
    hdib.Decay_Z = static_cast<float>(kf_hdib.DV.Z());
    hdib.Px = static_cast<float>(kf_hdib.Px());
    hdib.Py = static_cast<float>(kf_hdib.Py());
    hdib.Pz = static_cast<float>(kf_hdib.Pz());
    hdib.Energy = static_cast<float>(kf_hdib.E());
    hdib.AntiChannel = anti_channel;
    // -- lambda1
    hdib.Lambda1_PCAwrtDV_X = static_cast<float>(kf_hdib.Lambda1_PCAwrtDV->X());
    hdib.Lambda1_PCAwrtDV_Y = static_cast<float>(kf_hdib.Lambda1_PCAwrtDV->Y());
    hdib.Lambda1_PCAwrtDV_Z = static_cast<float>(kf_hdib.Lambda1_PCAwrtDV->Z());
    hdib.Lambda1_PCAwrtDV_Px = static_cast<float>(kf_hdib.Lambda1_PCAwrtDV->Px());
    hdib.Lambda1_PCAwrtDV_Py = static_cast<float>(kf_hdib.Lambda1_PCAwrtDV->Py());
    hdib.Lambda1_PCAwrtDV_Pz = static_cast<float>(kf_hdib.Lambda1_PCAwrtDV->Pz());
    // -- lambda2
    hdib.Lambda2_PCAwrtDV_X = static_cast<float>(kf_hdib.Lambda2_PCAwrtDV->X());
    hdib.Lambda2_PCAwrtDV_Y = static_cast<float>(kf_hdib.Lambda2_PCAwrtDV->Y());
    hdib.Lambda2_PCAwrtDV_Z = static_cast<float>(kf_hdib.Lambda2_PCAwrtDV->Z());
    hdib.Lambda2_PCAwrtDV_Px = static_cast<float>(kf_hdib.Lambda2_PCAwrtDV->Px());
    hdib.Lambda2_PCAwrtDV_Py = static_cast<float>(kf_hdib.Lambda2_PCAwrtDV->Py());
    hdib.Lambda2_PCAwrtDV_Pz = static_cast<float>(kf_hdib.Lambda2_PCAwrtDV->Pz());

    return hdib;
}

// ## END OF CYCLES ## //

// Only fill events that have h-dibaryon candidates.
void Verifier::EndOfEvent() {
    // if data, don't keep event with no candidates
    if (!fSettings.IsMC && fOutput.Hdibaryon.empty()) return;
    // if mc, keep event with injected or reconstructed candidates
    if (fSettings.IsMC && fOutput.Hdibaryon.empty() && fOutput.Injected.empty()) return;
    // fill schema
    fWriter->Fill();
    // clear schema
    fOutput.Clear(fSettings.IsMC);
}

void Verifier::EndOfAnalysis() {

    Logger::Info(__FUNCTION__, "The following objects have been written into TFile {}:", fSettings.PathOutputFile);

    Logger::Info(__FUNCTION__, "- RNTuple \"{}\"", R2DS::Name_FoundHdibaryonRNT);

    // write histograms //

    fOutput_File->cd();

    // -- event counter
    fHist_EventCounter->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_EventCounter->GetName());

    // -- cut flows
    fHist_CutFlow->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow->GetName());
    fHist_CutFlow_AntiChannel->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow_AntiChannel->GetName());

    Logger::Info(__FUNCTION__, "All done.");
}

}  // namespace R2DS
