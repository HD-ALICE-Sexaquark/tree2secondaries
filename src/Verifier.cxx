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
    // data kind must be MC //
    if (!fInput.McParticle.has_value()) return;
    const auto& mc_collection = fInput.McParticle.value();
    // init new injected vector //
    fOutput.InjectedHdib.emplace();
    // loop over mc particles //
    for (const auto& mc : mc_collection) {
        // select only injected h-dibaryons //
        if (std::abs(mc.PdgCode) != DB::Particles::Particle("Hdibaryon").pdg_code) continue;
        // create new injected h-dibaryon //
        POD::McParticle new_hdib = mc;
        // fill some extra info //
        new_hdib.SignalID = mc.Status;
        std::tie(new_hdib.Decay_X, new_hdib.Decay_Y, new_hdib.Decay_Z) = MC::GetDecayVertex(mc, mc_collection);
        // store //
        fOutput.InjectedHdib->emplace_back(new_hdib);
    }
}

// ## Single On-The-Fly Lambda ZONE ## //

void Verifier::BuildMcInfo(POD::OnTheFlyLambda& new_lambda, const HD::DecayTree& decay_pid) {
    // -- data kind must be MC
    if (!fInput.McParticle.has_value()) return;
    const auto& mc_collection = fInput.McParticle.value();  // alias
    // -- negative daughter
    if (new_lambda.Neg_McEntry.has_value()) {
        new_lambda.Neg_MC = mc_collection[new_lambda.Neg_McEntry.value()];
    }
    auto& mc_neg = new_lambda.Neg_MC.value();  // alias
    MC::Apply(mc_neg, MC::HdibaryonRules::Classify(mc_neg, mc_collection, decay_pid, decay_pid.neg.pdg_code, false));
    // -- positive daughter
    if (new_lambda.Pos_McEntry.has_value()) {
        new_lambda.Pos_MC = mc_collection[new_lambda.Pos_McEntry.value()];
    }
    auto& mc_pos = new_lambda.Pos_MC.value();  // alias
    MC::Apply(mc_pos, MC::HdibaryonRules::Classify(mc_pos, mc_collection, decay_pid, decay_pid.pos.pdg_code, false));
    // -- find V0's `McEntry`
    auto mc_entry = MC::FindCommonMotherMcEntry(mc_neg, mc_pos);
    if (!mc_entry.has_value()) return;
    // -- classify and store
    new_lambda.MC = mc_collection[mc_entry.value()];
    auto& mc_lambda = new_lambda.MC.value();  // alias
    MC::Apply(mc_lambda, MC::HdibaryonRules::Classify(mc_lambda, mc_collection, decay_pid, decay_pid.lambda.pdg_code, true));
    mc_lambda.IsHybrid = mc_neg.IsTrueSignal != mc_pos.IsTrueSignal;
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

            // store //
            POD::LambdaPair new_hdib;
            // -- reconstructed info
            BuildRecInfo(new_hdib, kf_hdib, anti_channel);
            // -- linked mc info
            BuildMcInfo(new_hdib, decay_pid);
            // fill //
            fOutput.LambdaPair.emplace_back(new_hdib);
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

void Verifier::BuildMcInfo(POD::LambdaPair& new_hdib, const HD::DecayTree& decay_pid) {
    // -- data kind must be MC
    if (!fInput.McParticle.has_value()) return;
    const auto& mc_collection = fInput.McParticle.value();  // alias
    // lambdas
    BuildMcInfo(new_hdib.Lambda1, decay_pid);
    BuildMcInfo(new_hdib.Lambda2, decay_pid);
    // injected
    if (!new_hdib.Lambda1.MC.has_value() || !new_hdib.Lambda2.MC.has_value()) return;
    const auto& mc_lambda1 = new_hdib.Lambda1.MC.value();  // alias
    const auto& mc_lambda2 = new_hdib.Lambda2.MC.value();  // alias
    // -- find common entry
    auto mc_entry = MC::FindCommonMotherMcEntry(mc_lambda1, mc_lambda2);
    if (!mc_entry.has_value()) return;
    // copy values
    new_hdib.MC = mc_collection[mc_entry.value()];
    auto& mc_hdib = new_hdib.MC.value();  // alias
    if (mc_hdib.PdgCode != decay_pid.hdibaryon.pdg_code) return;
    // fill some extra info //
    mc_hdib.SignalID = mc_hdib.Status;
    mc_hdib.Decay_X = mc_lambda1.Origin_X;
    mc_hdib.Decay_Y = mc_lambda1.Origin_Y;
    mc_hdib.Decay_Z = mc_lambda1.Origin_Z;
    new_hdib.IsHybrid =
        mc_lambda1.IsHybrid.value_or(false) || mc_lambda2.IsHybrid.value_or(false) || mc_lambda1.IsTrueSignal != mc_lambda2.IsTrueSignal;
}

void Verifier::BuildRecInfo(POD::LambdaPair& new_hdib, const KF::LambdaPair& kf_hdib, bool anti_channel) {
    // -- candidate info
    new_hdib.Decay_X = static_cast<float>(kf_hdib.DV.X());
    new_hdib.Decay_Y = static_cast<float>(kf_hdib.DV.Y());
    new_hdib.Decay_Z = static_cast<float>(kf_hdib.DV.Z());
    new_hdib.Px = static_cast<float>(kf_hdib.Px());
    new_hdib.Py = static_cast<float>(kf_hdib.Py());
    new_hdib.Pz = static_cast<float>(kf_hdib.Pz());
    new_hdib.Energy = static_cast<float>(kf_hdib.E());
    new_hdib.AntiChannel = anti_channel;
    // -- lambda1
    new_hdib.Lambda1 = *kf_hdib.Lambda1;
    new_hdib.Lambda1_PCAwrtDV_X = static_cast<float>(kf_hdib.Lambda1_PCAwrtDV->X());
    new_hdib.Lambda1_PCAwrtDV_Y = static_cast<float>(kf_hdib.Lambda1_PCAwrtDV->Y());
    new_hdib.Lambda1_PCAwrtDV_Z = static_cast<float>(kf_hdib.Lambda1_PCAwrtDV->Z());
    new_hdib.Lambda1_PCAwrtDV_Px = static_cast<float>(kf_hdib.Lambda1_PCAwrtDV->Px());
    new_hdib.Lambda1_PCAwrtDV_Py = static_cast<float>(kf_hdib.Lambda1_PCAwrtDV->Py());
    new_hdib.Lambda1_PCAwrtDV_Pz = static_cast<float>(kf_hdib.Lambda1_PCAwrtDV->Pz());
    // -- lambda2
    new_hdib.Lambda2 = *kf_hdib.Lambda2;
    new_hdib.Lambda2_PCAwrtDV_X = static_cast<float>(kf_hdib.Lambda2_PCAwrtDV->X());
    new_hdib.Lambda2_PCAwrtDV_Y = static_cast<float>(kf_hdib.Lambda2_PCAwrtDV->Y());
    new_hdib.Lambda2_PCAwrtDV_Z = static_cast<float>(kf_hdib.Lambda2_PCAwrtDV->Z());
    new_hdib.Lambda2_PCAwrtDV_Px = static_cast<float>(kf_hdib.Lambda2_PCAwrtDV->Px());
    new_hdib.Lambda2_PCAwrtDV_Py = static_cast<float>(kf_hdib.Lambda2_PCAwrtDV->Py());
    new_hdib.Lambda2_PCAwrtDV_Pz = static_cast<float>(kf_hdib.Lambda2_PCAwrtDV->Pz());
}

// ## END OF CYCLES ## //

void Verifier::EndOfEvent() {
    // fill schema
    fWriter.Fill();
    // clear schema vectors
    fOutput.Clear();
}

void Verifier::EndOfAnalysis() {

    Logger::Info(__FUNCTION__, "The following objects have been written into TFile {}:", fSettings.PathOutputFile);

    Logger::Info(__FUNCTION__, "- RNTuple \"{}\"", R2DS::Name_FoundHdibRNT);

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
