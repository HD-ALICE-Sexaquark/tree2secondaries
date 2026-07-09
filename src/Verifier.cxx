#include <cstddef>
#include <memory>
#include <tuple>

#include "common/Cached_Hdibaryon.hpp"
#include "common/Cached_PreFoundLambda.hpp"
#include "common/Constants.hpp"
#include "common/DB_Particles.hpp"
#include "common/HD_Library.hpp"
#include "common/MC_Helpers.hpp"
#include "common/POD_InjectedHdib.hpp"
#include "common/POD_LambdaPair.hpp"
#include "common/POD_McParticle.hpp"
#include "common/POD_PreFoundLambda.hpp"

#include "common/Math.hpp"
namespace CMath = Common::Math;

#include "KalmanFitter/BaseKalmanFitter.hxx"
#include "KalmanFitter/KalmanFitterParticle.hxx"
#include "Seeder/BaseSeeder.hxx"
#include "Seeder/SeederHelixHelix.hxx"
#include "Seeder/SeederLineLine.hxx"

#include "Verifier/Verifier.hxx"

namespace T2DS {

// ## OUTPUT ZONE ## //

void Verifier::PrepareOutputHistograms() {
    // event counter
    fHist_EventCounter = std::make_unique<TH1D>("N_Events", ";N_Events;", 1, 0., 1.);
    // cut flows
    constexpr int x_nbins = 20;
    constexpr float x_min = 0.;
    constexpr float x_max = 20.;
    constexpr const char* hist_title = ";Cut N;N Passed Cut";
    // -- for (anti)lambdas
    fHist_CutFlow_AntiLambda = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", DB::Particles::Particle("AntiLambda").acronym).c_str(), hist_title, x_nbins, x_min, x_max);
    fHist_CutFlow_Lambda = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", DB::Particles::Particle("Lambda").acronym).c_str(), hist_title, x_nbins, x_min, x_max);
    // -- for (anti)h-dibaryons
    fHist_CutFlow_AntiHdibaryon = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", DB::Particles::Particle("AntiHdibaryon").acronym).c_str(), hist_title, x_nbins, x_min, x_max);
    fHist_CutFlow_Hdibaryon = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", DB::Particles::Particle("Hdibaryon").acronym).c_str(), hist_title, x_nbins, x_min, x_max);
}

// ## Event ZONE ## //

void Verifier::ProcessEvent() {
    // copy event info
    fOutput.Event = static_cast<POD::Event&>(fInput.Event);
    if (fSettings.IsMC) {
        fOutput.MC_Event = fInput.MC_Event;
    }
    // update event counter
    fHist_EventCounter->Fill(0.);
    // cache pv
    fPrimaryVertex.SetCoordinates(fOutput.Event.PV_X, fOutput.Event.PV_Y, fOutput.Event.PV_Z);
}

// ## Injected ZONE ## //

void Verifier::ProcessInjected() {
    // loop over mc particles //
    for (const auto& mc : fInput.McParticle) {
        // select only injected h-dibaryons //
        if (std::abs(mc.PdgCode) != DB::Particles::Particle("Hdibaryon").pdg_code) continue;
        // store //
        fOutput.Injected.emplace_back(BuildInjectedHdibaryon(mc));
    }
}

POD::InjectedHdib Verifier::BuildInjectedHdibaryon(const POD::McParticle& mc) {
    POD::InjectedHdib inj;
    inj.SignalID = static_cast<int>(mc.StatusCode);
    inj.IsAntiChannel = mc.PdgCode == DB::Particles::Particle("AntiHdibaryon").pdg_code;
    std::tie(inj.Decay_X, inj.Decay_Y, inj.Decay_Z) = MC::GetDecayVertex(mc, fInput.McParticle);
    inj.Px = mc.Px;
    inj.Py = mc.Py;
    inj.Pz = mc.Pz;
    inj.Energy = mc.Energy;
    // lambda 1
    if (mc.FirstDau_McEntry > 0) {
        auto& mc_l1 = fInput.McParticle[static_cast<std::size_t>(mc.FirstDau_McEntry)];
        std::tie(inj.Lambda1_Decay_X, inj.Lambda1_Decay_Y, inj.Lambda1_Decay_Z) = MC::GetDecayVertex(mc_l1, fInput.McParticle);
        inj.Lambda1_Px = mc_l1.Px;
        inj.Lambda1_Py = mc_l1.Py;
        inj.Lambda1_Pz = mc_l1.Pz;
        inj.Lambda1_Energy = mc_l1.Energy;
        // lambda 1's proton daughter
        auto entry_proton =
            MC::FindDaughterMcEntry(mc_l1, fInput.McParticle,
                                    inj.IsAntiChannel ? DB::Particles::Particle("AntiProton").pdg_code : DB::Particles::Particle("Proton").pdg_code);
        if (entry_proton.has_value()) {
            auto& mc_proton = fInput.McParticle[entry_proton.value()];
            inj.Proton1_Px = mc_proton.Px;
            inj.Proton1_Py = mc_proton.Py;
            inj.Proton1_Pz = mc_proton.Pz;
            inj.Proton1_Energy = mc_proton.Energy;
        }
        // lambda 1's pion daughter
        auto entry_pion = MC::FindDaughterMcEntry(
            mc_l1, fInput.McParticle, inj.IsAntiChannel ? DB::Particles::Particle("PiPlus").pdg_code : DB::Particles::Particle("PiMinus").pdg_code);
        if (entry_pion.has_value()) {
            auto& mc_pion = fInput.McParticle[entry_pion.value()];
            inj.Pion1_Px = mc_pion.Px;
            inj.Pion1_Py = mc_pion.Py;
            inj.Pion1_Pz = mc_pion.Pz;
            inj.Pion1_Energy = mc_pion.Energy;
        }
    }
    // lambda 2
    if (mc.LastDau_McEntry > 0) {
        auto& mc_l2 = fInput.McParticle[static_cast<std::size_t>(mc.LastDau_McEntry)];
        std::tie(inj.Lambda2_Decay_X, inj.Lambda2_Decay_Y, inj.Lambda2_Decay_Z) = MC::GetDecayVertex(mc_l2, fInput.McParticle);
        inj.Lambda2_Px = mc_l2.Px;
        inj.Lambda2_Py = mc_l2.Py;
        inj.Lambda2_Pz = mc_l2.Pz;
        inj.Lambda2_Energy = mc_l2.Energy;
        // lambda 2's proton daughter
        auto entry_proton =
            MC::FindDaughterMcEntry(mc_l2, fInput.McParticle,
                                    inj.IsAntiChannel ? DB::Particles::Particle("AntiProton").pdg_code : DB::Particles::Particle("Proton").pdg_code);
        if (entry_proton.has_value()) {
            auto& mc_proton = fInput.McParticle[entry_proton.value()];
            inj.Proton2_Px = mc_proton.Px;
            inj.Proton2_Py = mc_proton.Py;
            inj.Proton2_Pz = mc_proton.Pz;
            inj.Proton2_Energy = mc_proton.Energy;
        }
        // lambda 2's pion daughter
        auto entry_pion = MC::FindDaughterMcEntry(
            mc_l2, fInput.McParticle, inj.IsAntiChannel ? DB::Particles::Particle("PiPlus").pdg_code : DB::Particles::Particle("PiMinus").pdg_code);
        if (entry_pion.has_value()) {
            auto& mc_pion = fInput.McParticle[entry_pion.value()];
            inj.Pion2_Px = mc_pion.Px;
            inj.Pion2_Py = mc_pion.Py;
            inj.Pion2_Pz = mc_pion.Pz;
            inj.Pion2_Energy = mc_pion.Energy;
        }
    }

    return inj;
}

// ## Charged Tracks ZONE ## //

POD::Extended::McParticle Verifier::BuildMcTrack(unsigned int track_mc_entry, const HD::DecayTree& decay_pid, int pdg_code_hypothesis) {
    // copy linked mc info //
    POD::Extended::McParticle new_mc(fInput.McParticle[track_mc_entry]);
    MC::Apply(new_mc, MC::HdibaryonRules::ClassifyDownstream(new_mc, fInput.McParticle, decay_pid, pdg_code_hypothesis, false));
    return new_mc;
}

POD::Track Verifier::ExtractTrack(const POD::PreFoundLambda& pod_lambda, short charge) const {
    POD::Track new_track;
    new_track.EsdEntry = charge < 0 ? pod_lambda.Neg_EsdEntry : pod_lambda.Pos_EsdEntry;
    new_track.X = charge < 0 ? pod_lambda.Neg_State[0] : pod_lambda.Pos_State[0];
    new_track.Y = charge < 0 ? pod_lambda.Neg_State[1] : pod_lambda.Pos_State[1];
    new_track.Z = charge < 0 ? pod_lambda.Neg_State[2] : pod_lambda.Pos_State[2];
    new_track.Px = charge < 0 ? pod_lambda.Neg_State[3] : pod_lambda.Pos_State[3];
    new_track.Py = charge < 0 ? pod_lambda.Neg_State[4] : pod_lambda.Pos_State[4];
    new_track.Pz = charge < 0 ? pod_lambda.Neg_State[5] : pod_lambda.Pos_State[5];
    new_track.Charge = charge;
    new_track.PreDCAxy = charge < 0 ? pod_lambda.Neg_PreDCAxy : pod_lambda.Pos_PreDCAxy;
    new_track.PreDCAz = charge < 0 ? pod_lambda.Neg_PreDCAz : pod_lambda.Pos_PreDCAz;
    new_track.TPC_Signal = Common::DummyFloat;
    new_track.NSigmasPion = charge < 0 ? pod_lambda.Neg_NSigmasPion : pod_lambda.Pos_NSigmasPion;
    new_track.NSigmasKaon = charge < 0 ? pod_lambda.Neg_NSigmasKaon : pod_lambda.Pos_NSigmasKaon;
    new_track.NSigmasProton = charge < 0 ? pod_lambda.Neg_NSigmasProton : pod_lambda.Pos_NSigmasProton;
    new_track.CovMatrix = charge < 0 ? pod_lambda.Neg_CovMatrix : pod_lambda.Pos_CovMatrix;
    new_track.TPC_FirstRow = Common::DummyInt;

    return new_track;
}

// ## Single Pre-Found Lambda ZONE ## //

void Verifier::ProcessPreFoundLambda() {

    // loop over pre-found (anti)lambdas //
    for (std::size_t entry_lambda = 0; entry_lambda < fInput.PreFoundLambda.size(); ++entry_lambda) {

        const auto& in_lambda = fInput.PreFoundLambda[entry_lambda];

        auto track_neg = ExtractTrack(in_lambda, -1);
        auto track_pos = ExtractTrack(in_lambda, +1);

        // PCAs //
        auto [seed_neg, seed_pos, pca_cache] = Seeder::HelixHelix::FastCorrectPCAs(track_neg, track_pos, fOutput.Event.MagneticField);

        // apply cuts (1) //
        if (!FastCuts_Lambda(seed_neg.pca, seed_pos.pca)) continue;

        // PCAs derivatives //
        auto [deriv_neg, deriv_pos] = Seeder::HelixHelix::ComputeDerivatives(seed_neg, seed_pos, pca_cache);

        // loop over lambda / anti-lambda hypotheses //
        for (auto anti_channel : {false, true}) {

            const auto& decay_tree = HD::GetDecayTree(anti_channel);
            auto* hist_cut_flow = anti_channel ? fHist_CutFlow_AntiLambda.get() : fHist_CutFlow_Lambda.get();

            // fit vertex //
            auto fit = KF::FitVertex(track_neg, track_pos, decay_tree.neg.mass, decay_tree.pos.mass, {seed_neg, deriv_neg}, {seed_pos, deriv_pos},
                                     fOutput.Event.MagneticField);

            // create storage+computation (anti)lambda //
            POD::Extended::PreFoundLambda new_lambda =
                CreateExtendedPreFoundLambda(in_lambda, fit, seed_neg.pca, seed_pos.pca, decay_tree.neg.mass, decay_tree.pos.mass);
            Cached::PreFoundLambda c_lambda(new_lambda, decay_tree.neg.mass, decay_tree.pos.mass, fPrimaryVertex);

            // apply more cuts (2) //
            if (!SlowCuts_Lambda(c_lambda, hist_cut_flow)) continue;

            // store reconstructed //
            if (anti_channel) {
                fTemp_AntiLambda.emplace_back(new_lambda);
            } else {
                fTemp_Lambda.emplace_back(new_lambda);
            }

            // store mc //
            if (fSettings.IsMC) {
                // -- build mc particles
                auto mc_lambda_neg = BuildMcTrack(fInput.PreFoundLambda_Neg_McEntry[entry_lambda], decay_tree, decay_tree.neg.pdg_code);
                auto mc_lambda_pos = BuildMcTrack(fInput.PreFoundLambda_Pos_McEntry[entry_lambda], decay_tree, decay_tree.pos.pdg_code);
                auto mc_lambda = BuildMcPreFoundLambda(mc_lambda_neg, mc_lambda_pos, decay_tree);
                // -- store
                if (anti_channel) {
                    fTemp_MC_AntiLambda.emplace_back(mc_lambda);
                    fTemp_MC_AntiLambda_Neg.emplace_back(mc_lambda_neg);
                    fTemp_MC_AntiLambda_Pos.emplace_back(mc_lambda_pos);
                } else {
                    fTemp_MC_Lambda.emplace_back(mc_lambda);
                    fTemp_MC_Lambda_Neg.emplace_back(mc_lambda_neg);
                    fTemp_MC_Lambda_Pos.emplace_back(mc_lambda_pos);
                }
            }
        }
    }  // end of loop over pre-found (anti)lambdas
}

bool Verifier::FastCuts_Lambda(const Seeder::PCA& pca_neg, const Seeder::PCA& pca_pos) {
    // NOTE: yes, fill both histograms at this stage.
    fHist_CutFlow_AntiLambda->Fill(0.);
    fHist_CutFlow_Lambda->Fill(0.);

    // PENDING //

    return true;
}

bool Verifier::SlowCuts_Lambda(const Cached::PreFoundLambda& lambda, TH1D* hist_cut_flow) {
    hist_cut_flow->Fill(1.);

    // PENDING //

    return true;
}

POD::Extended::McParticle Verifier::BuildMcPreFoundLambda(const POD::Extended::McParticle& mc_neg, const POD::Extended::McParticle& mc_pos,
                                                          const HD::DecayTree& decay_pid) {
    POD::Extended::McParticle mc_lambda;

    // -- fill hybridness, independently of no common mother
    mc_lambda.IsHybrid = mc_neg.IsTrueSignal != mc_pos.IsTrueSignal;

    auto mc_entry = MC::FindCommonMotherMcEntry(mc_neg, mc_pos);
    if (!mc_entry.has_value()) return mc_lambda;

    // fill rest of values //
    static_cast<POD::McParticle&>(mc_lambda) = fInput.McParticle[mc_entry.value()];
    MC::Apply(mc_lambda, MC::HdibaryonRules::ClassifyDownstream(mc_lambda, fInput.McParticle, decay_pid, decay_pid.lambda.pdg_code, true));

    return mc_lambda;
}

POD::Extended::PreFoundLambda Verifier::CreateExtendedPreFoundLambda(const POD::PreFoundLambda& old_lambda, const KF::Particle& fit,
                                                                     const Seeder::PCA& pca_neg, const Seeder::PCA& pca_pos, double mass_neg,
                                                                     double mass_pos) {
    // profit from initialization list to extend data
    POD::Extended::PreFoundLambda new_lambda{old_lambda,
                                             static_cast<float>(fit.Px()),
                                             static_cast<float>(fit.Py()),
                                             static_cast<float>(fit.Pz()),
                                             static_cast<float>(fit.E()),
                                             fit.Cov<7>(),
                                             static_cast<float>(fit.Chi2NDF()),
                                             static_cast<float>(pca_neg.X()),
                                             static_cast<float>(pca_neg.Y()),
                                             static_cast<float>(pca_neg.Z()),
                                             static_cast<float>(CMath::Hypot4(pca_neg.Px(), pca_neg.Py(), pca_neg.Pz(), mass_neg)),
                                             static_cast<float>(pca_pos.X()),
                                             static_cast<float>(pca_pos.Y()),
                                             static_cast<float>(pca_pos.Z()),
                                             static_cast<float>(CMath::Hypot4(pca_pos.Px(), pca_pos.Py(), pca_pos.Pz(), mass_pos))};
    // override previous info
    new_lambda.Decay_X = static_cast<float>(fit.X());
    new_lambda.Decay_Y = static_cast<float>(fit.Y());
    new_lambda.Decay_Z = static_cast<float>(fit.Z());
    new_lambda.DcaV0Daughters = static_cast<float>(CMath::Distance(pca_neg.xyz, pca_pos.xyz));
    // -- negative daughter
    new_lambda.Neg_PCAwrtV0_Px = static_cast<float>(pca_neg.Px());
    new_lambda.Neg_PCAwrtV0_Py = static_cast<float>(pca_neg.Py());
    new_lambda.Neg_PCAwrtV0_Pz = static_cast<float>(pca_neg.Pz());
    // -- positive daughter
    new_lambda.Pos_PCAwrtV0_Px = static_cast<float>(pca_pos.Px());
    new_lambda.Pos_PCAwrtV0_Py = static_cast<float>(pca_pos.Py());
    new_lambda.Pos_PCAwrtV0_Pz = static_cast<float>(pca_pos.Pz());

    return new_lambda;
}

// ## H-dibaryon ZONE ## //

void Verifier::VerifyLambdaPair(bool anti_channel) {

    // determine rules based on reconstruction of lambdas or anti-lambdas //
    const auto& decay_pid = HD::GetDecayTree(anti_channel);
    auto* hist_cut_flow = anti_channel ? fHist_CutFlow_AntiHdibaryon.get() : fHist_CutFlow_Hdibaryon.get();
    const auto& input_lambdas = anti_channel ? fTemp_AntiLambda : fTemp_Lambda;
    const auto& input_mc_lambdas = anti_channel ? fTemp_MC_AntiLambda : fTemp_MC_Lambda;
    const auto& input_mc_lambdas_neg = anti_channel ? fTemp_MC_AntiLambda_Neg : fTemp_MC_Lambda_Neg;
    const auto& input_mc_lambdas_pos = anti_channel ? fTemp_MC_AntiLambda_Pos : fTemp_MC_Lambda_Pos;
    const std::size_t n_lambdas = input_lambdas.size();

    // loop over all possible pairs of (anti)lambdas //
    for (std::size_t entry_lambda1 = 0; entry_lambda1 + 1 < n_lambdas; ++entry_lambda1) {
        const auto& lambda1 = input_lambdas[entry_lambda1];  // cache index lookup

        for (std::size_t entry_lambda2 = entry_lambda1 + 1; entry_lambda2 < n_lambdas; ++entry_lambda2) {
            const auto& lambda2 = input_lambdas[entry_lambda2];  // cache index lookup
            // NOTE: sanity check not needed, because loops don't intersect

            // PCAs //
            Seeder::LineLine::Cache pca_cache;
            auto [seed_lambda1, seed_lambda2] = Seeder::LineLine::FastPCAs(lambda1, lambda2, &pca_cache);

            // apply cuts (1) //
            if (!FastCuts_Hdibaryon(seed_lambda1.pca, seed_lambda2.pca, hist_cut_flow)) continue;

            // PCAs derivatives //
            auto [deriv_lambda1, deriv_lambda2] = Seeder::LineLine::ComputeDerivatives(pca_cache);

            // fit vertex //
            auto fit = KF::FitVertex(lambda1, lambda2, {seed_lambda1, deriv_lambda1}, {seed_lambda2, deriv_lambda2});

            // create storage+computation units //
            POD::LambdaPair hdib = CreateLambdaPair(fit, seed_lambda1.pca, seed_lambda2.pca, anti_channel);
            Cached::Hdibaryon c_hdib(hdib, lambda1, lambda2, fPrimaryVertex);

            // apply cuts (2) //
            if (!SlowCuts_Hdibaryon(c_hdib, hist_cut_flow)) continue;

            // store reconstructed //
            fOutput.Hdibaryon.emplace_back(hdib);
            fOutput.Lambda1.emplace_back(lambda1);
            fOutput.Lambda2.emplace_back(lambda2);

            // store mc //
            if (fSettings.IsMC) {
                fOutput.MC_Hdibaryon.emplace_back(BuildMcHdibaryon(input_mc_lambdas[entry_lambda1], input_mc_lambdas[entry_lambda2], decay_pid));
                fOutput.MC_Lambda1.emplace_back(input_mc_lambdas[entry_lambda1]);
                fOutput.MC_Lambda1_Neg.emplace_back(input_mc_lambdas_neg[entry_lambda1]);
                fOutput.MC_Lambda1_Pos.emplace_back(input_mc_lambdas_pos[entry_lambda1]);
                fOutput.MC_Lambda2.emplace_back(input_mc_lambdas[entry_lambda2]);
                fOutput.MC_Lambda2_Neg.emplace_back(input_mc_lambdas_neg[entry_lambda2]);
                fOutput.MC_Lambda2_Pos.emplace_back(input_mc_lambdas_pos[entry_lambda2]);
            }
        }  // end of loop over lambda2
    }  // end of loop over lambda1
}

bool Verifier::FastCuts_Hdibaryon(const Seeder::PCA& pca_lambda1, const Seeder::PCA& pca_lambda2, TH1D* hist_cut_flow) {
    hist_cut_flow->Fill(0.);

    // PENDING //

    return true;
}

bool Verifier::SlowCuts_Hdibaryon(const Cached::Hdibaryon& c_hdib, TH1D* hist_cut_flow) {
    hist_cut_flow->Fill(1.);

    // PENDING //

    return true;
}

POD::Extended::McParticle Verifier::BuildMcHdibaryon(const POD::Extended::McParticle& mc_lambda1, const POD::Extended::McParticle& mc_lambda2,
                                                     const HD::DecayTree& decay_pid) {
    POD::Extended::McParticle mc_hdib;

    // -- fill hybridness, independently of no common mother
    mc_hdib.IsHybrid = mc_lambda1.IsHybrid || mc_lambda2.IsHybrid || mc_lambda1.IsTrueSignal != mc_lambda2.IsTrueSignal;

    auto mc_entry = MC::FindCommonMotherMcEntry(mc_lambda1, mc_lambda2);
    if (!mc_entry.has_value()) return mc_hdib;

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

    return mc_hdib;
}

POD::LambdaPair Verifier::CreateLambdaPair(const KF::Particle& fit, const Seeder::PCA& pca_lambda1, const Seeder::PCA& pca_lambda2,
                                           bool anti_channel) {
    POD::LambdaPair new_hdib;
    new_hdib.Decay_X = static_cast<float>(fit.X());
    new_hdib.Decay_Y = static_cast<float>(fit.Y());
    new_hdib.Decay_Z = static_cast<float>(fit.Z());
    new_hdib.Px = static_cast<float>(fit.Px());
    new_hdib.Py = static_cast<float>(fit.Py());
    new_hdib.Pz = static_cast<float>(fit.Pz());
    new_hdib.Energy = static_cast<float>(fit.E());
    new_hdib.Chi2NDF = static_cast<float>(fit.Chi2NDF());
    new_hdib.IsAntiChannel = anti_channel;
    // -- lambda1
    new_hdib.Lambda1_PCAwrtDV_X = static_cast<float>(pca_lambda1.X());
    new_hdib.Lambda1_PCAwrtDV_Y = static_cast<float>(pca_lambda1.Y());
    new_hdib.Lambda1_PCAwrtDV_Z = static_cast<float>(pca_lambda1.Z());
    new_hdib.Lambda1_PCAwrtDV_Px = static_cast<float>(pca_lambda1.Px());
    new_hdib.Lambda1_PCAwrtDV_Py = static_cast<float>(pca_lambda1.Py());
    new_hdib.Lambda1_PCAwrtDV_Pz = static_cast<float>(pca_lambda1.Pz());
    // -- lambda2
    new_hdib.Lambda2_PCAwrtDV_X = static_cast<float>(pca_lambda2.X());
    new_hdib.Lambda2_PCAwrtDV_Y = static_cast<float>(pca_lambda2.Y());
    new_hdib.Lambda2_PCAwrtDV_Z = static_cast<float>(pca_lambda2.Z());
    new_hdib.Lambda2_PCAwrtDV_Px = static_cast<float>(pca_lambda2.Px());
    new_hdib.Lambda2_PCAwrtDV_Py = static_cast<float>(pca_lambda2.Py());
    new_hdib.Lambda2_PCAwrtDV_Pz = static_cast<float>(pca_lambda2.Pz());

    return new_hdib;
}

// ## END OF CYCLES ## //

// Only fill events that have h-dibaryon candidates.
void Verifier::EndOfEvent() {
    // clear temporary vectors
    fTemp_AntiLambda.clear();
    fTemp_Lambda.clear();
    if (fSettings.IsMC) {
        fTemp_MC_AntiLambda.clear();
        fTemp_MC_AntiLambda_Neg.clear();
        fTemp_MC_AntiLambda_Pos.clear();
        fTemp_MC_Lambda.clear();
        fTemp_MC_Lambda_Neg.clear();
        fTemp_MC_Lambda_Pos.clear();
    }

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

    Logger::Info(__FUNCTION__, "The following objects have been written into TFile \"{}\":", fSettings.PathOutputFile);

    Logger::Info(__FUNCTION__, "- RNTuple \"{}\"", T2DS::Name_FoundHdibaryonRNT);

    // write histograms //

    fOutput_File->cd();

    // -- event counter
    fHist_EventCounter->Write();
    Logger::Info(__FUNCTION__, "- TH1D \"{}\"", fHist_EventCounter->GetName());
    // -- cut flow for (anti)lambdas
    fHist_CutFlow_AntiLambda->Write();
    Logger::Info(__FUNCTION__, "- TH1D \"{}\"", fHist_CutFlow_AntiLambda->GetName());
    fHist_CutFlow_Lambda->Write();
    Logger::Info(__FUNCTION__, "- TH1D \"{}\"", fHist_CutFlow_Lambda->GetName());
    // -- cut flow for (anti)h-dibaryons
    fHist_CutFlow_AntiHdibaryon->Write();
    Logger::Info(__FUNCTION__, "- TH1D \"{}\"", fHist_CutFlow_AntiHdibaryon->GetName());
    fHist_CutFlow_Hdibaryon->Write();
    Logger::Info(__FUNCTION__, "- TH1D \"{}\"", fHist_CutFlow_Hdibaryon->GetName());

    Logger::Info(__FUNCTION__, "All done.");
}

}  // namespace T2DS
