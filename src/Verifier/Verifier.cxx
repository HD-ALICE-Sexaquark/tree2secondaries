#include <cmath>
#include <cstddef>
#include <format>
#include <memory>
#include <optional>
#include <tuple>

#include "common/Cached_Hdibaryon.hpp"
#include "common/Cached_PreFoundLambda.hpp"
#include "common/Constants.hpp"
#include "common/Cuts_T2DS_Verifier.hpp"
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
#include "Seeder/SeederHelixHelix.hxx"
#include "Seeder/SeederLineLine.hxx"
#include "Seeder/SeederTypes.hxx"

#include "Verifier/Verifier.hxx"

namespace T2DS {

// ## OUTPUT ZONE ## //

void Verifier::PrepareOutputHistograms() {
    // event counter
    fHist_EventCounter = std::make_unique<TH1D>("N_Events", ";N_Events;", 1, 0., 1.);

    // cut flows
    constexpr const char* hist_title = ";;N Passed Cut";

    // -- for (anti)lambdas
    fHist_CutFlow_AntiLambda = std::make_unique<TH1D>(                                                 //
        std::format("CutFlow_{}", DB::Particles::Particle("AntiLambda").acronym).c_str(), hist_title,  //
        static_cast<int>(EPreFoundLambda::kNPreFoundLambdaCuts), 0., static_cast<double>(EPreFoundLambda::kNPreFoundLambdaCuts));
    fHist_CutFlow_Lambda = std::make_unique<TH1D>(                                                 //
        std::format("CutFlow_{}", DB::Particles::Particle("Lambda").acronym).c_str(), hist_title,  //
        static_cast<int>(EPreFoundLambda::kNPreFoundLambdaCuts), 0., static_cast<double>(EPreFoundLambda::kNPreFoundLambdaCuts));

    // -- define bin labels
    for (auto* hist_lambda : {fHist_CutFlow_AntiLambda.get(), fHist_CutFlow_Lambda.get()}) {
        auto* x_axis = hist_lambda->GetXaxis();
        x_axis->SetBinLabel(static_cast<int>(EPreFoundLambda::kAllPreFoundLambdas) + 1, "AllPreFoundLambdas");
        // pre-seed cuts
        x_axis->SetBinLabel(static_cast<int>(EPreFoundLambda::kPasses_DiffDaughters_Logical) + 1, "Passes_DiffDaughters_Logical");
        x_axis->SetBinLabel(static_cast<int>(EPreFoundLambda::kPasses_DiffDaughters_Physical) + 1, "Passes_DiffDaughters_Physical");
        // post-seed cuts
        x_axis->SetBinLabel(static_cast<int>(EPreFoundLambda::kPasses_Max_DCAbtwDaughters) + 1, "Passes_Max_DCAbtwDaughters");
        // post-fit cuts
        x_axis->SetBinLabel(static_cast<int>(EPreFoundLambda::kPasses_AbsMax_Pz) + 1, "Passes_AbsMax_Pz");
        x_axis->SetBinLabel(static_cast<int>(EPreFoundLambda::kPasses_Max_Pt) + 1, "Passes_Max_Pt");
        x_axis->SetBinLabel(static_cast<int>(EPreFoundLambda::kPasses_Min_Pt) + 1, "Passes_Min_Pt");
        x_axis->SetBinLabel(static_cast<int>(EPreFoundLambda::kPasses_AbsMax_Rapidity) + 1, "Passes_AbsMax_Rapidity");
        x_axis->SetBinLabel(static_cast<int>(EPreFoundLambda::kPasses_Min_CPAwrtPV) + 1, "Passes_Min_CPAwrtPV");
        x_axis->SetBinLabel(static_cast<int>(EPreFoundLambda::kPasses_AbsMax_ArmRadiusDev) + 1, "Passes_AbsMax_ArmRadiusDev");
        x_axis->SetBinLabel(static_cast<int>(EPreFoundLambda::kPasses_Max_DCAwrtPV) + 1, "Passes_Max_DCAwrtPV");
        x_axis->SetBinLabel(static_cast<int>(EPreFoundLambda::kPasses_Max_Chi2NDF) + 1, "Passes_Max_Chi2NDF");
        // -- depend on (anti)protons
        x_axis->SetBinLabel(static_cast<int>(EPreFoundLambda::kPasses_AbsMax_Pz_Proton) + 1, "Passes_AbsMax_Pz_Proton");
        x_axis->SetBinLabel(static_cast<int>(EPreFoundLambda::kPasses_Max_Pt_Proton) + 1, "Passes_Max_Pt_Proton");
        x_axis->SetBinLabel(static_cast<int>(EPreFoundLambda::kPasses_Min_Pt_Proton) + 1, "Passes_Min_Pt_Proton");
        // -- depend on pi(minus/plus)
        x_axis->SetBinLabel(static_cast<int>(EPreFoundLambda::kPasses_AbsMax_Pz_Pion) + 1, "Passes_AbsMax_Pz_Pion");
        x_axis->SetBinLabel(static_cast<int>(EPreFoundLambda::kPasses_Max_Pt_Pion) + 1, "Passes_Max_Pt_Pion");
        x_axis->SetBinLabel(static_cast<int>(EPreFoundLambda::kPasses_Min_Pt_Pion) + 1, "Passes_Min_Pt_Pion");
    }

    // -- for (anti)h-dibaryons + mixed lambda pairs
    fHist_CutFlow_AntiHdibaryon =
        std::make_unique<TH1D>(std::format("CutFlow_{}", DB::Particles::Particle("AntiHdibaryon").acronym).c_str(), hist_title,  //
                               static_cast<int>(ELambdaPair::kNLambdaPairCuts), 0., static_cast<double>(ELambdaPair::kNLambdaPairCuts));
    fHist_CutFlow_Hdibaryon =
        std::make_unique<TH1D>(std::format("CutFlow_{}", DB::Particles::Particle("Hdibaryon").acronym).c_str(), hist_title,  //
                               static_cast<int>(ELambdaPair::kNLambdaPairCuts), 0., static_cast<double>(ELambdaPair::kNLambdaPairCuts));
    fHist_CutFlow_MixedLambdaPair =
        std::make_unique<TH1D>("CutFlow_MixedLL", hist_title,  //
                               static_cast<int>(ELambdaPair::kNLambdaPairCuts), 0., static_cast<double>(ELambdaPair::kNLambdaPairCuts));

    // -- define bin labels
    for (auto* hist_hdib : {fHist_CutFlow_AntiHdibaryon.get(), fHist_CutFlow_Hdibaryon.get(), fHist_CutFlow_MixedLambdaPair.get()}) {
        auto* x_axis = hist_hdib->GetXaxis();
        x_axis->SetBinLabel(static_cast<int>(ELambdaPair::kAllCombinations) + 1, "AllCombinations");
        // pre-seed cuts
        x_axis->SetBinLabel(static_cast<int>(ELambdaPair::kPasses_DiffLambdas_Logical) + 1, "Passes_DiffLambdas_Logical");
        x_axis->SetBinLabel(static_cast<int>(ELambdaPair::kPasses_DiffTracks_Logical) + 1, "Passes_DiffTracks_Logical");
        x_axis->SetBinLabel(static_cast<int>(ELambdaPair::kPasses_DiffTracks_Physical) + 1, "Passes_DiffTracks_Physical");
        x_axis->SetBinLabel(static_cast<int>(ELambdaPair::kPasses_DiffLambdas_Physical) + 1, "Passes_DiffLambdas_Physical");
        // post-seed cuts
        x_axis->SetBinLabel(static_cast<int>(ELambdaPair::kPasses_Max_DCAbtwDau) + 1, "Passes_Max_DCAbtwDau");
        // post-fit cuts
        x_axis->SetBinLabel(static_cast<int>(ELambdaPair::kPasses_AbsMax_Pz) + 1, "Passes_AbsMax_Pz");
        x_axis->SetBinLabel(static_cast<int>(ELambdaPair::kPasses_Max_Pt) + 1, "Passes_Max_Pt");
        x_axis->SetBinLabel(static_cast<int>(ELambdaPair::kPasses_Min_Pt) + 1, "Passes_Min_Pt");
        x_axis->SetBinLabel(static_cast<int>(ELambdaPair::kPasses_Min_Mass) + 1, "Passes_Min_Mass");
        x_axis->SetBinLabel(static_cast<int>(ELambdaPair::kPasses_Max_Mass) + 1, "Passes_Max_Mass");
        x_axis->SetBinLabel(static_cast<int>(ELambdaPair::kPasses_AbsMax_Rapidity) + 1, "Passes_AbsMax_Rapidity");
        x_axis->SetBinLabel(static_cast<int>(ELambdaPair::kPasses_Max_DecayLength) + 1, "Passes_Max_DecayLength");
        x_axis->SetBinLabel(static_cast<int>(ELambdaPair::kPasses_Min_CPAwrtPV) + 1, "Passes_Min_CPAwrtPV");
        x_axis->SetBinLabel(static_cast<int>(ELambdaPair::kPasses_Max_Chi2NDF) + 1, "Passes_Max_Chi2NDF");
        x_axis->SetBinLabel(static_cast<int>(ELambdaPair::kPasses_Max_Chi2CV) + 1, "Passes_Max_Chi2CV");
        // (anti)lambdas : depend on (anti)h-dibaryon decay vertex
        x_axis->SetBinLabel(static_cast<int>(ELambdaPair::kPasses_Max_L1_DecayLength) + 1, "Passes_Max_L1_DecayLength");
        x_axis->SetBinLabel(static_cast<int>(ELambdaPair::kPasses_Min_L1_DecayLength) + 1, "Passes_Min_L1_DecayLength");
        x_axis->SetBinLabel(static_cast<int>(ELambdaPair::kPasses_Min_L1_CPAwrtDV) + 1, "Passes_Min_L1_CPAwrtDV");
        x_axis->SetBinLabel(static_cast<int>(ELambdaPair::kPasses_Max_L2_DecayLength) + 1, "Passes_Max_L2_DecayLength");
        x_axis->SetBinLabel(static_cast<int>(ELambdaPair::kPasses_Min_L2_DecayLength) + 1, "Passes_Min_L2_DecayLength");
        x_axis->SetBinLabel(static_cast<int>(ELambdaPair::kPasses_Min_L2_CPAwrtDV) + 1, "Passes_Min_L2_CPAwrtDV");
    }
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
    fPrimaryVertex.SetCoordinates(static_cast<double>(fOutput.Event.PV_X), static_cast<double>(fOutput.Event.PV_Y),
                                  static_cast<double>(fOutput.Event.PV_Z));
    fPrimaryVertexKF = KF::Vertex::FromEvent(fOutput.Event);
    // cache bz
    fMagneticField = static_cast<double>(fOutput.Event.MagneticField);
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
    if (mc.FirstDau_McEntry > Common::DummyInt) {
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
    if (mc.LastDau_McEntry > Common::DummyInt) {
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

// Extract track information from a `POD::PreFoundLambda` as a `POD::Track`.
POD::Track Verifier::ExtractTrack(const POD::PreFoundLambda& pod_lambda, short charge) const {
    return {
        charge < 0 ? pod_lambda.Neg_EsdEntry : pod_lambda.Pos_EsdEntry,
        charge < 0 ? pod_lambda.Neg_State[0] : pod_lambda.Pos_State[0],
        charge < 0 ? pod_lambda.Neg_State[1] : pod_lambda.Pos_State[1],
        charge < 0 ? pod_lambda.Neg_State[2] : pod_lambda.Pos_State[2],
        charge < 0 ? pod_lambda.Neg_State[3] : pod_lambda.Pos_State[3],
        charge < 0 ? pod_lambda.Neg_State[4] : pod_lambda.Pos_State[4],
        charge < 0 ? pod_lambda.Neg_State[5] : pod_lambda.Pos_State[5],
        charge < 0 ? pod_lambda.Neg_PreDCAxy : pod_lambda.Pos_PreDCAxy,
        charge < 0 ? pod_lambda.Neg_PreDCAz : pod_lambda.Pos_PreDCAz,
        charge < 0 ? pod_lambda.Neg_NSigmasPion : pod_lambda.Pos_NSigmasPion,
        charge < 0 ? pod_lambda.Neg_NSigmasKaon : pod_lambda.Pos_NSigmasKaon,
        charge < 0 ? pod_lambda.Neg_NSigmasProton : pod_lambda.Pos_NSigmasProton,
        charge < 0 ? pod_lambda.Neg_CovMatrix : pod_lambda.Pos_CovMatrix,
        Common::DummyFloat,
        Common::DummyInt,
        charge < 0 ? pod_lambda.Neg_TPC_NCrossedRows : pod_lambda.Pos_TPC_NCrossedRows,
        charge < 0 ? pod_lambda.Neg_TPC_NClusters : pod_lambda.Pos_TPC_NClusters,
        charge < 0 ? pod_lambda.Neg_TPC_NClustersFindable : pod_lambda.Pos_TPC_NClustersFindable,
        charge < 0 ? pod_lambda.Neg_TPC_Chi2 : pod_lambda.Pos_TPC_Chi2,
        charge,
    };
}

// ## Single Pre-Found Lambda ZONE ## //

void Verifier::ProcessPreFoundLambda() {

    constexpr FitSetup setup = GetFitSetup();
    const KF::FitPolicy fit_policy{
        .pin_daughters = setup.pin_lambda_daughters,
        .daughters_already_pinned = false,
        .mother_mass = setup.pin_lambda_mass ? std::make_optional(DB::Particles::Particle("Lambda").mass) : std::nullopt,
        .prod_vertex = std::nullopt,
    };

    // loop over pre-found (anti)lambdas //
    for (std::size_t entry_lambda = 0; entry_lambda < fInput.PreFoundLambda.size(); ++entry_lambda) {

        const auto& in_lambda = fInput.PreFoundLambda[entry_lambda];

        // apply cuts (1) //
        if (!PreSeedCuts_Lambda(in_lambda)) continue;

        auto track_neg = ExtractTrack(in_lambda, -1);
        auto track_pos = ExtractTrack(in_lambda, +1);

        // PCAs //
        Seeder::HelixHelix::Cache pca_cache;
        Seeder::SeedsPair seeds_xy;
        auto [seed_neg, seed_pos] = Seeder::HelixHelix::FastPCAs(track_neg, track_pos, fMagneticField, pca_cache, seeds_xy);

        // apply cuts (2) //
        if (!PostSeedCuts_Lambda(seed_neg.pca, seed_pos.pca)) continue;

        // PCAs derivatives //
        auto [deriv_neg, deriv_pos] = Seeder::HelixHelix::ComputeDerivatives(seeds_xy, pca_cache);

        // loop over lambda / anti-lambda hypotheses //
        for (auto anti_channel : {false, true}) {

            const auto& decay_tree = HD::GetDecayTree(anti_channel);

            // fit vertex //
            auto fit = KF::FitVertex(track_neg, track_pos, decay_tree.neg, decay_tree.pos, {seed_neg, deriv_neg}, {seed_pos, deriv_pos},
                                     fMagneticField, fit_policy);

            // create storage+computation (anti)lambda //
            POD::Extended::PreFoundLambda new_lambda = CreateExtendedPreFoundLambda(in_lambda, fit, seed_neg.pca, seed_pos.pca, anti_channel);
            Cached::PreFoundLambda c_lambda(new_lambda, fPrimaryVertex);

            // apply more cuts (3) //
            if (!PostFitCuts_Lambda(c_lambda)) continue;

            // store reconstructed //
            if (anti_channel) {
                fTemp_AntiLambda.emplace_back(new_lambda);
            } else {
                fTemp_Lambda.emplace_back(new_lambda);
            }

            // store mc //
            if (!fSettings.IsMC) continue;

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
    }  // end of loop over pre-found (anti)lambdas
}

bool Verifier::PreSeedCuts_Lambda(const POD::PreFoundLambda& lambda) {
    FillHist(fHist_CutFlow_AntiLambda.get(), EPreFoundLambda::kAllPreFoundLambdas);
    FillHist(fHist_CutFlow_Lambda.get(), EPreFoundLambda::kAllPreFoundLambdas);

    if (HD::SameDaughterEntries(lambda)) return false;
    FillHist(fHist_CutFlow_AntiLambda.get(), EPreFoundLambda::kPasses_DiffDaughters_Logical);
    FillHist(fHist_CutFlow_Lambda.get(), EPreFoundLambda::kPasses_DiffDaughters_Logical);

    if (CMath::IsSameHelix(lambda.Neg_State, lambda.Pos_State, Cuts::PreFoundLambda::Max_TracksDeltaR, Cuts::PreFoundLambda::Max_TracksRelDeltaP)) {
        return false;
    }
    FillHist(fHist_CutFlow_AntiLambda.get(), EPreFoundLambda::kPasses_DiffDaughters_Physical);
    FillHist(fHist_CutFlow_Lambda.get(), EPreFoundLambda::kPasses_DiffDaughters_Physical);

    return true;
}

bool Verifier::PostSeedCuts_Lambda(const Seeder::PCA& pca_neg, const Seeder::PCA& pca_pos) {

    // if (CMath::Distance(pca_neg.xyz, pca_pos.xyz) > Cuts::PreFoundLambda::Max_DCAbtwDaughters) return false; // PENDING: temporarily turned off
    // FillHist(fHist_CutFlow_AntiLambda.get(), EPreFoundLambda::kPasses_Max_DCAbtwDaughters); // PENDING: temporarily turned off
    // FillHist(fHist_CutFlow_Lambda.get(), EPreFoundLambda::kPasses_Max_DCAbtwDaughters); // PENDING: temporarily turned off

    return true;
}

bool Verifier::PostFitCuts_Lambda(const Cached::PreFoundLambda& c_lambda) {
    auto* hist_cut_flow = c_lambda.IsAntiLambda ? fHist_CutFlow_AntiLambda.get() : fHist_CutFlow_Lambda.get();

    if (std::abs(static_cast<double>(c_lambda.Pz)) > Cuts::PreFoundLambda::AbsMax_Pz) return false;
    FillHist(hist_cut_flow, EPreFoundLambda::kPasses_AbsMax_Pz);

    if (c_lambda.Pt() > Cuts::PreFoundLambda::Max_Pt) return false;
    FillHist(hist_cut_flow, EPreFoundLambda::kPasses_Max_Pt);

    // if (c_lambda.Pt() < Cuts::PreFoundLambda::Min_Pt) return false; // PENDING: temporarily turned off
    // FillHist(hist_cut_flow, EPreFoundLambda::kPasses_Min_Pt); // PENDING: temporarily turned off

    if (std::abs(c_lambda.Rapidity()) > Cuts::PreFoundLambda::AbsMax_Rapidity) return false;
    FillHist(hist_cut_flow, EPreFoundLambda::kPasses_AbsMax_Rapidity);

    // if (c_lambda.Mass() < Cuts::PreFoundLambda::Min_Mass) return false; // PENDING: temporarily turned off
    // FillHist(hist_cut_flow, EPreFoundLambda::kPasses_Min_Mass); // PENDING: temporarily turned off

    // if (c_lambda.Mass() > Cuts::PreFoundLambda::Max_Mass) return false; // PENDING: temporarily turned off
    // FillHist(hist_cut_flow, EPreFoundLambda::kPasses_Max_Mass); // PENDING: temporarily turned off

    if (c_lambda.CPA_wrt_PV() < Cuts::PreFoundLambda::Min_CPAwrtPV) return false;
    FillHist(hist_cut_flow, EPreFoundLambda::kPasses_Min_CPAwrtPV);

    if (std::abs(c_lambda.ArmRadiusDev()) > Cuts::PreFoundLambda::AbsMax_ArmRadiusDev) return false;
    FillHist(hist_cut_flow, EPreFoundLambda::kPasses_AbsMax_ArmRadiusDev);

    // if (c_lambda.DCA_wrt_PV() > Cuts::PreFoundLambda::Max_DCAwrtPV) return false; // PENDING: temporarily turned off
    // FillHist(hist_cut_flow, EPreFoundLambda::kPasses_Max_DCAwrtPV); // PENDING: temporarily turned off

    if (static_cast<double>(c_lambda.Chi2NDF) > 2.5) return false;
    FillHist(hist_cut_flow, EPreFoundLambda::kPasses_Max_Chi2NDF);

    // depend on (anti)protons //

    if (std::abs(c_lambda.Pr_Pz()) > Cuts::PreFoundLambda::AbsMax_Pz_Proton) return false;
    FillHist(hist_cut_flow, EPreFoundLambda::kPasses_AbsMax_Pz_Proton);

    if (c_lambda.Pr_Pt() > Cuts::PreFoundLambda::Max_Pt_Proton) return false;
    FillHist(hist_cut_flow, EPreFoundLambda::kPasses_Max_Pt_Proton);

    // if (c_lambda.Pr_Pt() < Cuts::PreFoundLambda::Min_Pt_Proton) return false; // PENDING: temporarily turned off
    // FillHist(hist_cut_flow, EPreFoundLambda::kPasses_Min_Pt_Proton); // PENDING: temporarily turned off

    // depend on pi(minus/plus) //

    if (std::abs(c_lambda.Pi_Pz()) > Cuts::PreFoundLambda::AbsMax_Pz_Pion) return false;
    FillHist(hist_cut_flow, EPreFoundLambda::kPasses_AbsMax_Pz_Pion);

    if (c_lambda.Pi_Pt() > Cuts::PreFoundLambda::Max_Pt_Pion) return false;
    FillHist(hist_cut_flow, EPreFoundLambda::kPasses_Max_Pt_Pion);

    // if (c_lambda.Pi_Pt() < Cuts::PreFoundLambda::Min_Pt_Pion) return false; // PENDING: temporarily turned off
    // FillHist(hist_cut_flow, EPreFoundLambda::kPasses_Min_Pt_Pion); // PENDING: temporarily turned off

    return true;
}

POD::Extended::McParticle Verifier::BuildMcPreFoundLambda(const POD::Extended::McParticle& mc_neg, const POD::Extended::McParticle& mc_pos,
                                                          const HD::DecayTree& decay_pid) {
    POD::Extended::McParticle mc_lambda;

    // -- fill hybridness, independently of no common mother
    mc_lambda.IsHybrid = mc_neg.IsTrueSignal != mc_pos.IsTrueSignal || mc_neg.SignalID != mc_pos.SignalID;

    auto mc_entry = MC::FindCommonMotherMcEntry(mc_neg, mc_pos);
    if (!mc_entry.has_value()) return mc_lambda;

    // fill rest of values //
    static_cast<POD::McParticle&>(mc_lambda) = fInput.McParticle[mc_entry.value()];
    MC::Apply(mc_lambda, MC::HdibaryonRules::ClassifyDownstream(mc_lambda, fInput.McParticle, decay_pid, decay_pid.lambda.pdg_code, true));

    return mc_lambda;
}

POD::Extended::PreFoundLambda Verifier::CreateExtendedPreFoundLambda(const POD::PreFoundLambda& old_lambda, const KF::FitResult& fit,
                                                                     const Seeder::PCA& pca_neg, const Seeder::PCA& pca_pos, bool anti_lambda) {
    // profit from initialization list to extend data
    POD::Extended::PreFoundLambda new_lambda{old_lambda,
                                             static_cast<float>(fit.mother.Px()),
                                             static_cast<float>(fit.mother.Py()),
                                             static_cast<float>(fit.mother.Pz()),
                                             static_cast<float>(fit.mother.E()),
                                             fit.mother.Cov<7>(),
                                             static_cast<float>(fit.mother.Chi2NDF()),
                                             static_cast<float>(pca_neg.X()),
                                             static_cast<float>(pca_neg.Y()),
                                             static_cast<float>(pca_neg.Z()),
                                             static_cast<float>(fit.Dau1_Px()),
                                             static_cast<float>(fit.Dau1_Py()),
                                             static_cast<float>(fit.Dau1_Pz()),
                                             static_cast<float>(fit.Dau1_E()),
                                             static_cast<float>(pca_pos.X()),
                                             static_cast<float>(pca_pos.Y()),
                                             static_cast<float>(pca_pos.Z()),
                                             static_cast<float>(fit.Dau2_Px()),
                                             static_cast<float>(fit.Dau2_Py()),
                                             static_cast<float>(fit.Dau2_Pz()),
                                             static_cast<float>(fit.Dau2_E()),
                                             anti_lambda};
    // override previous info
    new_lambda.Decay_X = static_cast<float>(fit.mother.X());
    new_lambda.Decay_Y = static_cast<float>(fit.mother.Y());
    new_lambda.Decay_Z = static_cast<float>(fit.mother.Z());
    new_lambda.DcaV0Daughters = static_cast<float>(CMath::Distance(pca_neg.xyz, pca_pos.xyz));
    new_lambda.Neg_PCAwrtV0_Px = static_cast<float>(pca_neg.Px());
    new_lambda.Neg_PCAwrtV0_Py = static_cast<float>(pca_neg.Py());
    new_lambda.Neg_PCAwrtV0_Pz = static_cast<float>(pca_neg.Pz());
    new_lambda.Pos_PCAwrtV0_Px = static_cast<float>(pca_pos.Px());
    new_lambda.Pos_PCAwrtV0_Py = static_cast<float>(pca_pos.Py());
    new_lambda.Pos_PCAwrtV0_Pz = static_cast<float>(pca_pos.Pz());

    return new_lambda;
}

// ## H-dibaryon ZONE ## //

void Verifier::VerifyLambdaPair(bool anti_channel_l1, bool anti_channel_l2) {

    // determine rules based on the species assigned to each (anti)lambda //
    // -- the mixed channel is the lambda + anti-lambda background, which no (anti)h-dibaryon can decay into
    const bool mixed_channel = anti_channel_l1 != anti_channel_l2;
    const bool anti_channel = anti_channel_l1 && anti_channel_l2;

    const auto& decay_pid_l1 = HD::GetDecayTree(anti_channel_l1);
    const auto& decay_pid_l2 = HD::GetDecayTree(anti_channel_l2);
    // NOTE: "Unknown" carries a pdg code no true particle can match, so the mixed channel is never labelled as signal
    const int pdg_code_hypothesis = mixed_channel ? DB::Particles::Particle("Unknown").pdg_code : decay_pid_l1.hdibaryon.pdg_code;

    auto* hist_cut_flow = mixed_channel  ? fHist_CutFlow_MixedLambdaPair.get()
                          : anti_channel ? fHist_CutFlow_AntiHdibaryon.get()
                                         : fHist_CutFlow_Hdibaryon.get();

    const auto& input_lambdas_l1 = anti_channel_l1 ? fTemp_AntiLambda : fTemp_Lambda;
    const auto& input_mc_lambdas_l1 = anti_channel_l1 ? fTemp_MC_AntiLambda : fTemp_MC_Lambda;
    const auto& input_mc_lambdas_l1_neg = anti_channel_l1 ? fTemp_MC_AntiLambda_Neg : fTemp_MC_Lambda_Neg;
    const auto& input_mc_lambdas_l1_pos = anti_channel_l1 ? fTemp_MC_AntiLambda_Pos : fTemp_MC_Lambda_Pos;
    const auto& input_lambdas_l2 = anti_channel_l2 ? fTemp_AntiLambda : fTemp_Lambda;
    const auto& input_mc_lambdas_l2 = anti_channel_l2 ? fTemp_MC_AntiLambda : fTemp_MC_Lambda;
    const auto& input_mc_lambdas_l2_neg = anti_channel_l2 ? fTemp_MC_AntiLambda_Neg : fTemp_MC_Lambda_Neg;
    const auto& input_mc_lambdas_l2_pos = anti_channel_l2 ? fTemp_MC_AntiLambda_Pos : fTemp_MC_Lambda_Pos;

    // -- the (anti)h-dibaryon mass is never pinned, because it's the property to study
    constexpr FitSetup setup = GetFitSetup();
    const KF::FitPolicy fit_policy{
        .pin_daughters = setup.pin_lambdas,
        .daughters_already_pinned = setup.lambdas_on_shell,
        .mother_mass = std::nullopt,
        .prod_vertex = setup.pin_hdib_to_pv ? std::make_optional(fPrimaryVertexKF) : std::nullopt,
    };

    // loop over all possible pairs of (anti)lambdas //
    for (std::size_t entry_lambda1 = 0; entry_lambda1 < input_lambdas_l1.size(); ++entry_lambda1) {
        const auto& lambda1 = input_lambdas_l1[entry_lambda1];  // cache index lookup

        // -- in the signal channels both legs come from a single collection and the pair is unordered, so only the
        //    upper triangle is scanned; in the mixed one the two collections are disjoint, so it's the full product
        for (std::size_t entry_lambda2 = mixed_channel ? 0 : entry_lambda1 + 1; entry_lambda2 < input_lambdas_l2.size(); ++entry_lambda2) {
            const auto& lambda2 = input_lambdas_l2[entry_lambda2];  // cache index lookup

            // apply cuts (1) //
            if (!PreSeedCuts_Hdibaryon(lambda1, lambda2, hist_cut_flow)) continue;

            // PCAs //
            Seeder::LineLine::Cache pca_cache;
            auto [seed_lambda1, seed_lambda2] = Seeder::LineLine::FastPCAs(lambda1, lambda2, pca_cache);

            // apply cuts (2) //
            if (!PostSeedCuts_Hdibaryon(seed_lambda1.pca, seed_lambda2.pca, hist_cut_flow)) continue;

            // PCAs derivatives //
            auto [deriv_lambda1, deriv_lambda2] = Seeder::LineLine::ComputeDerivatives(pca_cache);

            // fit vertex //
            auto fit = KF::FitVertex(lambda1, lambda2, decay_pid_l1.lambda, decay_pid_l2.lambda, {seed_lambda1, deriv_lambda1},
                                     {seed_lambda2, deriv_lambda2}, fit_policy);

            // create storage+computation units //
            POD::LambdaPair hdib = CreateLambdaPair(fit, seed_lambda1.pca, seed_lambda2.pca);
            Cached::Hdibaryon c_hdib(hdib, lambda1, lambda2, fPrimaryVertex);

            // apply cuts (3) //
            if (!PostFitCuts_Hdibaryon(c_hdib, hist_cut_flow)) continue;

            // store reconstructed //
            fOutput.Hdibaryon.emplace_back(hdib);
            fOutput.Lambda1.emplace_back(lambda1);
            fOutput.Lambda2.emplace_back(lambda2);

            // store mc //
            if (fSettings.IsMC) {
                fOutput.MC_Hdibaryon.emplace_back(
                    BuildMcHdibaryon(input_mc_lambdas_l1[entry_lambda1], input_mc_lambdas_l2[entry_lambda2], pdg_code_hypothesis));
                fOutput.MC_Lambda1.emplace_back(input_mc_lambdas_l1[entry_lambda1]);
                fOutput.MC_Lambda1_Neg.emplace_back(input_mc_lambdas_l1_neg[entry_lambda1]);
                fOutput.MC_Lambda1_Pos.emplace_back(input_mc_lambdas_l1_pos[entry_lambda1]);
                fOutput.MC_Lambda2.emplace_back(input_mc_lambdas_l2[entry_lambda2]);
                fOutput.MC_Lambda2_Neg.emplace_back(input_mc_lambdas_l2_neg[entry_lambda2]);
                fOutput.MC_Lambda2_Pos.emplace_back(input_mc_lambdas_l2_pos[entry_lambda2]);
            }
        }  // end of loop over lambda2
    }  // end of loop over lambda1
}

bool Verifier::PreSeedCuts_Hdibaryon(const POD::Extended::PreFoundLambda& lambda1, const POD::Extended::PreFoundLambda& lambda2,
                                     TH1D* hist_cut_flow) {
    FillHist(hist_cut_flow, ELambdaPair::kAllCombinations);

    if (HD::SameLambdasEntries(lambda1, lambda2)) return false;  // order is important; apply this before `SameDaughterEntries()`
    FillHist(hist_cut_flow, ELambdaPair::kPasses_DiffLambdas_Logical);

    if (HD::SameDaughterEntries(lambda1, lambda2)) return false;
    FillHist(hist_cut_flow, ELambdaPair::kPasses_DiffTracks_Logical);

    if (CMath::IsSameHelix(lambda1.Neg_State, lambda2.Neg_State, Cuts::PreFoundLambda::Max_TracksDeltaR, Cuts::PreFoundLambda::Max_TracksRelDeltaP) ||
        CMath::IsSameHelix(lambda1.Pos_State, lambda2.Pos_State, Cuts::PreFoundLambda::Max_TracksDeltaR, Cuts::PreFoundLambda::Max_TracksRelDeltaP)) {
        return false;
    }
    FillHist(hist_cut_flow, ELambdaPair::kPasses_DiffTracks_Physical);

    if (CMath::SquaredDistance(lambda1.Decay_X, lambda1.Decay_Y, lambda1.Decay_Z, lambda2.Decay_X, lambda2.Decay_Y, lambda2.Decay_Z) <
        Cuts::LambdaPair::Min_DistBtwLambdaDVs * Cuts::LambdaPair::Min_DistBtwLambdaDVs) {
        return false;
    }
    FillHist(hist_cut_flow, ELambdaPair::kPasses_DiffLambdas_Physical);

    return true;
}

bool Verifier::PostSeedCuts_Hdibaryon(const Seeder::PCA& pca_lambda1, const Seeder::PCA& pca_lambda2, TH1D* hist_cut_flow) {

    // if (CMath::Distance(pca_lambda1.xyz, pca_lambda2.xyz) > Cuts::LambdaPair::Max_DCAbtwDau) return false; // PENDING: temporarily turned off
    // FillHist(hist_cut_flow, ELambdaPair::kPasses_Max_DCAbtwDau); // PENDING: temporarily turned off

    return true;
}

bool Verifier::PostFitCuts_Hdibaryon(const Cached::Hdibaryon& c_hdib, TH1D* hist_cut_flow) {

    if (std::abs(static_cast<double>(c_hdib.Pz)) > Cuts::LambdaPair::AbsMax_Pz) return false;
    FillHist(hist_cut_flow, ELambdaPair::kPasses_AbsMax_Pz);

    if (c_hdib.Pt() > Cuts::LambdaPair::Max_Pt) return false;
    FillHist(hist_cut_flow, ELambdaPair::kPasses_Max_Pt);

    if (c_hdib.Pt() < Cuts::LambdaPair::Min_Pt) return false;
    FillHist(hist_cut_flow, ELambdaPair::kPasses_Min_Pt);

    if (std::abs(c_hdib.Rapidity()) > Cuts::LambdaPair::AbsMax_Rapidity) return false;
    FillHist(hist_cut_flow, ELambdaPair::kPasses_AbsMax_Rapidity);

    // if (static_cast<double>(c_hdib.DecayLength) > Cuts::LambdaPair::Max_DecayLength) return false; // PENDING: temporarily turned off, need to
    // re-tune FillHist(hist_cut_flow, ELambdaPair::kPasses_Max_DecayLength); // PENDING: temporarily turned off, need to re-tune

    // if (c_hdib.CPA_wrt_PV() < Cuts::LambdaPair::Min_CPAwrtPV) return false; // PENDING: temporarily turned off
    // FillHist(hist_cut_flow, ELambdaPair::kPasses_Min_CPAwrtPV); // PENDING: temporarily turned off

    if (static_cast<double>(c_hdib.Chi2NDF) > Cuts::LambdaPair::Max_Chi2NDF) return false;
    FillHist(hist_cut_flow, ELambdaPair::kPasses_Max_Chi2NDF);

    if (static_cast<double>(c_hdib.Chi2CV) > 5.) return false;
    FillHist(hist_cut_flow, ELambdaPair::kPasses_Max_Chi2CV);

    // (anti)lambda : depend on (anti)h-dibaryon decay vertex //

    if (c_hdib.Lambda1_DecayLength() > Cuts::PreFoundLambda::Max_DecayLength) return false;
    FillHist(hist_cut_flow, ELambdaPair::kPasses_Max_L1_DecayLength);

    // if (c_hdib.Lambda1_DecayLength() < Cuts::PreFoundLambda::Min_DecayLength) return false; // PENDING: temporarily turned off
    // FillHist(hist_cut_flow, ELambdaPair::kPasses_Min_L1_DecayLength); // PENDING: temporarily turned off

    // if (c_hdib.Lambda1_CPA_wrt_DV() < Cuts::PreFoundLambda::Min_CPAwrtDV) return false; // PENDING: temporarily turned off
    // FillHist(hist_cut_flow, ELambdaPair::kPasses_Min_L1_CPAwrtDV); // PENDING: temporarily turned off

    if (c_hdib.Lambda2_DecayLength() > Cuts::PreFoundLambda::Max_DecayLength) return false;
    FillHist(hist_cut_flow, ELambdaPair::kPasses_Max_L2_DecayLength);

    // if (c_hdib.Lambda2_DecayLength() < Cuts::PreFoundLambda::Min_DecayLength) return false; // PENDING: temporarily turned off
    // FillHist(hist_cut_flow, ELambdaPair::kPasses_Min_L2_DecayLength); // PENDING: temporarily turned off

    // if (c_hdib.Lambda2_CPA_wrt_DV() < Cuts::PreFoundLambda::Min_CPAwrtDV) return false; // PENDING: temporarily turned off
    // FillHist(hist_cut_flow, ELambdaPair::kPasses_Min_L2_CPAwrtDV); // PENDING: temporarily turned off

    return true;
}

POD::Extended::McParticle Verifier::BuildMcHdibaryon(const POD::Extended::McParticle& mc_lambda1, const POD::Extended::McParticle& mc_lambda2,
                                                     int pdg_code_hypothesis) {
    POD::Extended::McParticle mc_hdib;

    // -- fill hybridness, independently of no common mother
    mc_hdib.IsHybrid = mc_lambda1.IsHybrid || mc_lambda2.IsHybrid || mc_lambda1.IsTrueSignal != mc_lambda2.IsTrueSignal ||
                       mc_lambda1.SignalID != mc_lambda2.SignalID;

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
    mc_hdib.IsTrue = mc_hdib.PdgCode == pdg_code_hypothesis;
    mc_hdib.IsGen1Signal = false;
    mc_hdib.IsGen2Signal = false;
    mc_hdib.IsTrueSignal = mc_hdib.IsTrue;
    mc_hdib.IsSecondary = mc_hdib.IsSecFromMat || mc_hdib.IsSecFromWeak;

    return mc_hdib;
}

POD::LambdaPair Verifier::CreateLambdaPair(const KF::FitResult& fit, const Seeder::PCA& pca_lambda1, const Seeder::PCA& pca_lambda2) {
    POD::LambdaPair new_hdib;
    // candidate info
    new_hdib.Decay_X = static_cast<float>(fit.mother.X());
    new_hdib.Decay_Y = static_cast<float>(fit.mother.Y());
    new_hdib.Decay_Z = static_cast<float>(fit.mother.Z());
    new_hdib.Px = static_cast<float>(fit.mother.Px());
    new_hdib.Py = static_cast<float>(fit.mother.Py());
    new_hdib.Pz = static_cast<float>(fit.mother.Pz());
    new_hdib.Energy = static_cast<float>(fit.mother.E());
    new_hdib.Chi2NDF = static_cast<float>(fit.mother.Chi2NDF());
    // available when fit constrained to a production vertex
    new_hdib.CV_X = static_cast<float>(fit.at_pv ? fit.at_pv->X() : Common::DummyDouble);
    new_hdib.CV_Y = static_cast<float>(fit.at_pv ? fit.at_pv->Y() : Common::DummyDouble);
    new_hdib.CV_Z = static_cast<float>(fit.at_pv ? fit.at_pv->Z() : Common::DummyDouble);
    new_hdib.CV_Px = static_cast<float>(fit.at_pv ? fit.at_pv->Px() : Common::DummyDouble);
    new_hdib.CV_Py = static_cast<float>(fit.at_pv ? fit.at_pv->Py() : Common::DummyDouble);
    new_hdib.CV_Pz = static_cast<float>(fit.at_pv ? fit.at_pv->Pz() : Common::DummyDouble);
    new_hdib.CV_Energy = static_cast<float>(fit.at_pv ? fit.at_pv->E() : Common::DummyDouble);
    new_hdib.CV_DecayLength = static_cast<float>(fit.at_pv ? fit.at_pv->DecayLength() : Common::DummyDouble);
    new_hdib.CV_DecayLengthErr = static_cast<float>(fit.at_pv ? fit.at_pv->DecayLengthErr().value_or(Common::DummyDouble) : Common::DummyDouble);
    new_hdib.Chi2CV = static_cast<float>(fit.at_pv ? fit.at_pv->Chi2() - fit.mother.Chi2() : Common::DummyDouble);
    // (anti)lambda 1
    new_hdib.Lambda1_PCAwrtDV_X = static_cast<float>(pca_lambda1.X());
    new_hdib.Lambda1_PCAwrtDV_Y = static_cast<float>(pca_lambda1.Y());
    new_hdib.Lambda1_PCAwrtDV_Z = static_cast<float>(pca_lambda1.Z());
    new_hdib.Lambda1_Fit_Px = static_cast<float>(fit.Dau1_Px());
    new_hdib.Lambda1_Fit_Py = static_cast<float>(fit.Dau1_Py());
    new_hdib.Lambda1_Fit_Pz = static_cast<float>(fit.Dau1_Pz());
    new_hdib.Lambda1_Fit_Energy = static_cast<float>(fit.Dau1_E());
    // (anti)lambda 2
    new_hdib.Lambda2_PCAwrtDV_X = static_cast<float>(pca_lambda2.X());
    new_hdib.Lambda2_PCAwrtDV_Y = static_cast<float>(pca_lambda2.Y());
    new_hdib.Lambda2_PCAwrtDV_Z = static_cast<float>(pca_lambda2.Z());
    new_hdib.Lambda2_Fit_Px = static_cast<float>(fit.Dau2_Px());
    new_hdib.Lambda2_Fit_Py = static_cast<float>(fit.Dau2_Py());
    new_hdib.Lambda2_Fit_Pz = static_cast<float>(fit.Dau2_Pz());
    new_hdib.Lambda2_Fit_Energy = static_cast<float>(fit.Dau2_E());

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

bool Verifier::EndOfAnalysis() {

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
    fHist_CutFlow_MixedLambdaPair->Write();
    Logger::Info(__FUNCTION__, "- TH1D \"{}\"", fHist_CutFlow_MixedLambdaPair->GetName());

    Logger::Info(__FUNCTION__, "All done.");

    return fHist_EventCounter->GetEntries() != 0;
}

}  // namespace T2DS
