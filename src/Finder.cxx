#include "Finder/Finder.hxx"

#include <cstddef>
#include <filesystem>
#include <memory>
#include <string_view>

#include "App/DB_Particles.hxx"
#include "App/DB_ReactionChannels.hxx"
#include "App/Logger.hxx"
#include "Finder/FinderCuts.hxx"
#include "KalmanFitter/BaseKalmanFitter.hxx"
#include "KalmanFitter/KalmanFitterChannelA.hxx"
#include "KalmanFitter/KalmanFitterChannelD.hxx"
#include "Math/BaseMath.hxx"
#include "Math/Constants.hxx"
#include "MonteCarlo/MonteCarloChannelA.hxx"
#include "MonteCarlo/MonteCarloChannelD.hxx"
#include "Schema/SchemaFlatPrimitives.hxx"
#include "Schema/SchemaVectorMcTracks.hxx"
#include "Schema/SchemaVectorPrimitives.hxx"
#include "Seeder/SeederHelixLine.hxx"
#include "Seeder/SeederLineLine.hxx"
#include "View/ViewVectorInjected.hxx"
#include "View/ViewVectorMcTracks.hxx"
#include "View/ViewVectorMcV0s.hxx"
#include "View/ViewVectorTracks.hxx"

namespace Tree2Secondaries {

namespace KF = KalmanFitter;
namespace MC = MonteCarlo;

bool Finder::Initialize() {

    fInputChain_PackedEvents = std::make_unique<TChain>(Const::TreeName_PackedEvents.data());
    for (const auto& path : fSettings.PathInputFiles) {
        if (fInputChain_PackedEvents->Add(path.c_str()) == 0) {
            Logger::Error(__FUNCTION__, "Couldn't add TFile {}", path);
        }
    }
    if (!fInputChain_PackedEvents->GetEntries()) {
        Logger::Error(__FUNCTION__, "Couldn't manage to read any entry.");
        return false;
    }
    Logger::Info(__FUNCTION__, "TChain \"{}\" loaded successfully with {} trees and {} total entries.", fInputChain_PackedEvents->GetName(),
                 fInputChain_PackedEvents->GetNtrees(), fInputChain_PackedEvents->GetEntries());

    ReadInputBranches();

    if (!PrepareOutputFile()) return false;

    PrepareOutputHistograms();

    if (!PrepareOutputTree()) return false;
    CreateOutputBranches();

    if (fSettings.IsMC) {
        if (!Injected_PrepareOutputTree()) return false;
        fOutput_Injected.CreateBranches(fOutputTree_Injected.get());
    }

    Logger::Info(__FUNCTION__, "Finder initialized successfully.");

    return true;
}

// ## INPUT ZONE ## //

void Finder::ReadInputBranches() {

    // by default, turn off all branches
    fInputChain_PackedEvents->SetBranchStatus("*", false);
    // connect input branches to memory
    fInput_Event.ReadBranches(fInputChain_PackedEvents.get(), fSettings.IsMC);
    if (fSettings.IsMC) fInput_Injected.ReadBranches(fInputChain_PackedEvents.get(), true);

    // -- depending on reaction channels
    fInput_AntiLambdas.ReadBranches(fInputChain_PackedEvents.get(), Particles::Particle("AntiLambda").acronym);
    fInput_Lambdas.ReadBranches(fInputChain_PackedEvents.get(), Particles::Particle("Lambda").acronym);
    fInput_KaonsZeroShort.ReadBranches(fInputChain_PackedEvents.get(), Particles::Particle("KaonZeroShort").acronym);

    fInput_NegKaons.ReadBranches(fInputChain_PackedEvents.get(), false, Particles::Particle("NegKaon").acronym);
    fInput_PosKaons.ReadBranches(fInputChain_PackedEvents.get(), false, Particles::Particle("PosKaon").acronym);

    if (fSettings.IsMC) {
        fInput_MC_AntiLambdas.ReadBranches(fInputChain_PackedEvents.get(), Particles::Particle("AntiLambda").acronym);
        fInput_MC_AntiLambdas_Neg.ReadBranches(fInputChain_PackedEvents.get(), false,
                                               std::format("{}_Neg", Particles::Particle("AntiLambda").acronym));
        fInput_MC_AntiLambdas_Pos.ReadBranches(fInputChain_PackedEvents.get(), false,
                                               std::format("{}_Pos", Particles::Particle("AntiLambda").acronym));

        fInput_MC_Lambdas.ReadBranches(fInputChain_PackedEvents.get(), Particles::Particle("Lambda").acronym);
        fInput_MC_Lambdas_Neg.ReadBranches(fInputChain_PackedEvents.get(), false, std::format("{}_Neg", Particles::Particle("Lambda").acronym));
        fInput_MC_Lambdas_Pos.ReadBranches(fInputChain_PackedEvents.get(), false, std::format("{}_Pos", Particles::Particle("Lambda").acronym));

        fInput_MC_KaonsZeroShort.ReadBranches(fInputChain_PackedEvents.get(), Particles::Particle("KaonZeroShort").acronym);
        fInput_MC_KaonsZeroShort_Neg.ReadBranches(fInputChain_PackedEvents.get(), false,
                                                  std::format("{}_Neg", Particles::Particle("KaonZeroShort").acronym));
        fInput_MC_KaonsZeroShort_Pos.ReadBranches(fInputChain_PackedEvents.get(), false,
                                                  std::format("{}_Pos", Particles::Particle("KaonZeroShort").acronym));

        fInput_MC_NegKaons.ReadBranches(fInputChain_PackedEvents.get(), true, Particles::Particle("NegKaon").acronym);
        fInput_MC_PosKaons.ReadBranches(fInputChain_PackedEvents.get(), true, Particles::Particle("PosKaon").acronym);
    }
}

// ## OUTPUT ZONE ## //

bool Finder::PrepareOutputFile() {

    const std::filesystem::path output_path(fSettings.PathOutputFile);
    if (output_path.has_parent_path()) std::filesystem::create_directories(output_path.parent_path());

    fOutputFile = std::unique_ptr<TFile>(TFile::Open(fSettings.PathOutputFile.c_str(), "RECREATE"));
    if (!fOutputFile) {
        Logger::Error(__FUNCTION__, "Couldn't create TFile {}", fSettings.PathOutputFile);
        return false;
    }

    return true;
}

void Finder::PrepareOutputHistograms() {
    // event counter //
    fHist_EventCounter = std::make_unique<TH1D>("N_Events", ";;N_Events", 1, 0, 1);
    // cut flows //
    constexpr int x_nbins = 20;
    constexpr float x_min = 0.;
    constexpr float x_max = 20.;
    std::string_view hist_title = ";Cut N;N Passed Cut";
    fHist_CutFlow = std::make_unique<TH1D>("CutFlow", hist_title.data(), x_nbins, x_min, x_max);
    fHist_CutFlow_AntiChannel = std::make_unique<TH1D>("CutFlow_Anti", hist_title.data(), x_nbins, x_min, x_max);
}

bool Finder::PrepareOutputTree() {

    auto tree_name = std::format("FoundChannel{}", fSettings.ReactionChannel.name);

    fOutputTree = std::make_unique<TTree>(tree_name.data(), "");
    if (!fOutputTree) {
        Logger::Error(__FUNCTION__, "Couldn't create TTree \"{}\"", tree_name);
        return false;
    }

    return true;
}

void Finder::CreateOutputBranches() {
    switch (fSettings.ReactionChannel.name) {
        case 'A': {
            fOutput_ChannelA.CreateBranches(fOutputTree.get(), fSettings.IsMC);
            if (fSettings.IsMC) {
                fOutput_MC_ChannelA.CreateBranches(fOutputTree.get());
                fOutput_MC_ChannelA_V0A.CreateBranches(fOutputTree.get(), "V0A");
                fOutput_MC_ChannelA_V0A_Neg.CreateBranches(fOutputTree.get(), false, "V0A_Neg");
                fOutput_MC_ChannelA_V0A_Pos.CreateBranches(fOutputTree.get(), false, "V0A_Pos");
                fOutput_MC_ChannelA_V0B.CreateBranches(fOutputTree.get(), "V0B");
                fOutput_MC_ChannelA_V0B_Neg.CreateBranches(fOutputTree.get(), false, "V0B_Neg");
                fOutput_MC_ChannelA_V0B_Pos.CreateBranches(fOutputTree.get(), false, "V0B_Pos");
            }
            break;
        }
        case 'D': {
            fOutput_ChannelD.CreateBranches(fOutputTree.get(), fSettings.IsMC);
            if (fSettings.IsMC) {
                fOutput_MC_ChannelD.CreateBranches(fOutputTree.get());
                fOutput_MC_ChannelD_V0.CreateBranches(fOutputTree.get(), "V0");
                fOutput_MC_ChannelD_V0_Neg.CreateBranches(fOutputTree.get(), false, "V0_Neg");
                fOutput_MC_ChannelD_V0_Pos.CreateBranches(fOutputTree.get(), false, "V0_Pos");
                fOutput_MC_ChannelD_Kaon.CreateBranches(fOutputTree.get(), true, "Kaon");
            }
            break;
        }
        default:
            break;
    }
}

// ## Event ZONE ## //

void Finder::ProcessEvent() {  //
    fHist_EventCounter->Fill(0.);
}

// ## OUTPUT / Injected ZONE ## //

bool Finder::Injected_PrepareOutputTree() {

    fOutputTree_Injected = std::make_unique<TTree>(Const::TreeName_Injected.data(), "");
    if (fOutputTree_Injected == nullptr) {
        Logger::Error(__FUNCTION__, "Couldn't create TTree \"{}\"", Const::TreeName_Injected);
        return false;
    }

    return true;
}

// Flatten the injected anti-sexaquark vector tree and store it in TTree `Injected`.
void Finder::ProcessInjected() {
    for (std::size_t idx = 0; idx < NumberInjected(); ++idx) {
        View::VecInjected v(&fInput_Injected, idx);
        // assign values //
        fOutput_Injected.run_number = fInput_Event.run_number;
        fOutput_Injected.dir_number = fInput_Event.dir_number;
        fOutput_Injected.event_number = fInput_Event.event_number;
        fOutput_Injected.reaction_id = v.ReactionID();
        fOutput_Injected.sv.x = v.SV_X();
        fOutput_Injected.sv.y = v.SV_Y();
        fOutput_Injected.sv.z = v.SV_Z();
        fOutput_Injected.lv.px = v.Px();
        fOutput_Injected.lv.py = v.Py();
        fOutput_Injected.lv.pz = v.Pz();
        fOutput_Injected.lv.energy = v.Energy(fSettings.SexaquarkMass);
        fOutput_Injected.lv_nucleon.px = v.Nucleon_Px();
        fOutput_Injected.lv_nucleon.py = v.Nucleon_Py();
        fOutput_Injected.lv_nucleon.pz = v.Nucleon_Pz();
        fOutput_Injected.lv_nucleon.energy = v.Nucleon_Energy(Particles::FindParticle(fSettings.ReactionChannel.nucleon_pdg)->mass);
        // fill tree //
        fOutputTree_Injected->Fill();
    }
}

// ## Helpers ## //

void Finder::Assign(const View::VecTracks& v, Schema::FlatTrack& out) {
    out.esd_entry = v.EsdEntry();
    out.charge = v.Charge<char>();
    out.state = {v.X(), v.Y(), v.Z(), v.Px(), v.Py(), v.Pz()};
    out.pre_dca_xy = v.PreDCAxy();
    out.pre_dca_z = v.PreDCAz();
    out.tpc_signal = v.TPC_Signal();
    out.n_sigma_pion = v.NSigmaPion();
    out.n_sigma_kaon = v.NSigmaKaon();
    out.n_sigma_proton = v.NSigmaProton();
}

void Finder::Assign(const View::VecV0s& v, Schema::FlatV0& out) {
    out.decay = {v.X(), v.Y(), v.Z()};
    out.lv = {v.Px(), v.Py(), v.Pz(), v.Energy()};
    out.chi2ndf = v.Chi2NDF();
    Assign(v.neg, out.neg);
    out.neg_pca_v0 = {v.Neg_PCAwrtV0_X(), v.Neg_PCAwrtV0_Y(), v.Neg_PCAwrtV0_Z(), v.Neg_PCAwrtV0_Px(), v.Neg_PCAwrtV0_Py(), v.Neg_PCAwrtV0_Pz()};
    Assign(v.pos, out.pos);
    out.pos_pca_v0 = {v.Pos_PCAwrtV0_X(), v.Pos_PCAwrtV0_Y(), v.Pos_PCAwrtV0_Z(), v.Pos_PCAwrtV0_Px(), v.Pos_PCAwrtV0_Py(), v.Pos_PCAwrtV0_Pz()};
}

void Finder::Assign(const View::VecMcTracks& v, Schema::FlatMcTrack& out, bool ascendants_info) {
    out.mc_entry = v.McEntry();
    out.pdg_code = v.PdgCode();
    if (ascendants_info) out.origin = {v.Origin_X(), v.Origin_Y(), v.Origin_Z()};
    out.lv = {v.Px(), v.Py(), v.Pz(), v.Energy()};
    out.is_true = v.IsTrue();
    out.is_secondary = v.IsSecondary();
    out.is_signal = v.IsSignal();
    out.reaction_id = v.ReactionID();
    if (ascendants_info) {
        out.mother_mc_entry = v.Mother_McEntry();
        out.mother_pdg_code = v.Mother_PdgCode();
        out.gm_mc_entry = v.GrandMother_McEntry();
        out.gm_pdg_code = v.GrandMother_PdgCode();
    }
}

void Finder::Assign(const View::VecMcV0s& v, Schema::FlatMcV0& out) {
    out.mc_entry = v.McEntry();
    out.pdg_code = v.PdgCode();
    out.origin = {v.Origin_X(), v.Origin_Y(), v.Origin_Z()};
    out.decay = {v.Decay_X(), v.Decay_Y(), v.Decay_Z()};
    out.lv = {v.Px(), v.Py(), v.Pz(), v.Energy()};
    out.reaction_id = v.ReactionID();
    out.is_true = v.IsTrue();
    out.is_signal = v.IsSignal();
    out.is_secondary = v.IsSecondary();
    out.is_hybrid = v.IsHybrid();
    out.mother_mc_entry = v.Mother_McEntry();
    out.mother_pdg_code = v.Mother_PdgCode();
}

// ## Channel A ZONE ## //

void Finder::FindSexaquarks_ChannelA(bool anti_channel) {

    // determine V0A based on anti-channel or not
    // -- reconstructed
    Schema::VecV0s* Input_V0A = anti_channel ? &fInput_Lambdas : &fInput_AntiLambdas;
    // -- mc
    Schema::VecMcV0s* Input_MC_V0A = nullptr;
    Schema::VecMcTracks* Input_MC_V0A_Neg = nullptr;
    Schema::VecMcTracks* Input_MC_V0A_Pos = nullptr;
    if (fSettings.IsMC) {
        Input_MC_V0A = &fInput_MC_AntiLambdas;
        Input_MC_V0A_Neg = &fInput_MC_AntiLambdas_Neg;
        Input_MC_V0A_Pos = &fInput_MC_AntiLambdas_Pos;
        if (anti_channel) {
            Input_MC_V0A = &fInput_MC_Lambdas;
            Input_MC_V0A_Neg = &fInput_MC_Lambdas_Neg;
            Input_MC_V0A_Pos = &fInput_MC_Lambdas_Pos;
        }
    }

    // -- cut flow hist
    TH1D* hist = anti_channel ? fHist_CutFlow_AntiChannel.get() : fHist_CutFlow.get();

    // loop over all possible pairs of (anti)lambda + K0S //
    std::size_t n_v0a = Input_V0A->decay.x->size();
    std::size_t n_v0b = fInput_KaonsZeroShort.decay.x->size();
    for (std::size_t v0a_idx = 0; v0a_idx < n_v0a; ++v0a_idx) {
        for (std::size_t v0b_idx = 0; v0b_idx < n_v0b; ++v0b_idx) {

            // create views //
            View::VecV0s v0a(Input_V0A, v0a_idx);
            View::VecV0s v0b(&fInput_KaonsZeroShort, v0b_idx);

            // sanity check //
            if (v0a.neg.EsdEntry() == v0b.neg.EsdEntry() || v0a.neg.EsdEntry() == v0b.pos.EsdEntry() || v0a.pos.EsdEntry() == v0b.neg.EsdEntry() ||
                v0a.pos.EsdEntry() == v0b.pos.EsdEntry()) {
                continue;
            }

            // PCAs //
            Seeder::LineLine::Cache pca_cache;
            auto [seed_v0a, seed_v0b] = Seeder::LineLine::FastPCAs(v0a, v0b, &pca_cache);

            // apply cuts (1) //
            if (!FastCuts_ChannelA(seed_v0a, seed_v0b, hist)) continue;

            // PCAs derivatives //
            auto [deriv_v0a, deriv_v0b] = Seeder::LineLine::ComputeDerivatives(pca_cache);

            // fit vertex //
            auto fit = KF::FitVertex(v0a, v0b, {seed_v0a, deriv_v0a}, {seed_v0b, deriv_v0b});
            KF::ChannelA sexa(fit, seed_v0a.pca, seed_v0b.pca, v0a, v0b);

            // apply cuts (2) //
            if (!SlowCuts(sexa, hist)) continue;

            // store //
            Assign(sexa, anti_channel);

            // mc treatment //
            if (fSettings.IsMC) {
                View::VecMcV0s mc_v0a(Input_MC_V0A, v0a_idx);
                Assign(mc_v0a, fOutput_MC_ChannelA_V0A);
                View::VecMcTracks mc_v0a_neg(Input_MC_V0A_Neg, v0a_idx);
                Assign(mc_v0a_neg, fOutput_MC_ChannelA_V0A_Neg, false);
                View::VecMcTracks mc_v0a_pos(Input_MC_V0A_Pos, v0a_idx);
                Assign(mc_v0a_pos, fOutput_MC_ChannelA_V0A_Pos, false);

                View::VecMcV0s mc_v0b(&fInput_MC_KaonsZeroShort, v0b_idx);
                Assign(mc_v0b, fOutput_MC_ChannelA_V0B);
                View::VecMcTracks mc_v0b_neg(&fInput_MC_KaonsZeroShort_Neg, v0b_idx);
                Assign(mc_v0b_neg, fOutput_MC_ChannelA_V0B_Neg, false);
                View::VecMcTracks mc_v0b_pos(&fInput_MC_KaonsZeroShort_Pos, v0b_idx);
                Assign(mc_v0b_pos, fOutput_MC_ChannelA_V0B_Pos, false);

                auto mc_sexa = MC::ChannelA::Create(&fInput_Injected, mc_v0a, mc_v0b);
                Assign(mc_sexa);
            }

            // fill flat tree //
            fOutputTree->Fill();
        }
    }
}

bool Finder::FastCuts_ChannelA(const Seeder::Seed& seed_v0a, const Seeder::Seed& seed_v0b, TH1D* cut_flow_hist) const {

    cut_flow_hist->Fill(0.);

    if (Math::SquaredDistance(seed_v0a.pca.xyz, seed_v0b.pca.xyz) > Cuts::ChannelA::Max_DCAbtwV0s * Cuts::ChannelA::Max_DCAbtwV0s) return false;
    cut_flow_hist->Fill(1.);

    return true;
}

bool Finder::SlowCuts(const KF::ChannelA& sexa, TH1D* cut_flow_hist) const {

    if (sexa.SquaredRadius2D() < Cuts::ChannelA::Min_Radius2D * Cuts::ChannelA::Min_Radius2D) return false;
    cut_flow_hist->Fill(2.);

    if (sexa.SquaredDCA_V0A_SV() > Cuts::ChannelA::Max_DCALaSV * Cuts::ChannelA::Max_DCALaSV) return false;
    cut_flow_hist->Fill(3.);

    if (sexa.SquaredDCA_V0B_SV() > Cuts::ChannelA::Max_DCAK0SV * Cuts::ChannelA::Max_DCAK0SV) return false;
    cut_flow_hist->Fill(4.);

    // if (sexa.CPA_wrt(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) < Cuts::ChannelA::Min_CPAwrtPV) return false;
    // cut_flow_hist->Fill(5.);

    // if (sexa.DCA_wrt(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) > Cuts::ChannelA::Max_DCAwrtPV) return false;
    // cut_flow_hist->Fill(6.);

    // if (sexa.CPA_V0A_SV() < Cuts::ChannelA::Min_La_CPAwrtSV) return false;
    // cut_flow_hist->Fill(7.);

    // if (sexa.CPA_V0B_SV() < Cuts::ChannelA::Min_K0S_CPAwrtSV) return false;
    // cut_flow_hist->Fill(8.);

    return true;
}

void Finder::Assign(const KF::ChannelA& sexa, bool anti_channel) {
    // -- event info
    fOutput_ChannelA.event = fInput_Event;
    // -- candidate info
    fOutput_ChannelA.sv.x = static_cast<float>(sexa.X());
    fOutput_ChannelA.sv.y = static_cast<float>(sexa.Y());
    fOutput_ChannelA.sv.z = static_cast<float>(sexa.Z());
    fOutput_ChannelA.lv.px = static_cast<float>(sexa.Px());
    fOutput_ChannelA.lv.py = static_cast<float>(sexa.Py());
    fOutput_ChannelA.lv.pz = static_cast<float>(sexa.Pz());
    fOutput_ChannelA.lv.energy = static_cast<float>(sexa.E());
    fOutput_ChannelA.chi2ndf = static_cast<float>(sexa.Chi2NDF());
    fOutput_ChannelA.e_minus_nucleon = static_cast<float>(sexa.E_MinusNucleon());
    fOutput_ChannelA.is_anti = anti_channel;
    // -- V0A
    Assign(sexa.V0A, fOutput_ChannelA.v0a);
    fOutput_ChannelA.v0a_pca_sv.x = static_cast<float>(sexa.V0A_PCAwrtSV.X());
    fOutput_ChannelA.v0a_pca_sv.y = static_cast<float>(sexa.V0A_PCAwrtSV.Y());
    fOutput_ChannelA.v0a_pca_sv.z = static_cast<float>(sexa.V0A_PCAwrtSV.Z());
    fOutput_ChannelA.v0a_pca_sv.px = static_cast<float>(sexa.V0A_PCAwrtSV.Px());
    fOutput_ChannelA.v0a_pca_sv.py = static_cast<float>(sexa.V0A_PCAwrtSV.Py());
    fOutput_ChannelA.v0a_pca_sv.pz = static_cast<float>(sexa.V0A_PCAwrtSV.Pz());
    // -- V0B
    Assign(sexa.V0B, fOutput_ChannelA.v0b);
    fOutput_ChannelA.v0b_pca_sv.x = static_cast<float>(sexa.V0B_PCAwrtSV.X());
    fOutput_ChannelA.v0b_pca_sv.y = static_cast<float>(sexa.V0B_PCAwrtSV.Y());
    fOutput_ChannelA.v0b_pca_sv.z = static_cast<float>(sexa.V0B_PCAwrtSV.Z());
    fOutput_ChannelA.v0b_pca_sv.px = static_cast<float>(sexa.V0B_PCAwrtSV.Px());
    fOutput_ChannelA.v0b_pca_sv.py = static_cast<float>(sexa.V0B_PCAwrtSV.Py());
    fOutput_ChannelA.v0b_pca_sv.pz = static_cast<float>(sexa.V0B_PCAwrtSV.Pz());
}

void Finder::Assign(const std::optional<MC::ChannelA>& mc_sexa) {
    if (!mc_sexa) {
        fOutput_MC_ChannelA = {};  // NOTE: dummify
        return;
    }
    fOutput_MC_ChannelA.before.px = mc_sexa->Px();
    fOutput_MC_ChannelA.before.py = mc_sexa->Py();
    fOutput_MC_ChannelA.before.pz = mc_sexa->Pz();
    fOutput_MC_ChannelA.before.energy = mc_sexa->Energy(fSettings.SexaquarkMass);
    fOutput_MC_ChannelA.after.px = mc_sexa->After_Px();
    fOutput_MC_ChannelA.after.py = mc_sexa->After_Py();
    fOutput_MC_ChannelA.after.pz = mc_sexa->After_Pz();
    fOutput_MC_ChannelA.after.energy = mc_sexa->After_Energy();
    fOutput_MC_ChannelA.nucleon.px = mc_sexa->Nucleon_Px();
    fOutput_MC_ChannelA.nucleon.py = mc_sexa->Nucleon_Py();
    fOutput_MC_ChannelA.nucleon.pz = mc_sexa->Nucleon_Pz();
    fOutput_MC_ChannelA.nucleon.energy = mc_sexa->Nucleon_Energy(Particles::FindParticle(fSettings.ReactionChannel.nucleon_pdg)->mass);
    fOutput_MC_ChannelA.sv.x = mc_sexa->SV_X();
    fOutput_MC_ChannelA.sv.y = mc_sexa->SV_Y();
    fOutput_MC_ChannelA.sv.z = mc_sexa->SV_Z();
    fOutput_MC_ChannelA.reaction_id = static_cast<int>(mc_sexa->ReactionID());
    fOutput_MC_ChannelA.is_signal = mc_sexa->IsSignal();
    fOutput_MC_ChannelA.is_hybrid = mc_sexa->IsHybrid();
}

// ## Channel D ZONE ## //

void Finder::FindSexaquarks_ChannelD(bool anti_channel) {

    // determine properties based on anti-channel or not
    // -- (anti)lambda
    //    -- reconstructed
    const Schema::VecV0s* Input_V0s = anti_channel ? &fInput_Lambdas : &fInput_AntiLambdas;
    //    -- mc
    const Schema::VecMcV0s* Input_MC_V0s = nullptr;
    const Schema::VecMcTracks* Input_MC_V0s_Neg = nullptr;
    const Schema::VecMcTracks* Input_MC_V0s_Pos = nullptr;
    if (fSettings.IsMC) {
        Input_MC_V0s = &fInput_MC_AntiLambdas;
        Input_MC_V0s_Neg = &fInput_MC_AntiLambdas_Neg;
        Input_MC_V0s_Pos = &fInput_MC_AntiLambdas_Pos;
        if (anti_channel) {
            Input_MC_V0s = &fInput_MC_Lambdas;
            Input_MC_V0s_Neg = &fInput_MC_Lambdas_Neg;
            Input_MC_V0s_Pos = &fInput_MC_Lambdas_Pos;
        }
    }
    // -- (neg/pos)kaon
    //    -- reconstructed
    const Schema::VecTracks* Input_Kaons = anti_channel ? &fInput_PosKaons : &fInput_NegKaons;
    //    -- mc
    const Schema::VecMcTracks* Input_MC_Kaons = nullptr;
    if (fSettings.IsMC) {
        Input_MC_Kaons = anti_channel ? &fInput_MC_PosKaons : &fInput_MC_NegKaons;
    }
    // -- cut flow hist
    TH1D* hist = anti_channel ? fHist_CutFlow_AntiChannel.get() : fHist_CutFlow.get();

    // loop over all possible pairs of (anti)lambda + (pos/neg)kaon //
    std::size_t n_v0 = Input_V0s->decay.x->size();
    std::size_t n_kaons = Input_Kaons->state.x->size();
    for (std::size_t v0_idx = 0; v0_idx < n_v0; ++v0_idx) {
        for (std::size_t kaon_idx = 0; kaon_idx < n_kaons; ++kaon_idx) {

            // create views //
            View::VecV0s v0(Input_V0s, v0_idx);
            View::VecTracks ka(Input_Kaons, kaon_idx);

            // sanity check //
            if (v0.neg.EsdEntry() == ka.EsdEntry() || v0.pos.EsdEntry() == ka.EsdEntry()) continue;

            // PCAs (1) //
            auto [seed_ka, seed_v0, pca_cache] = Seeder::HelixLine::FastCorrectPCAs(ka, v0, fInput_Event.magnetic_field);

            // apply cuts (1) //
            if (!FastCuts_ChannelD(seed_ka, seed_v0, hist)) continue;

            // PCAs derivatives //
            auto [deriv_v0, deriv_ka] = Seeder::HelixLine::ComputeDerivatives(seed_ka, seed_v0, pca_cache);

            // fit vertex //
            auto fit = KF::FitVertex(ka, v0, Particles::Particle("PosKaon").mass,  //
                                     {seed_ka, deriv_ka}, {seed_v0, deriv_v0},     //
                                     fInput_Event.magnetic_field);
            KF::ChannelD sexa(fit, seed_ka.pca, seed_v0.pca, ka, v0);

            // apply cuts (2) //
            if (!SlowCuts(sexa, hist)) continue;

            // store //
            Assign(sexa, anti_channel);

            if (fSettings.IsMC) {
                View::VecMcV0s mc_v0(Input_MC_V0s, v0_idx);
                Assign(mc_v0, fOutput_MC_ChannelD_V0);
                View::VecMcTracks mc_v0_neg(Input_MC_V0s_Neg, v0_idx);
                Assign(mc_v0_neg, fOutput_MC_ChannelD_V0_Neg, false);
                View::VecMcTracks mc_v0_pos(Input_MC_V0s_Pos, v0_idx);
                Assign(mc_v0_pos, fOutput_MC_ChannelD_V0_Pos, false);

                View::VecMcTracks mc_kaon(Input_MC_Kaons, kaon_idx);
                Assign(mc_kaon, fOutput_MC_ChannelD_Kaon, false);

                auto mc_sexa = MC::ChannelD::Create(&fInput_Injected, mc_v0, mc_kaon);
                Assign(mc_sexa);
            }

            // fill flat tree //
            fOutputTree->Fill();
        }
    }
}

bool Finder::FastCuts_ChannelD(const Seeder::Seed& seed_ka, const Seeder::Seed& seed_v0, TH1D* cut_flow_hist) const {

    cut_flow_hist->Fill(0.);

    if (Math::SquaredDistance(seed_ka.pca.xyz, seed_v0.pca.xyz) > Cuts::ChannelD::Max_DCAKaLa * Cuts::ChannelD::Max_DCAKaLa) return false;
    cut_flow_hist->Fill(1.);

    return true;
}

bool Finder::SlowCuts(const KF::ChannelD& sexa, TH1D* cut_flow_hist) const {

    double sq_radius_2d = sexa.SquaredRadius2D();
    if (sq_radius_2d < Cuts::ChannelD::Min_Radius2D * Cuts::ChannelD::Min_Radius2D ||
        sq_radius_2d > Cuts::ChannelD::Max_Radius2D * Cuts::ChannelD::Max_Radius2D) {
        return false;
    }
    cut_flow_hist->Fill(2.);

    // if (sexa.AbsRapidity_MinusNucleon() > Cuts::ChannelD::AbsMax_Rapidity) return false;  // PENDING: kinematics, affected by Fermi motion
    // cut_flow_hist->Fill(3.);

    // if (sexa.CPA_Vertex(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) < Cuts::ChannelD::Min_CPAwrtPV ||
    // sexa.CPA_Vertex(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) > Cuts::ChannelD::Max_CPAwrtPV) {
    // return false;  // PENDING: kinematics, affected by Fermi motion
    // }
    // cut_flow_hist->Fill(3.);

    if (sexa.SquaredDCA_V0_SV() > Cuts::ChannelD::Max_DCALaSV * Cuts::ChannelD::Max_DCALaSV) return false;
    cut_flow_hist->Fill(4.);

    if (sexa.SquaredDCA_Kaon_SV() > Cuts::ChannelD::Max_DCAKaSV * Cuts::ChannelD::Max_DCAKaSV) return false;
    cut_flow_hist->Fill(5.);

    // if (sexa.DCA_V0Neg_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelD::Max_DCALaNegSV) return false;
    // cut_flow_hist->Fill(6.);

    // if (sexa.DCA_V0Pos_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelD::Max_DCALaPosSV) return false;
    // cut_flow_hist->Fill(7.);

    return true;
}

void Finder::Assign(const KF::ChannelD& sexa, bool anti_channel) {
    fOutput_ChannelD.event = fInput_Event;
    fOutput_ChannelD.sv.x = static_cast<float>(sexa.X());
    fOutput_ChannelD.sv.y = static_cast<float>(sexa.Y());
    fOutput_ChannelD.sv.z = static_cast<float>(sexa.Z());
    fOutput_ChannelD.lv.px = static_cast<float>(sexa.Px());
    fOutput_ChannelD.lv.py = static_cast<float>(sexa.Py());
    fOutput_ChannelD.lv.pz = static_cast<float>(sexa.Pz());
    fOutput_ChannelD.lv.energy = static_cast<float>(sexa.E());
    fOutput_ChannelD.chi2ndf = static_cast<float>(sexa.Chi2NDF());
    fOutput_ChannelD.e_minus_nucleon = static_cast<float>(sexa.E_MinusNucleon());
    fOutput_ChannelD.is_anti = anti_channel;
    // kaon
    Assign(sexa.Kaon, fOutput_ChannelD.kaon);
    fOutput_ChannelD.kaon_pca_sv.x = static_cast<float>(sexa.Kaon_PCAwrtSV.X());
    fOutput_ChannelD.kaon_pca_sv.y = static_cast<float>(sexa.Kaon_PCAwrtSV.Y());
    fOutput_ChannelD.kaon_pca_sv.z = static_cast<float>(sexa.Kaon_PCAwrtSV.Z());
    fOutput_ChannelD.kaon_pca_sv.px = static_cast<float>(sexa.Kaon_PCAwrtSV.Px());
    fOutput_ChannelD.kaon_pca_sv.py = static_cast<float>(sexa.Kaon_PCAwrtSV.Py());
    fOutput_ChannelD.kaon_pca_sv.pz = static_cast<float>(sexa.Kaon_PCAwrtSV.Pz());
    // V0
    Assign(sexa.V0, fOutput_ChannelD.v0);
    fOutput_ChannelD.v0_pca_sv.x = static_cast<float>(sexa.V0_PCAwrtSV.X());
    fOutput_ChannelD.v0_pca_sv.y = static_cast<float>(sexa.V0_PCAwrtSV.Y());
    fOutput_ChannelD.v0_pca_sv.z = static_cast<float>(sexa.V0_PCAwrtSV.Z());
    fOutput_ChannelD.v0_pca_sv.px = static_cast<float>(sexa.V0_PCAwrtSV.Px());
    fOutput_ChannelD.v0_pca_sv.py = static_cast<float>(sexa.V0_PCAwrtSV.Py());
    fOutput_ChannelD.v0_pca_sv.pz = static_cast<float>(sexa.V0_PCAwrtSV.Pz());
}

void Finder::Assign(const std::optional<MC::ChannelD>& sexa) {
    if (!sexa) {
        fOutput_MC_ChannelD = {};  // NOTE: dummify
        return;
    }
    fOutput_MC_ChannelD.before.px = sexa->Px();
    fOutput_MC_ChannelD.before.py = sexa->Py();
    fOutput_MC_ChannelD.before.pz = sexa->Pz();
    fOutput_MC_ChannelD.before.energy = sexa->Energy(fSettings.SexaquarkMass);
    fOutput_MC_ChannelD.after.px = sexa->After_Px();
    fOutput_MC_ChannelD.after.py = sexa->After_Py();
    fOutput_MC_ChannelD.after.pz = sexa->After_Pz();
    fOutput_MC_ChannelD.after.energy = sexa->After_Energy();
    fOutput_MC_ChannelD.nucleon.px = sexa->Nucleon_Px();
    fOutput_MC_ChannelD.nucleon.py = sexa->Nucleon_Py();
    fOutput_MC_ChannelD.nucleon.pz = sexa->Nucleon_Pz();
    fOutput_MC_ChannelD.nucleon.energy = sexa->Nucleon_Energy(Particles::FindParticle(fSettings.ReactionChannel.nucleon_pdg)->mass);
    fOutput_MC_ChannelD.sv.x = sexa->SV_X();
    fOutput_MC_ChannelD.sv.y = sexa->SV_Y();
    fOutput_MC_ChannelD.sv.z = sexa->SV_Z();
    fOutput_MC_ChannelD.reaction_id = static_cast<int>(sexa->ReactionID());
    fOutput_MC_ChannelD.is_signal = sexa->IsSignal();
    fOutput_MC_ChannelD.is_hybrid = sexa->IsHybrid();
}

// ## END OF CYCLES ## //

void Finder::EndOfAnalysis() {

    Logger::Info(__FUNCTION__, "The following objects have been written into TFile {}:", fSettings.PathOutputFile);

    // write trees //

    if (fSettings.IsMC) {
        fOutputTree_Injected->Write();
        Logger::Info(__FUNCTION__, "- TTree \"{}\"", fOutputTree_Injected->GetName());
    }
    fOutputTree->Write();
    Logger::Info(__FUNCTION__, "- TTree \"{}\"", fOutputTree->GetName());

    // write histograms //

    // -- event counter
    fHist_EventCounter->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_EventCounter->GetName());
    // -- cut flows
    fHist_CutFlow->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow->GetName());
    fHist_CutFlow_AntiChannel->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow_AntiChannel->GetName());

    // reset branch addresses //

    fInputChain_PackedEvents->ResetBranchAddresses();
    if (fSettings.IsMC) fOutputTree_Injected->ResetBranchAddresses();
    fOutputTree->ResetBranchAddresses();

    Logger::Info(__FUNCTION__, "All done.");
}

}  // namespace Tree2Secondaries
