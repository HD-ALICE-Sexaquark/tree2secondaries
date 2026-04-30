#include "Packager/Packager.hxx"

#include <cstddef>
#include <filesystem>

#include "App/DB_Particles.hxx"
#include "App/Logger.hxx"
#include "KalmanFitter/KalmanFitterV0.hxx"
#include "Math/Constants.hxx"
#include "MonteCarlo/MonteCarloParticle.hxx"
#include "MonteCarlo/MonteCarloV0.hxx"
#include "Packager/PackagerCuts.hxx"
#include "Schema/SchemaVectorMcTracks.hxx"
#include "View/ViewVectorTracks.hxx"
#if !T2S_LEGACY_KF
#include "KalmanFitter/BaseKalmanFitter.hxx"
#include "Seeder/SeederHelixHelix.hxx"
#else
#include "Legacy/LegacyFitter.hxx"
#include "Legacy/LegacyHelixHelix.hxx"
#endif

namespace Tree2Secondaries {

namespace KF = KalmanFitter;
namespace MC = MonteCarlo;

bool Packager::Initialize() {

    fInputChain_Events = std::make_unique<TChain>(Const::TreeName_Events.data());
    for (const auto& path : fSettings.PathInputFiles) {
        if (fInputChain_Events->Add(path.c_str()) == 0) {
            Logger::Error(__FUNCTION__, "Couldn't add TFile {}", path);
        }
    }
    if (!fInputChain_Events->GetEntries()) {
        Logger::Error(__FUNCTION__, "Couldn't read any entry.");
        return false;
    }
    Logger::Info(__FUNCTION__, "TChain \"{}\" loaded successfully with {} trees and {} total entries.", fInputChain_Events->GetName(),
                 fInputChain_Events->GetNtrees(), fInputChain_Events->GetEntries());

    ReadInputBranches();

    if (!PrepareOutputFile()) return false;

    PrepareOutputHistograms();

    if (!PrepareOutputTree()) return false;
    CreateOutputBranches();

    Logger::Info(__FUNCTION__, "Packager initialized successfully.");

    return true;
}

// ## INPUT ZONE ## //

void Packager::ReadInputBranches() {
    // by default, turn off all branches
    fInputChain_Events->SetBranchStatus("*", false);
    // connect input branches to memory
    fInput_Event.ReadBranches(fInputChain_Events.get(), fSettings.IsMC);
    if (fSettings.IsMC) {
        fInput_Injected.ReadBranches(fInputChain_Events.get(), false);
        fInput_MC.ReadBranches(fInputChain_Events.get());
    }
    fInput_Tracks.ReadBranches(fInputChain_Events.get(), fSettings.IsMC);
}

// ## OUTPUT ZONE ## //

bool Packager::PrepareOutputFile() {

    const std::filesystem::path output_path(fSettings.PathOutputFile);
    if (output_path.has_parent_path()) std::filesystem::create_directories(output_path.parent_path());

    fOutputFile = std::unique_ptr<TFile>(TFile::Open(fSettings.PathOutputFile.c_str(), "RECREATE"));
    if (!fOutputFile) {
        Logger::Error(__FUNCTION__, "Couldn't create TFile {}", fSettings.PathOutputFile);
        return false;
    }

    return true;
}

bool Packager::PrepareOutputTree() {

    fOutputTree = std::make_unique<TTree>(Const::TreeName_PackedEvents.data(), "Packed Events");
    if (fOutputTree == nullptr) {
        Logger::Error(__FUNCTION__, "Couldn't create TTree \"{}\"", Const::TreeName_PackedEvents);
        return false;
    }

    return true;
}

void Packager::CreateOutputBranches() {
    // event info //
    fOutput_Event.CreateBranches(fOutputTree.get(), fSettings.IsMC);
    // injected anti-sexaquark reaction info //
    if (fSettings.IsMC) fOutput_Injected.CreateBranches(fOutputTree.get(), true);
    // anti-lambdas //
    fOutput_AntiLambdas.CreateBranches(fOutputTree.get(), Particles::Particle("AntiLambda").acronym);
    if (fSettings.IsMC) {
        fOutput_MC_AntiLambdas.CreateBranches(fOutputTree.get(), Particles::Particle("AntiLambda").acronym);
        fOutput_MC_AntiLambdas_Neg.CreateBranches(fOutputTree.get(), false, std::format("{}_Neg", Particles::Particle("AntiLambda").acronym));
        fOutput_MC_AntiLambdas_Pos.CreateBranches(fOutputTree.get(), false, std::format("{}_Pos", Particles::Particle("AntiLambda").acronym));
    }
    // lambdas //
    fOutput_Lambdas.CreateBranches(fOutputTree.get(), Particles::Particle("Lambda").acronym);
    if (fSettings.IsMC) {
        fOutput_MC_Lambdas.CreateBranches(fOutputTree.get(), Particles::Particle("Lambda").acronym);
        fOutput_MC_Lambdas_Neg.CreateBranches(fOutputTree.get(), false, std::format("{}_Neg", Particles::Particle("Lambda").acronym));
        fOutput_MC_Lambdas_Pos.CreateBranches(fOutputTree.get(), false, std::format("{}_Pos", Particles::Particle("Lambda").acronym));
    }
    // K0S //
    fOutput_KaonsZeroShort.CreateBranches(fOutputTree.get(), Particles::Particle("KaonZeroShort").acronym);
    if (fSettings.IsMC) {
        fOutput_MC_KaonsZeroShort.CreateBranches(fOutputTree.get(), Particles::Particle("KaonZeroShort").acronym);
        fOutput_MC_KaonsZeroShort_Neg.CreateBranches(fOutputTree.get(), false, std::format("{}_Neg", Particles::Particle("KaonZeroShort").acronym));
        fOutput_MC_KaonsZeroShort_Pos.CreateBranches(fOutputTree.get(), false, std::format("{}_Pos", Particles::Particle("KaonZeroShort").acronym));
    }
    // neg kaons //
    fOutput_NegKaons.CreateBranches(fOutputTree.get(), false, Particles::Particle("NegKaon").acronym);
    if (fSettings.IsMC) fOutput_MC_NegKaons.CreateBranches(fOutputTree.get(), true, Particles::Particle("NegKaon").acronym);
    // pos kaons //
    fOutput_PosKaons.CreateBranches(fOutputTree.get(), false, Particles::Particle("PosKaon").acronym);
    if (fSettings.IsMC) fOutput_MC_PosKaons.CreateBranches(fOutputTree.get(), true, Particles::Particle("PosKaon").acronym);
}

void Packager::PrepareOutputHistograms() {

    // event counter //
    fHist_EventCounter = std::make_unique<TH1D>("N_Events", ";;N_Events", 1, 0, 1);

    constexpr int x_nbins = 20;
    constexpr float x_min = 0.;
    constexpr float x_max = 20.;
    std::string_view hist_title = ";Cut N;N Passed Cut";

    fHist_CutFlow_AntiProton = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", Particles::Particle("AntiProton").acronym).c_str(), hist_title.data(), x_nbins, x_min, x_max);
    fHist_CutFlow_Proton = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", Particles::Particle("Proton").acronym).c_str(), hist_title.data(), x_nbins, x_min, x_max);
    fHist_CutFlow_NegKaon = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", Particles::Particle("NegKaon").acronym).c_str(), hist_title.data(), x_nbins, x_min, x_max);
    fHist_CutFlow_PosKaon = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", Particles::Particle("PosKaon").acronym).c_str(), hist_title.data(), x_nbins, x_min, x_max);
    fHist_CutFlow_PiMinus = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", Particles::Particle("PiMinus").acronym).c_str(), hist_title.data(), x_nbins, x_min, x_max);
    fHist_CutFlow_PiPlus = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", Particles::Particle("PiPlus").acronym).c_str(), hist_title.data(), x_nbins, x_min, x_max);

    fHist_CutFlow_AntiLambda = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", Particles::Particle("AntiLambda").acronym).c_str(), hist_title.data(), x_nbins, x_min, x_max);
    fHist_CutFlow_Lambda = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", Particles::Particle("Lambda").acronym).c_str(), hist_title.data(), x_nbins, x_min, x_max);
    fHist_CutFlow_KaonZeroShort = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", Particles::Particle("KaonZeroShort").acronym).c_str(), hist_title.data(), x_nbins, x_min, x_max);
}

// ## Event ZONE ## //

void Packager::ProcessEvent() {
    fOutput_Event = fInput_Event;
    fHist_EventCounter->Fill(0.);
}

// ## MC/Injected ZONE ## //

// Loop over all MC particles.
// Select Primary Particles (particles with no mother), generated via the Sexaquark-Reaction Generator, and with valid Reaction ID;
// and store their origin vertex as the coordinates for this particular secondary vertex.
void Packager::ProcessInjected() {

    std::vector<float> sv_x(Const::NInjectedPerEvent, Const::DummyFloat);
    std::vector<float> sv_y(Const::NInjectedPerEvent, Const::DummyFloat);
    std::vector<float> sv_z(Const::NInjectedPerEvent, Const::DummyFloat);
    std::vector<char> sv_found(Const::NInjectedPerEvent, 0);

    for (std::size_t mc_idx = 0; mc_idx < NumberMC(); ++mc_idx) {

        auto mc = MC::Particle::Create(&fInput_MC, mc_idx, fSettings.ReactionChannel, Const::DummyInt);
        if (!mc || !mc->IsFirstGenSignal()) continue;

        auto reaction_idx = static_cast<std::size_t>(mc->ReactionID() - Const::ReactionID_Offset);
        if (sv_found[reaction_idx] == 1) continue;

        sv_x[reaction_idx] = mc->Origin_X();
        sv_y[reaction_idx] = mc->Origin_Y();
        sv_z[reaction_idx] = mc->Origin_Z();
        sv_found[reaction_idx] = 1;
    }

    // store //

    fOutput_Injected.mom = fInput_Injected.mom;
    fOutput_Injected.mom_nucleon = fInput_Injected.mom_nucleon;
    for (std::size_t i = 0; i < Const::NInjectedPerEvent; ++i) {
        fOutput_Injected.sv.x->push_back(sv_x[i]);
        fOutput_Injected.sv.y->push_back(sv_y[i]);
        fOutput_Injected.sv.z->push_back(sv_z[i]);
    }
    fOutput_Injected.reaction_id = fInput_Injected.reaction_id;
}

// ## Tracks ZONE ## //

// Filter and group tracks into indices vectors, according to their respective species.
void Packager::ProcessTracks() {
    for (std::size_t t_idx = 0; t_idx < NumberTracks(); ++t_idx) {
        // create track view //
        View::VecTracks track(&fInput_Tracks, t_idx);
        // PID and pre-selection //
        if (track.Charge<char>() < 0) {
            if (PassesProtonCuts(track, fHist_CutFlow_AntiProton.get())) {
                fIndices_AntiProtons.push_back(t_idx);
            }
            if (PassesKaonCuts(track, fHist_CutFlow_NegKaon.get())) {
                fIndices_NegKaons.push_back(t_idx);
            }
            if (PassesPionCuts(track, fHist_CutFlow_PiMinus.get())) {
                fIndices_PiMinus.push_back(t_idx);
            }
        }
        if (track.Charge<char>() > 0) {
            if (PassesProtonCuts(track, fHist_CutFlow_Proton.get())) {
                fIndices_Protons.push_back(t_idx);
            }
            if (PassesKaonCuts(track, fHist_CutFlow_PosKaon.get())) {
                fIndices_PosKaons.push_back(t_idx);
            }
            if (PassesPionCuts(track, fHist_CutFlow_PiPlus.get())) {
                fIndices_PiPlus.push_back(t_idx);
            }
        }
    }  // end of loop over tracks

#if T2S_DEBUG
    Logger::Debug(__FUNCTION__, "n_antiprotons = {}", fIndices_AntiProtons.size());
    Logger::Debug(__FUNCTION__, "n_protons     = {}", fIndices_Protons.size());
    Logger::Debug(__FUNCTION__, "n_negkaons    = {}", fIndices_NegKaons.size());
    Logger::Debug(__FUNCTION__, "n_poskaons    = {}", fIndices_PosKaons.size());
    Logger::Debug(__FUNCTION__, "n_piminus     = {}", fIndices_PiMinus.size());
    Logger::Debug(__FUNCTION__, "n_piplus      = {}", fIndices_PiPlus.size());
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
}

// NOTE: exclusive to secondary charged kaons.
void Packager::PackTracks(const Particles::Definition& pid) {

    // determine rules based on particle species //
    std::vector<std::size_t>* vec = nullptr;
    Schema::VecTracks* out = nullptr;
    Schema::VecMcTracks* mc_out = nullptr;
    switch (pid.pdg_code) {
        case Particles::Particle("NegKaon").pdg_code: {
            vec = &fIndices_NegKaons;
            out = &fOutput_NegKaons;
            if (fSettings.IsMC) mc_out = &fOutput_MC_NegKaons;
            break;
        }
        case Particles::Particle("PosKaon").pdg_code: {
            vec = &fIndices_PosKaons;
            out = &fOutput_PosKaons;
            if (fSettings.IsMC) mc_out = &fOutput_MC_PosKaons;
            break;
        }
        default:
            return;
    }

    // loop over selected tracks //
    for (const std::size_t& t_idx : *vec) {
        View::VecTracks track(&fInput_Tracks, t_idx);

        Store(track, *out);

        if (fSettings.IsMC) {
            auto mc = MC::Particle::Create(&fInput_MC, track.McEntry(), fSettings.ReactionChannel, pid.pdg_code);
            Store(mc, *mc_out, true);
        }
    }  // end of loop over selected tracks
}

void Packager::Store(const View::VecTracks& v, Schema::VecTracks& df) {
    df.esd_entry->push_back(v.EsdEntry());
    df.state.x->push_back(v.X());
    df.state.y->push_back(v.Y());
    df.state.z->push_back(v.Z());
    df.state.px->push_back(v.Px());
    df.state.py->push_back(v.Py());
    df.state.pz->push_back(v.Pz());
    df.charge->push_back(v.Charge<char>());
    df.pre_dca_xy->push_back(v.PreDCAxy());
    df.pre_dca_z->push_back(v.PreDCAz());
    df.tpc_signal->push_back(v.TPC_Signal());
    df.n_sigma_pion->push_back(v.NSigmaPion());
    df.n_sigma_kaon->push_back(v.NSigmaKaon());
    df.n_sigma_proton->push_back(v.NSigmaProton());
    df.cov.Push(v.Cov());
}

void Packager::Store(const std::optional<MC::Particle>& mc, Schema::VecMcTracks& df, bool ascendants_info) {
    if (!mc) {
        df.PushDummy(ascendants_info);
        return;
    }
    df.mc_entry->push_back(mc->McEntry());
    df.pdg_code->push_back(mc->PdgCode());
    if (ascendants_info) {
        df.origin.x->push_back(mc->Origin_X());
        df.origin.y->push_back(mc->Origin_Y());
        df.origin.z->push_back(mc->Origin_Z());
    }
    df.lv.px->push_back(mc->Px());
    df.lv.py->push_back(mc->Py());
    df.lv.pz->push_back(mc->Pz());
    df.lv.energy->push_back(mc->Energy());
    df.is_true->push_back(static_cast<char>(mc->IsTrue()));
    df.is_secondary->push_back(static_cast<char>(mc->IsSecondary()));
    df.is_signal->push_back(static_cast<char>(mc->IsSignal()));
    df.reaction_id->push_back(mc->ReactionID());
    if (ascendants_info) {
        df.mother_mc_entry->push_back(mc->Mother_McEntry());
        df.mother_pdg_code->push_back(mc->Mother_PdgCode());
        df.gm_mc_entry->push_back(mc->GrandMother_McEntry());
        df.gm_pdg_code->push_back(mc->GrandMother_PdgCode());
    }
}

bool Packager::PassesProtonCuts(const View::VecTracks& v, TH1D* cut_flow_hist) const {

    cut_flow_hist->Fill(0.);
    if (std::abs(v.NSigmaProton()) > Cuts::Proton::AbsMax_NSigmaProton) return false;
    cut_flow_hist->Fill(1.);

    return true;
}

bool Packager::PassesKaonCuts(const View::VecTracks& v, TH1D* cut_flow_hist) const {

    cut_flow_hist->Fill(0.);
    if (std::abs(v.NSigmaKaon()) > Cuts::Kaon::AbsMax_NSigmaKaon) return false;
    cut_flow_hist->Fill(1.);

    return true;
}

bool Packager::PassesPionCuts(const View::VecTracks& v, TH1D* cut_flow_hist) const {

    cut_flow_hist->Fill(0.);
    if (std::abs(v.NSigmaPion()) > Cuts::Pion::AbsMax_NSigmaPion) return false;
    cut_flow_hist->Fill(1.);

    return true;
}

// ## V0s ZONE ## //

void Packager::FindV0s(const Particles::Definition& pid) {

    // determine rules based on V0 species //
    Schema::VecV0s* out = nullptr;
    Schema::VecMcV0s* mc_out = nullptr;
    Schema::VecMcTracks* mc_neg_out = nullptr;
    Schema::VecMcTracks* mc_pos_out = nullptr;
    const std::vector<std::size_t>* vec_neg = &fIndices_PiMinus;
    const std::vector<std::size_t>* vec_pos = &fIndices_PiPlus;
    switch (pid.pdg_code) {
        case Particles::Particle("AntiLambda").pdg_code: {
            out = &fOutput_AntiLambdas;
            if (fSettings.IsMC) {
                mc_out = &fOutput_MC_AntiLambdas;
                mc_neg_out = &fOutput_MC_AntiLambdas_Neg;
                mc_pos_out = &fOutput_MC_AntiLambdas_Pos;
            }
            vec_neg = &fIndices_AntiProtons;
            break;
        }
        case Particles::Particle("Lambda").pdg_code: {
            out = &fOutput_Lambdas;
            if (fSettings.IsMC) {
                mc_out = &fOutput_MC_Lambdas;
                mc_neg_out = &fOutput_MC_Lambdas_Neg;
                mc_pos_out = &fOutput_MC_Lambdas_Pos;
            }
            vec_pos = &fIndices_Protons;
            break;
        }
        case Particles::Particle("KaonZeroShort").pdg_code: {
            out = &fOutput_KaonsZeroShort;
            if (fSettings.IsMC) {
                mc_out = &fOutput_MC_KaonsZeroShort;
                mc_neg_out = &fOutput_MC_KaonsZeroShort_Neg;
                mc_pos_out = &fOutput_MC_KaonsZeroShort_Pos;
            }
            break;
        }
        default: {
            Logger::Error(__FUNCTION__, "Invalid PID {} for a V0.", pid.name);
            return;
        }
    }
    auto neg_pid = Particles::FindParticle(pid.daughters_pdg[0]);
    auto pos_pid = Particles::FindParticle(pid.daughters_pdg[1]);

    // loop over all possible pairs of tracks //
    for (const auto& neg_idx : *vec_neg) {
        for (const auto& pos_idx : *vec_pos) {

            // sanity check //
            if (neg_idx == pos_idx) continue;

#if T2S_DEBUG
            if (pid == PID_V0::AntiLambda && (esd_neg != 725 || esd_pos != 739)) continue;
            if (pid == PID_V0::Lambda) continue;
            if (pid == PID_V0::KaonZeroShort && (esd_neg != 129 || esd_pos != 425)) continue;
#endif

            // get views //
            View::VecTracks neg(&fInput_Tracks, neg_idx);
            View::VecTracks pos(&fInput_Tracks, pos_idx);

#if T2S_LEGACY_KF
            // PCAs //
            auto [seed_neg, seed_pos, deriv_neg, deriv_pos] = Legacy::HelixHelix::FullPCAs(neg, pos, mass_neg, mass_pos, fInput_Event.magnetic_field);

            // fit vertex //
            auto l_fit = Legacy::Fit(neg, pos, mass_neg, mass_pos, fInput_Event.magnetic_field);
            auto fit = KF::Particle::FromLegacy(l_fit);
#else
            // PCAs //
            auto [seed_neg, seed_pos, pca_cache] = Seeder::HelixHelix::FastCorrectPCAs(neg, pos, fInput_Event.magnetic_field);

            // apply cuts (1) //
            if (!FastCuts(seed_neg, seed_pos, pid)) continue;

            // PCAs derivatives //
            auto [deriv_neg, deriv_pos] = Seeder::HelixHelix::ComputeDerivatives(seed_neg, seed_pos, pca_cache);

            // fit vertex //
            auto fit =
                KF::FitVertex(neg, pos, neg_pid->mass, pos_pid->mass, {seed_neg, deriv_neg}, {seed_pos, deriv_pos}, fInput_Event.magnetic_field);

            // create composite particle //
            KF::V0 v0(fit, seed_neg.pca, seed_pos.pca, neg, pos);

            // apply cuts (2) //
            if (!SlowCuts(v0, pid)) continue;

            Store(v0, *out);
#endif
            // store mc //
            if (fSettings.IsMC) {
                auto mc_neg = MC::Particle::Create(&fInput_MC, neg.McEntry(), fSettings.ReactionChannel, neg_pid->pdg_code);
                Store(mc_neg, *mc_neg_out, false);
                auto mc_pos = MC::Particle::Create(&fInput_MC, pos.McEntry(), fSettings.ReactionChannel, pos_pid->pdg_code);
                Store(mc_pos, *mc_pos_out, false);
                auto mc_v0 = MC::V0::Create(mc_neg, mc_pos, fSettings.ReactionChannel, pid.pdg_code);
                Store(mc_v0, *mc_out);
            }
        }  // end of loop over pos
    }  // end of loop over neg
}

bool Packager::FastCuts_SecondaryLambdas(const Seeder::Seed& seed_neg, const Seeder::Seed& seed_pos, TH1D* cut_flow_hist) const {
    cut_flow_hist->Fill(0.);

    if (Math::SquaredDistance(seed_neg.pca.GetXYZ_AsROOT(), seed_pos.pca.GetXYZ_AsROOT()) > Cuts::Lambda::Max_DCAbtwDau * Cuts::Lambda::Max_DCAbtwDau)
        return false;
    cut_flow_hist->Fill(1.);

    return true;
}

bool Packager::FastCuts_KaonZeroShort(const Seeder::Seed& seed_neg, const Seeder::Seed& seed_pos, TH1D* cut_flow_hist) const {
    cut_flow_hist->Fill(0.);

    if (Math::SquaredDistance(seed_neg.pca.GetXYZ_AsROOT(), seed_pos.pca.GetXYZ_AsROOT()) >
        Cuts::KaonZeroShort::Max_DCAbtwDau * Cuts::KaonZeroShort::Max_DCAbtwDau)
        return false;
    cut_flow_hist->Fill(1.);

    return true;
}

bool Packager::SlowCuts_SecondaryLambdas(const KF::V0& v0, TH1D* cut_flow_hist) const {

    double mass = v0.Mass().value_or(Const::DummyDouble);  // cached
    if (mass < Cuts::Lambda::Min_Mass || mass > Cuts::Lambda::Max_Mass) return false;
    cut_flow_hist->Fill(1.);

    if (v0.SquaredRadius2D() < Cuts::Lambda::Min_Radius2D * Cuts::Lambda::Min_Radius2D) return false;
    cut_flow_hist->Fill(3.);

    if (v0.SquaredDCA_Neg_V0() > Cuts::Lambda::Max_DCAnegV0 * Cuts::Lambda::Max_DCAnegV0) return false;
    cut_flow_hist->Fill(4.);

    if (v0.SquaredDCA_Pos_V0() > Cuts::Lambda::Max_DCAposV0 * Cuts::Lambda::Max_DCAposV0) return false;
    cut_flow_hist->Fill(5.);

    // if (v0.Pt() < Cuts::Lambda::Min_Pt) return false; // PENDING
    // cut_flow_hist->Fill(6.); // PENDING

    if (std::abs(v0.Rapidity()) > Cuts::Lambda::AbsMax_Rapidity) return false;
    cut_flow_hist->Fill(7.);

    if (v0.AbsArmQtOverAlpha() > Cuts::Lambda::AbsMax_ArmQtOverAlpha) return false;
    cut_flow_hist->Fill(8.);

    double cpa_wrt_pv = v0.CPA_Vertex(fInput_Event.pv.x, fInput_Event.pv.y, fInput_Event.pv.z);  // cached
    if (cpa_wrt_pv < Cuts::Lambda::Min_CPAwrtPV || cpa_wrt_pv > Cuts::Lambda::Max_CPAwrtPV) return false;
    cut_flow_hist->Fill(9.);

    if (v0.DCA_Vertex(fInput_Event.pv.x, fInput_Event.pv.y, fInput_Event.pv.z) < Cuts::Lambda::Min_DCAwrtPV) return false;
    cut_flow_hist->Fill(10.);

    return true;
}

bool Packager::SlowCuts_KaonZeroShort(const KF::V0& v0, TH1D* cut_flow_hist) const {

    // if (v0.Pt() < Cuts::KaonZeroShort::Min_Pt) return false; // PENDING
    // cut_flow_hist->Fill(2.); // PENDING

    double mass = v0.Mass().value_or(Const::DummyDouble);  // cached
    if (mass < Cuts::KaonZeroShort::Min_Mass || mass > Cuts::KaonZeroShort::Max_Mass) return false;
    cut_flow_hist->Fill(3.);

    if (std::abs(v0.Rapidity()) > Cuts::KaonZeroShort::AbsMax_Rapidity) return false;
    cut_flow_hist->Fill(4.);

    if (v0.SquaredRadius2D() < Cuts::KaonZeroShort::Min_Radius2D * Cuts::KaonZeroShort::Min_Radius2D) return false;
    cut_flow_hist->Fill(5.);

    if (v0.SquaredDCA_Neg_V0() > Cuts::KaonZeroShort::Max_DCAnegV0 * Cuts::KaonZeroShort::Max_DCAnegV0) return false;
    cut_flow_hist->Fill(6.);

    if (v0.SquaredDCA_Pos_V0() > Cuts::KaonZeroShort::Max_DCAposV0 * Cuts::KaonZeroShort::Max_DCAposV0) return false;
    cut_flow_hist->Fill(7.);

    double cpa_wrt_pv = v0.CPA_Vertex(fInput_Event.pv.x, fInput_Event.pv.y, fInput_Event.pv.z);  // cached
    if (cpa_wrt_pv < Cuts::KaonZeroShort::Min_CPAwrtPV || cpa_wrt_pv > Cuts::KaonZeroShort::Max_CPAwrtPV) return false;
    cut_flow_hist->Fill(8.);

    if (v0.DCA_Vertex(fInput_Event.pv.x, fInput_Event.pv.y, fInput_Event.pv.z) < Cuts::KaonZeroShort::Min_DCAwrtPV) return false;
    cut_flow_hist->Fill(9.);

    return true;
}

void Packager::Store(const KF::V0& v0, Schema::VecV0s& df) {
    df.decay.x->push_back(static_cast<float>(v0.X()));
    df.decay.y->push_back(static_cast<float>(v0.Y()));
    df.decay.z->push_back(static_cast<float>(v0.Z()));
    df.lv.px->push_back(static_cast<float>(v0.Px()));
    df.lv.py->push_back(static_cast<float>(v0.Py()));
    df.lv.pz->push_back(static_cast<float>(v0.Pz()));
    df.lv.energy->push_back(static_cast<float>(v0.E()));
    df.chi2ndf->push_back(static_cast<float>(v0.Chi2NDF()));
    df.cov.Push(v0.Cov<7>());
    // -- neg daughter
    Store(v0.Neg, df.neg);
    df.neg_pca_v0.x->push_back(static_cast<float>(v0.Neg_PCAwrtV0.X()));
    df.neg_pca_v0.y->push_back(static_cast<float>(v0.Neg_PCAwrtV0.Y()));
    df.neg_pca_v0.z->push_back(static_cast<float>(v0.Neg_PCAwrtV0.Z()));
    df.neg_pca_v0.px->push_back(static_cast<float>(v0.Neg_PCAwrtV0.Px()));
    df.neg_pca_v0.py->push_back(static_cast<float>(v0.Neg_PCAwrtV0.Py()));
    df.neg_pca_v0.pz->push_back(static_cast<float>(v0.Neg_PCAwrtV0.Pz()));
    // -- pos daughter
    Store(v0.Pos, df.pos);
    df.pos_pca_v0.x->push_back(static_cast<float>(v0.Pos_PCAwrtV0.X()));
    df.pos_pca_v0.y->push_back(static_cast<float>(v0.Pos_PCAwrtV0.Y()));
    df.pos_pca_v0.z->push_back(static_cast<float>(v0.Pos_PCAwrtV0.Z()));
    df.pos_pca_v0.px->push_back(static_cast<float>(v0.Pos_PCAwrtV0.Px()));
    df.pos_pca_v0.py->push_back(static_cast<float>(v0.Pos_PCAwrtV0.Py()));
    df.pos_pca_v0.pz->push_back(static_cast<float>(v0.Pos_PCAwrtV0.Pz()));
}

void Packager::Store(const std::optional<MC::V0>& mc_v0, Schema::VecMcV0s& df) {
    if (!mc_v0) {
        df.PushDummy();
        return;
    }
    df.mc_entry->push_back(mc_v0->McEntry());
    df.pdg_code->push_back(mc_v0->PdgCode());
    df.origin.x->push_back(mc_v0->Origin_X());
    df.origin.y->push_back(mc_v0->Origin_Y());
    df.origin.z->push_back(mc_v0->Origin_Z());
    df.decay.x->push_back(mc_v0->Decay_X());
    df.decay.y->push_back(mc_v0->Decay_Y());
    df.decay.z->push_back(mc_v0->Decay_Z());
    df.lv.px->push_back(mc_v0->Px());
    df.lv.py->push_back(mc_v0->Py());
    df.lv.pz->push_back(mc_v0->Pz());
    df.lv.energy->push_back(mc_v0->Energy());
    df.is_true->push_back(static_cast<char>(mc_v0->IsTrue()));
    df.is_secondary->push_back(static_cast<char>(mc_v0->IsSecondary()));
    df.is_signal->push_back(static_cast<char>(mc_v0->IsSignal()));
    df.is_hybrid->push_back(static_cast<char>(mc_v0->IsHybrid()));
    df.reaction_id->push_back(mc_v0->ReactionID());
    df.mother_mc_entry->push_back(mc_v0->Mother_McEntry());
    df.mother_pdg_code->push_back(mc_v0->Mother_PdgCode());
}

// ## END OF CYCLES ## //

void Packager::EndOfEvent() {

    // fill tree//

    fOutputTree->Fill();

    // clear temporary indices containers //

    fIndices_AntiProtons.clear();
    fIndices_Protons.clear();
    fIndices_NegKaons.clear();
    fIndices_PosKaons.clear();
    fIndices_PiMinus.clear();
    fIndices_PiPlus.clear();

    // clear output vector branches //

    if (fSettings.IsMC) fOutput_Injected.ClearBranches(true);

    fOutput_AntiLambdas.ClearBranches();
    fOutput_Lambdas.ClearBranches();
    fOutput_KaonsZeroShort.ClearBranches();

    fOutput_NegKaons.ClearBranches(false);
    fOutput_PosKaons.ClearBranches(false);

    if (fSettings.IsMC) {
        fOutput_MC_AntiLambdas.ClearBranches();
        fOutput_MC_AntiLambdas_Neg.ClearBranches(false);
        fOutput_MC_AntiLambdas_Pos.ClearBranches(false);

        fOutput_MC_Lambdas.ClearBranches();
        fOutput_MC_Lambdas_Neg.ClearBranches(false);
        fOutput_MC_Lambdas_Pos.ClearBranches(false);

        fOutput_MC_KaonsZeroShort.ClearBranches();
        fOutput_MC_KaonsZeroShort_Neg.ClearBranches(false);
        fOutput_MC_KaonsZeroShort_Pos.ClearBranches(false);

        fOutput_MC_NegKaons.ClearBranches(true);
        fOutput_MC_PosKaons.ClearBranches(true);
    }
}

void Packager::EndOfAnalysis() {

    Logger::Info(__FUNCTION__, "The following objects have been written into TFile {}:", fSettings.PathOutputFile);

    // write tree //

    fOutputTree->Write();
    Logger::Info(__FUNCTION__, "- TTree \"{}\"", fOutputTree->GetName());

    // write histograms //

    // -- event counter
    fHist_EventCounter->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_EventCounter->GetName());

    // -- selected tracks
    fHist_CutFlow_AntiProton->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow_AntiProton->GetName());
    fHist_CutFlow_Proton->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow_Proton->GetName());
    fHist_CutFlow_NegKaon->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow_NegKaon->GetName());
    fHist_CutFlow_PosKaon->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow_PosKaon->GetName());
    fHist_CutFlow_PiMinus->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow_PiMinus->GetName());
    fHist_CutFlow_PiPlus->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow_PiPlus->GetName());

    // -- found v0s
    fHist_CutFlow_AntiLambda->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow_AntiLambda->GetName());
    fHist_CutFlow_Lambda->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow_Lambda->GetName());
    fHist_CutFlow_KaonZeroShort->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow_KaonZeroShort->GetName());

    fInputChain_Events->ResetBranchAddresses();
    fOutputTree->ResetBranchAddresses();

    Logger::Info(__FUNCTION__, "All done.");
}

}  // namespace Tree2Secondaries
