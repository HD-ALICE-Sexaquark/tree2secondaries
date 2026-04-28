#include "Packager/Packager.hxx"

#include <array>
#include <cstddef>
#include <filesystem>

#include "App/Logger.hxx"
#include "KalmanFitter/KalmanFitterV0.hxx"
#include "Math/Constants.hxx"
#include "Packager/PackagerCuts.hxx"
#include "Storage/ReadWrite/ReadWriteSchemaFlat.hxx"
#include "Storage/ReadWrite/ReadWriteSchemaVector.hxx"
#include "Truth/TruthTrack.hxx"
#include "Truth/TruthV0.hxx"
#include "View/BaseView.hxx"
#include "View/MC/ViewMcParticle.hxx"
#include "View/Reconstructed/ViewRecTrack.hxx"
#if T2S_LEGACY_KF
#include "Legacy/LegacyFitter.hxx"
#include "Legacy/LegacyHelixHelix.hxx"
#else
#include "KalmanFitter/BaseKalmanFitter.hxx"
#include "Seeder/SeederHelixHelix.hxx"
#endif

namespace Tree2Secondaries {

namespace KF = KalmanFitter;

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
    IO::ReadBranches(fInputChain_Events.get(), fInput_Event, IsMC());
    if (IsMC()) {
        IO::ReadBranches(fInputChain_Events.get(), fInput_Injected, false);
        IO::ReadBranches(fInputChain_Events.get(), fInput_MC);
    }
    IO::ReadBranches(fInputChain_Events.get(), fInput_Tracks, IsMC());
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

    IO::CreateBranches(fOutputTree.get(), fOutput_Event, IsMC());
    if (IsMC()) IO::CreateBranches(fOutputTree.get(), fOutput_Injected, true);

    switch (GetReactionChannel()) {
        // standard channels //
        case EReactionChannel::A: {
            IO::CreateBranches(fOutputTree.get(), fOutput_AntiLambdas, Const::V0_Acronym[PID_V0::AntiLambda]);
            IO::CreateBranches(fOutputTree.get(), fOutput_Lambdas, Const::V0_Acronym[PID_V0::Lambda]);
            IO::CreateBranches(fOutputTree.get(), fOutput_KaonsZeroShort, Const::V0_Acronym[PID_V0::KaonZeroShort]);
            if (IsMC()) {
                IO::CreateBranches(fOutputTree.get(), fOutput_MC_AntiLambdas, Const::V0_Acronym[PID_V0::AntiLambda]);
                IO::CreateBranches(fOutputTree.get(), fOutput_MC_Lambdas, Const::V0_Acronym[PID_V0::Lambda]);
                IO::CreateBranches(fOutputTree.get(), fOutput_MC_KaonsZeroShort, Const::V0_Acronym[PID_V0::KaonZeroShort]);
            }
            break;
        }
        case EReactionChannel::D: {
            IO::CreateBranches(fOutputTree.get(), fOutput_AntiLambdas, Const::V0_Acronym[PID_V0::AntiLambda]);
            IO::CreateBranches(fOutputTree.get(), fOutput_Lambdas, Const::V0_Acronym[PID_V0::Lambda]);
            IO::CreateBranches(fOutputTree.get(), fOutput_NegKaons, false, Const::Particle_Acronym[PID_StableParticle::NegKaon]);
            IO::CreateBranches(fOutputTree.get(), fOutput_PosKaons, false, Const::Particle_Acronym[PID_StableParticle::PosKaon]);
            if (IsMC()) {
                IO::CreateBranches(fOutputTree.get(), fOutput_MC_AntiLambdas, Const::V0_Acronym[PID_V0::AntiLambda]);
                IO::CreateBranches(fOutputTree.get(), fOutput_MC_Lambdas, Const::V0_Acronym[PID_V0::Lambda]);
                IO::CreateBranches(fOutputTree.get(), fOutput_MC_NegKaons, true, Const::Particle_Acronym[PID_StableParticle::NegKaon]);
                IO::CreateBranches(fOutputTree.get(), fOutput_MC_PosKaons, true, Const::Particle_Acronym[PID_StableParticle::PosKaon]);
            }
            break;
        }
        case EReactionChannel::E: {
            IO::CreateBranches(fOutputTree.get(), fOutput_AntiLambdas, Const::V0_Acronym[PID_V0::AntiLambda]);
            IO::CreateBranches(fOutputTree.get(), fOutput_Lambdas, Const::V0_Acronym[PID_V0::Lambda]);
            IO::CreateBranches(fOutputTree.get(), fOutput_NegKaons, false, Const::Particle_Acronym[PID_StableParticle::NegKaon]);
            IO::CreateBranches(fOutputTree.get(), fOutput_PosKaons, false, Const::Particle_Acronym[PID_StableParticle::PosKaon]);
            IO::CreateBranches(fOutputTree.get(), fOutput_PiMinus, false, Const::Particle_Acronym[PID_StableParticle::PiMinus]);
            IO::CreateBranches(fOutputTree.get(), fOutput_PiPlus, false, Const::Particle_Acronym[PID_StableParticle::PiPlus]);
            if (IsMC()) {
                IO::CreateBranches(fOutputTree.get(), fOutput_MC_AntiLambdas, Const::V0_Acronym[PID_V0::AntiLambda]);
                IO::CreateBranches(fOutputTree.get(), fOutput_MC_Lambdas, Const::V0_Acronym[PID_V0::Lambda]);
                IO::CreateBranches(fOutputTree.get(), fOutput_MC_NegKaons, true, Const::Particle_Acronym[PID_StableParticle::NegKaon]);
                IO::CreateBranches(fOutputTree.get(), fOutput_MC_PosKaons, true, Const::Particle_Acronym[PID_StableParticle::PosKaon]);
                IO::CreateBranches(fOutputTree.get(), fOutput_MC_PiMinus, true, Const::Particle_Acronym[PID_StableParticle::PiMinus]);
                IO::CreateBranches(fOutputTree.get(), fOutput_MC_PiPlus, true, Const::Particle_Acronym[PID_StableParticle::PiPlus]);
            }
            break;
        }
        case EReactionChannel::H: {
            IO::CreateBranches(fOutputTree.get(), fOutput_NegKaons, false, Const::Particle_Acronym[PID_StableParticle::NegKaon]);
            IO::CreateBranches(fOutputTree.get(), fOutput_PosKaons, false, Const::Particle_Acronym[PID_StableParticle::PosKaon]);
            if (IsMC()) {
                IO::CreateBranches(fOutputTree.get(), fOutput_MC_NegKaons, true, Const::Particle_Acronym[PID_StableParticle::NegKaon]);
                IO::CreateBranches(fOutputTree.get(), fOutput_MC_PosKaons, true, Const::Particle_Acronym[PID_StableParticle::PosKaon]);
            }
            break;
        }
        default: {
            break;
        }
    }  // end of switch statement
}

void Packager::PrepareOutputHistograms() {

    // event counter //
    fHist_EventCounter = std::make_unique<TH1D>("N_Events", ";;N_Events", 1, 0, 1);

    constexpr int x_nbins = 20;
    constexpr float x_min = 0.;
    constexpr float x_max = 20.;
    std::string_view hist_title = ";Cut N;N Passed Cut";

    fHist_CutFlow_AntiProton = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", Const::Particle_Acronym[PID_StableParticle::AntiProton]).c_str(), hist_title.data(), x_nbins, x_min, x_max);
    fHist_CutFlow_Proton = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", Const::Particle_Acronym[PID_StableParticle::Proton]).c_str(), hist_title.data(), x_nbins, x_min, x_max);
    fHist_CutFlow_NegKaon = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", Const::Particle_Acronym[PID_StableParticle::NegKaon]).c_str(), hist_title.data(), x_nbins, x_min, x_max);
    fHist_CutFlow_PosKaon = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", Const::Particle_Acronym[PID_StableParticle::PosKaon]).c_str(), hist_title.data(), x_nbins, x_min, x_max);
    fHist_CutFlow_PiMinus = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", Const::Particle_Acronym[PID_StableParticle::PiMinus]).c_str(), hist_title.data(), x_nbins, x_min, x_max);
    fHist_CutFlow_PiPlus = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", Const::Particle_Acronym[PID_StableParticle::PiPlus]).c_str(), hist_title.data(), x_nbins, x_min, x_max);

    switch (GetReactionChannel()) {
        case EReactionChannel::A:
            fHist_CutFlow_AntiLambda = std::make_unique<TH1D>(  //
                std::format("CutFlow_{}", Const::V0_Acronym[PID_V0::AntiLambda]).c_str(), hist_title.data(), x_nbins, x_min, x_max);
            fHist_CutFlow_Lambda = std::make_unique<TH1D>(  //
                std::format("CutFlow_{}", Const::V0_Acronym[PID_V0::Lambda]).c_str(), hist_title.data(), x_nbins, x_min, x_max);
            fHist_CutFlow_KaonZeroShort = std::make_unique<TH1D>(  //
                std::format("CutFlow_{}", Const::V0_Acronym[PID_V0::KaonZeroShort]).c_str(), hist_title.data(), x_nbins, x_min, x_max);
            break;
        case EReactionChannel::D:
        case EReactionChannel::E:
            fHist_CutFlow_AntiLambda = std::make_unique<TH1D>(  //
                std::format("CutFlow_{}", Const::V0_Acronym[PID_V0::AntiLambda]).c_str(), hist_title.data(), x_nbins, x_min, x_max);
            fHist_CutFlow_Lambda = std::make_unique<TH1D>(  //
                std::format("CutFlow_{}", Const::V0_Acronym[PID_V0::Lambda]).c_str(), hist_title.data(), x_nbins, x_min, x_max);
            break;
        default:
            break;
    }
}

// ## Event ZONE ## //

void Packager::ProcessEvent() {
    fOutput_Event = fInput_Event;
    fHist_EventCounter->Fill(0.);
}

// ## MC/Injected ZONE ## //

// Loop again over all MC particles.
// Select Primary Particles (particles with no mother), generated via the Sexaquark-Reaction Generator, and with valid Reaction ID;
// and store their origin vertex as the coordinates for this particular secondary vertex.
void Packager::Injected_GetSecondaryVertex() {

    fVec_SV_X.resize(NumberInjected(), Const::DummyFloat);
    fVec_SV_Y.resize(NumberInjected(), Const::DummyFloat);
    fVec_SV_Z.resize(NumberInjected(), Const::DummyFloat);

    if (NumberInjected() != Const::NInjectedPerEvent) {
        Logger::Error(__FUNCTION__, "N Sexaquark Reactions is not standard!");
        return;
    }

    std::array<bool, Const::NInjectedPerEvent> sv_found{};

    int reaction_id_lower = Const::ReactionID_Offset;                          // [600,
    int reaction_id_upper = reaction_id_lower + Const::NInjectedPerEvent - 1;  // 619]

    for (std::size_t mc_idx = 0; mc_idx < NumberMC(); ++mc_idx) {

        if ((*fInput_MC.mother_mc_entry)[mc_idx] != Const::DummyInt) continue;
        if ((*fInput_MC.generator)[mc_idx] != Const::SignalGeneratorID) continue;

        int status = (*fInput_MC.status)[mc_idx];
        if (status < reaction_id_lower || status > reaction_id_upper) continue;

        auto reaction_idx = static_cast<std::size_t>(status - Const::ReactionID_Offset);

        if (sv_found[reaction_idx]) continue;

        fVec_SV_X[reaction_idx] = (*fInput_MC.origin.x)[mc_idx];
        fVec_SV_Y[reaction_idx] = (*fInput_MC.origin.y)[mc_idx];
        fVec_SV_Z[reaction_idx] = (*fInput_MC.origin.z)[mc_idx];
        sv_found[reaction_idx] = true;
    }
}

void Packager::Injected_Store() {
    fOutput_Injected = fInput_Injected;
    fOutput_Injected.sv.x = &fVec_SV_X;
    fOutput_Injected.sv.y = &fVec_SV_Y;
    fOutput_Injected.sv.z = &fVec_SV_Z;
}

// ## Tracks ZONE ## //

// Filter and group tracks into indices vectors, according to their respective species.
void Packager::ProcessTracks() {

    std::size_t n_reserve = NumberTracks() / 3;
    fIndices_AntiProtons.reserve(n_reserve);
    fIndices_NegKaons.reserve(n_reserve);
    fIndices_PiMinus.reserve(n_reserve);
    fIndices_Protons.reserve(n_reserve);
    fIndices_PosKaons.reserve(n_reserve);
    fIndices_PiPlus.reserve(n_reserve);

    for (unsigned int esd_idx = 0; esd_idx < NumberTracks(); ++esd_idx) {
        // create track reference //
        View::Rec::Track track{&fInput_Tracks, esd_idx};
        // PID and pre-selection //
        if (track.Charge<int>() < 0) {
            if (Cuts_Proton(track, fHist_CutFlow_AntiProton.get())) {
                fIndices_AntiProtons.emplace_back(esd_idx);
            }
            if (Cuts_Kaon(track, fHist_CutFlow_NegKaon.get())) {
                fIndices_NegKaons.emplace_back(esd_idx);
            }
            if (Cuts_Pion(track, fHist_CutFlow_PiMinus.get())) {
                fIndices_PiMinus.emplace_back(esd_idx);
            }
        }
        if (track.Charge<int>() > 0) {
            if (Cuts_Proton(track, fHist_CutFlow_Proton.get())) {
                fIndices_Protons.emplace_back(esd_idx);
            }
            if (Cuts_Kaon(track, fHist_CutFlow_PosKaon.get())) {
                fIndices_PosKaons.emplace_back(esd_idx);
            }
            if (Cuts_Pion(track, fHist_CutFlow_PiPlus.get())) {
                fIndices_PiPlus.emplace_back(esd_idx);
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

// NOTE: intended for light particles only, i.e., kaons and pions.
void Packager::PackTracks(PID_StableParticle pid) {

    // determine rules based on particle species //

    std::vector<unsigned int>* vec = nullptr;
    Schema::Vector::Tracks* out = nullptr;
    Schema::Vector::MC_Tracks* mc_out = nullptr;

    switch (pid) {
        case PID_StableParticle::NegKaon:
            vec = &fIndices_NegKaons;
            out = &fOutput_NegKaons;
            if (IsMC()) mc_out = &fOutput_MC_NegKaons;
            break;
        case PID_StableParticle::PosKaon:
            vec = &fIndices_PosKaons;
            out = &fOutput_PosKaons;
            if (IsMC()) mc_out = &fOutput_MC_PosKaons;
            break;
        case PID_StableParticle::PiMinus:
            vec = &fIndices_PiMinus;
            out = &fOutput_PiMinus;
            if (IsMC()) mc_out = &fOutput_MC_PiMinus;
            break;
        case PID_StableParticle::PiPlus:
            vec = &fIndices_PiPlus;
            out = &fOutput_PiPlus;
            if (IsMC()) mc_out = &fOutput_MC_PiPlus;
            break;
        default:
            return;
    }

    // loop over selected tracks //
    for (const auto& esd_idx : *vec) {

        // store //
        View::Rec::Track track(&fInput_Tracks, esd_idx);
        Store(track, *out);

        if (IsMC()) {
            View::MC::Particle mc{&fInput_MC, track.McEntry()};
            if (View::IsValid(mc)) {
                StoreMC(mc, *mc_out, pid);
            } else {
                StoreDummyMC(*mc_out);
            }
        }
    }  // end of loop over selected tracks
}

void Packager::Store(const View::Rec::Track& track, Schema::Vector::Tracks& df) {
    df.esd_entry->push_back(track.EsdEntry());
    df.state.x->push_back(track.X());
    df.state.y->push_back(track.Y());
    df.state.z->push_back(track.Z());
    df.state.px->push_back(track.Px());
    df.state.py->push_back(track.Py());
    df.state.pz->push_back(track.Pz());
    df.charge->push_back(track.Charge<char>());
    df.pre_dca_xy->push_back(track.PreDCAxy());
    df.pre_dca_z->push_back(track.PreDCAz());
    df.tpc_signal->push_back(track.TPC_Signal());
    df.n_sigma_pion->push_back(track.NSigmaPion());
    df.n_sigma_kaon->push_back(track.NSigmaKaon());
    df.n_sigma_proton->push_back(track.NSigmaProton());
    track.AppendCov(*df.cov.mat);
}

void Packager::StoreMC(const View::MC::Particle& mc, Schema::Vector::MC_Tracks& df, PID_StableParticle pid) {
    df.origin.x->push_back(mc.Origin_X());
    df.origin.y->push_back(mc.Origin_Y());
    df.origin.z->push_back(mc.Origin_Z());
    df.lv.px->push_back(mc.Px());
    df.lv.py->push_back(mc.Py());
    df.lv.pz->push_back(mc.Pz());
    df.lv.energy->push_back(mc.Energy());

    View::MC::Particle mother{mc.Source, mc.Mother_McEntry()};
    bool valid_mother = View::IsValid(mother);  // used again below for grandmother
    if (valid_mother) {
        df.mother_mc_entry->push_back(mother.Entry);
        df.mother_pdg_code->push_back(mother.PdgCode());
    } else {
        df.mother_mc_entry->push_back(Const::DummyInt);
        df.mother_pdg_code->push_back(Const::DummyInt);
    }

    df.mc_entry->push_back(mc.Entry);
    df.pdg_code->push_back(mc.PdgCode());
    df.reaction_id->push_back(Truth::Track::ReactionID(mc, mother, pid));
    df.is_true->push_back(static_cast<char>(Truth::Track::IsTrue(mc, pid)));
    df.is_signal->push_back(static_cast<char>(Truth::Track::IsSignal(mc, pid)));
    df.is_secondary->push_back(static_cast<char>(Truth::Track::IsSecondary(mc, pid)));

    if (!valid_mother) return;  // early return
    View::MC::Particle grandmother{mother.Source, mother.Mother_McEntry()};
    if (View::IsValid(grandmother)) {
        df.gm_mc_entry->push_back(grandmother.Entry);
        df.gm_pdg_code->push_back(grandmother.PdgCode());
    } else {
        df.gm_mc_entry->push_back(Const::DummyInt);
        df.gm_pdg_code->push_back(Const::DummyInt);
    }
}

void Packager::StoreDummyMC(Schema::Vector::MC_Tracks& df) {

    df.origin.x->push_back(Const::DummyFloat);
    df.origin.y->push_back(Const::DummyFloat);
    df.origin.z->push_back(Const::DummyFloat);
    df.lv.px->push_back(Const::DummyFloat);
    df.lv.py->push_back(Const::DummyFloat);
    df.lv.pz->push_back(Const::DummyFloat);
    df.lv.energy->push_back(Const::DummyFloat);

    df.mother_mc_entry->push_back(Const::DummyInt);
    df.mother_pdg_code->push_back(Const::DummyInt);

    df.mc_entry->push_back(Const::DummyInt);
    df.pdg_code->push_back(Const::DummyInt);
    df.reaction_id->push_back(Const::DummyInt);
    df.is_true->push_back(0);
    df.is_signal->push_back(0);
    df.is_secondary->push_back(0);

    df.gm_mc_entry->push_back(Const::DummyInt);
    df.gm_pdg_code->push_back(Const::DummyInt);
}

bool Packager::Cuts_Proton(const View::Rec::Track& track, TH1D* cut_flow_hist) const {

    cut_flow_hist->Fill(0.);
    if (std::abs(track.NSigmaProton()) > Cuts::Proton::AbsMax_NSigmaProton) return false;
    cut_flow_hist->Fill(1.);
    // if (std::abs(track.DCAxy()) < Cuts::Proton::AbsMin_DCAxy) return false;  // PENDING
    // cut_flow_hist->Fill(2.);                                                 // PENDING

    return true;
}

bool Packager::Cuts_Kaon(const View::Rec::Track& track, TH1D* cut_flow_hist) const {

    cut_flow_hist->Fill(0.);
    if (std::abs(track.NSigmaKaon()) > Cuts::Kaon::AbsMax_NSigmaKaon) return false;
    cut_flow_hist->Fill(1.);
    // if (track.TPCSignal() < Cuts::Kaon::Min_TPCSignal) return false; // PENDING
    cut_flow_hist->Fill(2.);

    return true;
}

bool Packager::Cuts_Pion(const View::Rec::Track& track, TH1D* cut_flow_hist) const {

    cut_flow_hist->Fill(0.);
    if (std::abs(track.NSigmaPion()) > Cuts::Pion::AbsMax_NSigmaPion) return false;
    cut_flow_hist->Fill(1.);
    // if (std::abs(track.DCAxy()) < Cuts::Pion::AbsMin_DCAxy) return false;  // PENDING
    // cut_flow_hist->Fill(2.);                                               // PENDING

    return true;
}

// ## V0s ZONE ## //

void Packager::FindV0s(PID_V0 pid) {

    // determine rules based on V0 species //

    Schema::Vector::V0s* out = nullptr;
    Schema::Vector::MC_V0s* mc_out = nullptr;
    const std::vector<unsigned int>* vec_neg = &fIndices_PiMinus;
    const std::vector<unsigned int>* vec_pos = &fIndices_PiPlus;
    double mass_neg = Const::Particle_Mass[Const::V0_NegativePID[pid]];
    double mass_pos = Const::Particle_Mass[Const::V0_PositivePID[pid]];

    switch (pid) {
        case PID_V0::AntiLambda:
            out = &fOutput_AntiLambdas;
            if (IsMC()) mc_out = &fOutput_MC_AntiLambdas;
            vec_neg = &fIndices_AntiProtons;
            break;
        case PID_V0::Lambda:
            out = &fOutput_Lambdas;
            if (IsMC()) mc_out = &fOutput_MC_Lambdas;
            vec_pos = &fIndices_Protons;
            break;
        case PID_V0::KaonZeroShort:
            out = &fOutput_KaonsZeroShort;
            if (IsMC()) mc_out = &fOutput_MC_KaonsZeroShort;
            break;
        default:
            Logger::Error(__FUNCTION__, "Invalid PID {} for a V0.", Const::V0_Acronym[pid]);
            return;
    }

    // loop over all possible pairs of tracks //
    for (const auto& esd_neg : *vec_neg) {
        for (const auto& esd_pos : *vec_pos) {

            // sanity check //
            if (esd_neg == esd_pos) continue;

#if T2S_DEBUG
            if (pid == PID_V0::AntiLambda && (esd_neg != 725 || esd_pos != 739)) continue;
            if (pid == PID_V0::Lambda) continue;
            if (pid == PID_V0::KaonZeroShort && (esd_neg != 129 || esd_pos != 425)) continue;
#endif

            // get views //
            View::Rec::Track neg(&fInput_Tracks, esd_neg);
            View::Rec::Track pos(&fInput_Tracks, esd_pos);

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
            auto fit = KF::FitVertex(neg, pos, mass_neg, mass_pos, {seed_neg, deriv_neg}, {seed_pos, deriv_pos}, fInput_Event.magnetic_field);
#endif

            // create composite particle //
            KF::V0 v0(fit, seed_neg.pca, seed_pos.pca, neg, pos);

#if T2S_DEBUG
            Logger::Debug(__FUNCTION__, "{}", fit);
            Logger::Debug(__FUNCTION__, "  neg,pos={},{}", esd_neg, esd_pos);
            Logger::Debug(__FUNCTION__, "  x,y,z(neg)={},{},{}", v0.Neg_PCAwrtV0.X(), v0.Neg_PCAwrtV0.Y(), v0.Neg_PCAwrtV0.Z());
            Logger::Debug(__FUNCTION__, "  x,y,z(pos)={},{},{}", v0.Pos_PCAwrtV0.X(), v0.Pos_PCAwrtV0.Y(), v0.Pos_PCAwrtV0.Z());
#if T2S_LEGACY_KF == OFF
            Logger::Debug(__FUNCTION__, "  dca_dau={}", v0.DCA_Daughters());
            Logger::Debug(__FUNCTION__, "  dca_neg={}", v0.DCA_Neg_V0());
            Logger::Debug(__FUNCTION__, "  dca_pos={}", v0.DCA_Pos_V0());
            Logger::Debug(__FUNCTION__, "  qt={}", v0.ArmenterosQt());
            Logger::Debug(__FUNCTION__, "  alpha={}", v0.ArmenterosAlpha());
            Logger::Debug(__FUNCTION__, "  cpa_pv={}", v0.CPA_Vertex(fInput_Event.pv.x, fInput_Event.pv.y, fInput_Event.pv.z));
            Logger::Debug(__FUNCTION__, "  dca_pv={}", v0.DCA_Vertex(fInput_Event.pv.x, fInput_Event.pv.y, fInput_Event.pv.z));
#endif
#else
            // apply cuts (2) //
            if (!SlowCuts(v0, pid)) continue;
#endif

            // store //
            Store(v0, *out);

            // store mc //
            if (IsMC()) {
                View::MC::Particle mc_neg(&fInput_MC, neg.McEntry());
                View::MC::Particle mc_pos(&fInput_MC, pos.McEntry());
                View::MC::V0 mc_v0(mc_neg, mc_pos);
                if (View::IsValid(mc_v0)) {
                    StoreMC(mc_v0, *mc_out, pid);
                } else {
                    StoreDummyMC(*mc_out);
                }
            }
        }  // end of loop over pos
    }  // end of loop over neg
}

bool Packager::FastCuts_Lambda(const Seeder::Seed& seed_neg, const Seeder::Seed& seed_pos, TH1D* cut_flow_hist) const {
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

bool Packager::SlowCuts_Lambda(const KF::V0& v0, TH1D* cut_flow_hist) const {

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

void Packager::Store(const KF::V0& v0, Schema::Vector::V0s& df) {
    df.decay.x->push_back(static_cast<float>(v0.X()));
    df.decay.y->push_back(static_cast<float>(v0.Y()));
    df.decay.z->push_back(static_cast<float>(v0.Z()));
    df.lv.px->push_back(static_cast<float>(v0.Px()));
    df.lv.py->push_back(static_cast<float>(v0.Py()));
    df.lv.pz->push_back(static_cast<float>(v0.Pz()));
    df.lv.energy->push_back(static_cast<float>(v0.E()));
    v0.AppendCov<7>(*df.cov.mat);

    // -- fit info
    df.chi2ndf->push_back(static_cast<float>(v0.Chi2NDF()));

    // -- neg daughter
    df.neg.esd_entry->push_back(v0.Neg.EsdEntry());
    df.neg.state.x->push_back(v0.Neg.X());
    df.neg.state.y->push_back(v0.Neg.Y());
    df.neg.state.z->push_back(v0.Neg.Z());
    df.neg.state.px->push_back(v0.Neg.Px());
    df.neg.state.py->push_back(v0.Neg.Py());
    df.neg.state.pz->push_back(v0.Neg.Pz());
    v0.Neg.AppendCov(*df.neg.cov.mat);
    df.neg.charge->push_back(v0.Neg.Charge<char>());
    df.neg.pre_dca_xy->push_back(v0.Neg.PreDCAxy());
    df.neg.pre_dca_z->push_back(v0.Neg.PreDCAz());
    df.neg.tpc_signal->push_back(v0.Neg.TPC_Signal());
    df.neg.n_sigma_pion->push_back(v0.Neg.NSigmaPion());
    df.neg.n_sigma_kaon->push_back(v0.Neg.NSigmaKaon());
    df.neg.n_sigma_proton->push_back(v0.Neg.NSigmaProton());

    df.neg_pca_v0.x->push_back(static_cast<float>(v0.Neg_PCAwrtV0.X()));
    df.neg_pca_v0.y->push_back(static_cast<float>(v0.Neg_PCAwrtV0.Y()));
    df.neg_pca_v0.z->push_back(static_cast<float>(v0.Neg_PCAwrtV0.Z()));
    df.neg_pca_v0.px->push_back(static_cast<float>(v0.Neg_PCAwrtV0.Px()));
    df.neg_pca_v0.py->push_back(static_cast<float>(v0.Neg_PCAwrtV0.Py()));
    df.neg_pca_v0.pz->push_back(static_cast<float>(v0.Neg_PCAwrtV0.Pz()));

    // -- pos daughter
    df.pos.esd_entry->push_back(v0.Pos.EsdEntry());
    df.pos.state.x->push_back(v0.Pos.X());
    df.pos.state.y->push_back(v0.Pos.Y());
    df.pos.state.z->push_back(v0.Pos.Z());
    df.pos.state.px->push_back(v0.Pos.Px());
    df.pos.state.py->push_back(v0.Pos.Py());
    df.pos.state.pz->push_back(v0.Pos.Pz());
    v0.Pos.AppendCov(*df.pos.cov.mat);
    df.pos.charge->push_back(v0.Pos.Charge<char>());
    df.pos.pre_dca_xy->push_back(v0.Pos.PreDCAxy());
    df.pos.pre_dca_z->push_back(v0.Pos.PreDCAz());
    df.pos.tpc_signal->push_back(v0.Pos.TPC_Signal());
    df.pos.n_sigma_pion->push_back(v0.Pos.NSigmaPion());
    df.pos.n_sigma_kaon->push_back(v0.Pos.NSigmaKaon());
    df.pos.n_sigma_proton->push_back(v0.Pos.NSigmaProton());

    df.pos_pca_v0.x->push_back(static_cast<float>(v0.Pos_PCAwrtV0.X()));
    df.pos_pca_v0.y->push_back(static_cast<float>(v0.Pos_PCAwrtV0.Y()));
    df.pos_pca_v0.z->push_back(static_cast<float>(v0.Pos_PCAwrtV0.Z()));
    df.pos_pca_v0.px->push_back(static_cast<float>(v0.Pos_PCAwrtV0.Px()));
    df.pos_pca_v0.py->push_back(static_cast<float>(v0.Pos_PCAwrtV0.Py()));
    df.pos_pca_v0.pz->push_back(static_cast<float>(v0.Pos_PCAwrtV0.Pz()));
}

void Packager::StoreMC(const View::MC::V0& mc_v0, Schema::Vector::MC_V0s& df, PID_V0 pid) {

    df.origin.x->push_back(mc_v0.Origin_X());
    df.origin.y->push_back(mc_v0.Origin_Y());
    df.origin.z->push_back(mc_v0.Origin_Z());

    df.lv.px->push_back(mc_v0.Px());
    df.lv.py->push_back(mc_v0.Py());
    df.lv.pz->push_back(mc_v0.Pz());
    df.lv.energy->push_back(mc_v0.Energy());

    View::MC::Particle mother{mc_v0.Source, mc_v0.Mother_McEntry()};
    if (View::IsValid(mother)) {
        df.mother_mc_entry->push_back(mother.Entry);
        df.mother_pdg_code->push_back(mother.PdgCode());
    } else {
        df.mother_mc_entry->push_back(Const::DummyInt);
        df.mother_pdg_code->push_back(Const::DummyInt);
    }

    df.mc_entry->push_back(mc_v0.Entry);
    df.pdg_code->push_back(mc_v0.PdgCode());
    df.reaction_id->push_back(Truth::V0::ReactionID(mc_v0, pid));
    df.is_true->push_back(static_cast<char>(Truth::V0::IsTrue(mc_v0, pid)));
    df.is_signal->push_back(static_cast<char>(Truth::V0::IsSignal(mc_v0, pid)));
    df.is_secondary->push_back(static_cast<char>(Truth::V0::IsSecondary(mc_v0, pid)));
    df.is_hybrid->push_back(static_cast<char>(Truth::V0::IsHybrid(mc_v0, pid)));

    df.decay.x->push_back(Truth::V0::Decay_X(mc_v0));
    df.decay.y->push_back(Truth::V0::Decay_Y(mc_v0));
    df.decay.z->push_back(Truth::V0::Decay_Z(mc_v0));

    df.neg.lv.px->push_back(mc_v0.Neg.Px());
    df.neg.lv.py->push_back(mc_v0.Neg.Py());
    df.neg.lv.pz->push_back(mc_v0.Neg.Pz());
    df.neg.lv.energy->push_back(mc_v0.Neg.Energy());

    df.neg.mc_entry->push_back(mc_v0.Neg.Entry);
    df.neg.pdg_code->push_back(mc_v0.Neg.PdgCode());
    df.neg.reaction_id->push_back(Truth::Track::ReactionID(mc_v0.Neg, mc_v0, Const::V0_NegativePID[pid]));
    df.neg.is_true->push_back(static_cast<char>(Truth::Track::IsTrue(mc_v0.Neg, Const::V0_NegativePID[pid])));
    df.neg.is_signal->push_back(static_cast<char>(Truth::Track::IsSignal(mc_v0.Neg, Const::V0_NegativePID[pid])));
    df.neg.is_secondary->push_back(static_cast<char>(Truth::Track::IsSecondary(mc_v0.Neg, Const::V0_NegativePID[pid])));

    df.pos.lv.px->push_back(mc_v0.Pos.Px());
    df.pos.lv.py->push_back(mc_v0.Pos.Py());
    df.pos.lv.pz->push_back(mc_v0.Pos.Pz());
    df.pos.lv.energy->push_back(mc_v0.Pos.Energy());

    df.pos.mc_entry->push_back(mc_v0.Pos.Entry);
    df.pos.pdg_code->push_back(mc_v0.Pos.PdgCode());
    df.pos.reaction_id->push_back(Truth::Track::ReactionID(mc_v0.Pos, mc_v0, Const::V0_PositivePID[pid]));
    df.pos.is_true->push_back(static_cast<char>(Truth::Track::IsTrue(mc_v0.Pos, Const::V0_PositivePID[pid])));
    df.pos.is_signal->push_back(static_cast<char>(Truth::Track::IsSignal(mc_v0.Pos, Const::V0_PositivePID[pid])));
    df.pos.is_secondary->push_back(static_cast<char>(Truth::Track::IsSecondary(mc_v0.Pos, Const::V0_PositivePID[pid])));
}

void Packager::StoreDummyMC(Schema::Vector::MC_V0s& df) {

    df.origin.x->push_back(Const::DummyFloat);
    df.origin.y->push_back(Const::DummyFloat);
    df.origin.z->push_back(Const::DummyFloat);
    df.lv.px->push_back(Const::DummyFloat);
    df.lv.py->push_back(Const::DummyFloat);
    df.lv.pz->push_back(Const::DummyFloat);
    df.lv.energy->push_back(Const::DummyFloat);

    df.mother_mc_entry->push_back(Const::DummyInt);
    df.mother_pdg_code->push_back(Const::DummyInt);

    df.mc_entry->push_back(Const::DummyInt);
    df.pdg_code->push_back(Const::DummyInt);
    df.reaction_id->push_back(Const::DummyInt);
    df.is_true->push_back(0);
    df.is_signal->push_back(0);
    df.is_secondary->push_back(0);

    df.decay.x->push_back(Const::DummyFloat);
    df.decay.y->push_back(Const::DummyFloat);
    df.decay.z->push_back(Const::DummyFloat);

    df.is_hybrid->push_back(0);

    df.neg.lv.px->push_back(Const::DummyFloat);
    df.neg.lv.py->push_back(Const::DummyFloat);
    df.neg.lv.pz->push_back(Const::DummyFloat);
    df.neg.lv.energy->push_back(Const::DummyFloat);

    df.neg.mc_entry->push_back(Const::DummyInt);
    df.neg.pdg_code->push_back(Const::DummyInt);
    df.neg.reaction_id->push_back(Const::DummyInt);
    df.neg.is_true->push_back(0);
    df.neg.is_signal->push_back(0);
    df.neg.is_secondary->push_back(0);

    df.pos.lv.px->push_back(Const::DummyFloat);
    df.pos.lv.py->push_back(Const::DummyFloat);
    df.pos.lv.pz->push_back(Const::DummyFloat);
    df.pos.lv.energy->push_back(Const::DummyFloat);

    df.pos.mc_entry->push_back(Const::DummyInt);
    df.pos.pdg_code->push_back(Const::DummyInt);
    df.pos.reaction_id->push_back(Const::DummyInt);
    df.pos.is_true->push_back(0);
    df.pos.is_signal->push_back(0);
    df.pos.is_secondary->push_back(0);
}

// ## END OF CYCLES ## //

void Packager::EndOfEvent() {

    // fill tree//

    fOutputTree->Fill();

    // clear temporary containers //

    fVec_SV_X.clear();
    fVec_SV_Y.clear();
    fVec_SV_Z.clear();

    fIndices_AntiProtons.clear();
    fIndices_Protons.clear();
    fIndices_NegKaons.clear();
    fIndices_PosKaons.clear();
    fIndices_PiMinus.clear();
    fIndices_PiPlus.clear();

    // clear output vector branches //

    if (IsMC()) IO::ClearBranches(fOutput_Injected, true);

    switch (GetReactionChannel()) {
        case EReactionChannel::A: {
            IO::ClearBranches(fOutput_AntiLambdas);
            IO::ClearBranches(fOutput_Lambdas);
            IO::ClearBranches(fOutput_KaonsZeroShort);
            if (IsMC()) {
                IO::ClearBranches(fOutput_MC_AntiLambdas);
                IO::ClearBranches(fOutput_MC_Lambdas);
                IO::ClearBranches(fOutput_MC_KaonsZeroShort);
            }
            break;
        }
        case EReactionChannel::D: {
            IO::ClearBranches(fOutput_AntiLambdas);
            IO::ClearBranches(fOutput_Lambdas);
            IO::ClearBranches(fOutput_NegKaons, false);
            IO::ClearBranches(fOutput_PosKaons, false);
            if (IsMC()) {
                IO::ClearBranches(fOutput_MC_AntiLambdas);
                IO::ClearBranches(fOutput_MC_Lambdas);
                IO::ClearBranches(fOutput_MC_NegKaons, true);
                IO::ClearBranches(fOutput_MC_PosKaons, true);
            }
            break;
        }
        case EReactionChannel::E: {
            IO::ClearBranches(fOutput_AntiLambdas);
            IO::ClearBranches(fOutput_Lambdas);
            IO::ClearBranches(fOutput_NegKaons, false);
            IO::ClearBranches(fOutput_PosKaons, false);
            IO::ClearBranches(fOutput_PiMinus, false);
            IO::ClearBranches(fOutput_PiPlus, false);
            if (IsMC()) {
                IO::ClearBranches(fOutput_MC_AntiLambdas);
                IO::ClearBranches(fOutput_MC_Lambdas);
                IO::ClearBranches(fOutput_MC_NegKaons, true);
                IO::ClearBranches(fOutput_MC_PosKaons, true);
                IO::ClearBranches(fOutput_MC_PiMinus, true);
                IO::ClearBranches(fOutput_MC_PiPlus, true);
            }
            break;
        }
        case EReactionChannel::H: {
            IO::ClearBranches(fOutput_NegKaons, false);
            IO::ClearBranches(fOutput_PosKaons, false);
            if (IsMC()) {
                IO::ClearBranches(fOutput_MC_NegKaons, true);
                IO::ClearBranches(fOutput_MC_PosKaons, true);
            }
            break;
        }
        default: {
            break;
        }
    }
}

void Packager::EndOfAnalysis() {

    Logger::Info(__FUNCTION__, "The following objects have been written into TFile {}:", fSettings.PathOutputFile);

    // write tree
    fOutputTree->Write();
    Logger::Info(__FUNCTION__, "- TTree \"{}\"", fOutputTree->GetName());

    // write histograms
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
    switch (GetReactionChannel()) {
        case EReactionChannel::A:
            fHist_CutFlow_AntiLambda->Write();
            Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow_AntiLambda->GetName());
            fHist_CutFlow_Lambda->Write();
            Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow_Lambda->GetName());
            fHist_CutFlow_KaonZeroShort->Write();
            Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow_KaonZeroShort->GetName());
            break;
        case EReactionChannel::D:
        case EReactionChannel::E:
            fHist_CutFlow_AntiLambda->Write();
            Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow_AntiLambda->GetName());
            fHist_CutFlow_Lambda->Write();
            Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow_Lambda->GetName());
            break;
        default:
            break;
    }

    fInputChain_Events->ResetBranchAddresses();
    fOutputTree->ResetBranchAddresses();

    Logger::Info(__FUNCTION__, "All done.");
}

}  // namespace Tree2Secondaries
