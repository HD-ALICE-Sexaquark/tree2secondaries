#include <filesystem>

#include "App/Logger.hxx"
#include "Fit/FitTrack.hxx"
#include "Math/Constants.hxx"
#include "Packager/Packager.hxx"
#include "Packager/PackagerCuts.hxx"
#include "Truth/TruthParticle.hxx"
#include "View/MC/ViewMcParticle.hxx"
#include "View/Reconstructed/ViewTrack.hxx"

namespace Tree2Secondaries {

bool Packager::Initialize() {

    fInputChain_Events = std::make_unique<TChain>(Const::TreeName_Events.c_str());
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
    fInput_Event.ReadBranches_FlatEvent(fInputChain_Events.get(), IsMC());
    if (IsMC()) {
        fInput_MC_PV.ReadBranches_FlatCoordinates(fInputChain_Events.get(), "MC_PV", "v");
        fInput_MC.ReadBranches_VectorMCParticles(fInputChain_Events.get());
        fInput_Injected.ReadBranches_VectorInjected(fInputChain_Events.get(), false);
    }
    fInput_Tracks.ReadBranches_VectorTracks(fInputChain_Events.get(), IsMC(), false, "Track");
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

    fOutputTree = std::make_unique<TTree>(Const::TreeName_PackedEvents.c_str(), "Packed Events");
    if (!fOutputTree) {
        Logger::Error(__FUNCTION__, "Couldn't create TTree \"{}\"", fOutputTree->GetName());
        return false;
    }

    return true;
}

void Packager::CreateOutputBranches() {

    fOutput_Event.CreateBranches_FlatEvent(fOutputTree.get(), IsMC());
    if (IsMC()) {
        fOutput_MC_PV.CreateBranches_FlatCoordinates(fOutputTree.get(), "MC_PV", "v");
        fOutput_Injected.CreateBranches_VectorInjected(fOutputTree.get(), true);
    }

    switch (GetReactionChannel()) {
        // standard channels //
        case EReactionChannel::A:
            fOutput_AntiLambdas.CreateBranches_VectorV0s(fOutputTree.get(), Const::Particle_Acronym[EParticle::AntiLambda]);
            fOutput_Lambdas.CreateBranches_VectorV0s(fOutputTree.get(), Const::Particle_Acronym[EParticle::Lambda]);
            fOutput_KaonsZeroShort.CreateBranches_VectorV0s(fOutputTree.get(), Const::Particle_Acronym[EParticle::KaonZeroShort]);
            if (IsMC()) {
                fOutput_MC_AntiLambdas.CreateBranches_VectorMC_V0s(fOutputTree.get(), Const::Particle_Acronym[EParticle::AntiLambda]);
                fOutput_MC_Lambdas.CreateBranches_VectorMC_V0s(fOutputTree.get(), Const::Particle_Acronym[EParticle::Lambda]);
                fOutput_MC_KaonsZeroShort.CreateBranches_VectorMC_V0s(fOutputTree.get(), Const::Particle_Acronym[EParticle::KaonZeroShort]);
            }
            break;
        case EReactionChannel::D:
            fOutput_AntiLambdas.CreateBranches_VectorV0s(fOutputTree.get(), Const::Particle_Acronym[EParticle::AntiLambda]);
            fOutput_Lambdas.CreateBranches_VectorV0s(fOutputTree.get(), Const::Particle_Acronym[EParticle::Lambda]);
            fOutput_NegKaons.CreateBranches_VectorTracks(fOutputTree.get(), false, true, Const::Particle_Acronym[NegKaon]);
            fOutput_PosKaons.CreateBranches_VectorTracks(fOutputTree.get(), false, true, Const::Particle_Acronym[PosKaon]);
            if (IsMC()) {
                fOutput_MC_AntiLambdas.CreateBranches_VectorMC_V0s(fOutputTree.get(), Const::Particle_Acronym[EParticle::AntiLambda]);
                fOutput_MC_Lambdas.CreateBranches_VectorMC_V0s(fOutputTree.get(), Const::Particle_Acronym[EParticle::Lambda]);
                fOutput_MC_NegKaons.CreateBranches_VectorMC_Tracks(fOutputTree.get(), Const::Particle_Acronym[EParticle::NegKaon]);
                fOutput_MC_PosKaons.CreateBranches_VectorMC_Tracks(fOutputTree.get(), Const::Particle_Acronym[EParticle::PosKaon]);
            }
            break;
        case EReactionChannel::E:
            fOutput_AntiLambdas.CreateBranches_VectorV0s(fOutputTree.get(), Const::Particle_Acronym[EParticle::AntiLambda]);
            fOutput_Lambdas.CreateBranches_VectorV0s(fOutputTree.get(), Const::Particle_Acronym[EParticle::Lambda]);
            fOutput_NegKaons.CreateBranches_VectorTracks(fOutputTree.get(), false, true, Const::Particle_Acronym[NegKaon]);
            fOutput_PosKaons.CreateBranches_VectorTracks(fOutputTree.get(), false, true, Const::Particle_Acronym[PosKaon]);
            fOutput_PiMinus.CreateBranches_VectorTracks(fOutputTree.get(), false, true, Const::Particle_Acronym[PiMinus]);
            fOutput_PiPlus.CreateBranches_VectorTracks(fOutputTree.get(), false, true, Const::Particle_Acronym[PiPlus]);
            if (IsMC()) {
                fOutput_MC_AntiLambdas.CreateBranches_VectorMC_V0s(fOutputTree.get(), Const::Particle_Acronym[EParticle::AntiLambda]);
                fOutput_MC_Lambdas.CreateBranches_VectorMC_V0s(fOutputTree.get(), Const::Particle_Acronym[EParticle::Lambda]);
                fOutput_MC_NegKaons.CreateBranches_VectorMC_Tracks(fOutputTree.get(), Const::Particle_Acronym[EParticle::NegKaon]);
                fOutput_MC_PosKaons.CreateBranches_VectorMC_Tracks(fOutputTree.get(), Const::Particle_Acronym[EParticle::PosKaon]);
                fOutput_MC_PiMinus.CreateBranches_VectorMC_Tracks(fOutputTree.get(), Const::Particle_Acronym[EParticle::PiMinus]);
                fOutput_MC_PiPlus.CreateBranches_VectorMC_Tracks(fOutputTree.get(), Const::Particle_Acronym[EParticle::PiPlus]);
            }
            break;
        case EReactionChannel::H:
            fOutput_NegKaons.CreateBranches_VectorTracks(fOutputTree.get(), false, true, Const::Particle_Acronym[NegKaon]);
            fOutput_PosKaons.CreateBranches_VectorTracks(fOutputTree.get(), false, true, Const::Particle_Acronym[PosKaon]);
            if (IsMC()) {
                fOutput_MC_NegKaons.CreateBranches_VectorMC_Tracks(fOutputTree.get(), Const::Particle_Acronym[EParticle::NegKaon]);
                fOutput_MC_PosKaons.CreateBranches_VectorMC_Tracks(fOutputTree.get(), Const::Particle_Acronym[EParticle::PosKaon]);
            }
            break;
    }  // end of switch statement
}

void Packager::PrepareOutputHistograms() {

    // event counter //
    fHist_EventCounter = std::make_unique<TH1D>("N_Events", ";;N_Events", 1, 0, 1);

    const int x_nbins{20};
    const float x_min{0.};
    const float x_max{20.};
    std::string hist_title{";Cut N;N Passed Cut"};

    fHist_CutFlow_AntiProton = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", Const::Particle_Acronym[EParticle::AntiProton]).c_str(), hist_title.c_str(), x_nbins, x_min, x_max);
    fHist_CutFlow_Proton = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", Const::Particle_Acronym[EParticle::Proton]).c_str(), hist_title.c_str(), x_nbins, x_min, x_max);
    fHist_CutFlow_NegKaon = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", Const::Particle_Acronym[EParticle::NegKaon]).c_str(), hist_title.c_str(), x_nbins, x_min, x_max);
    fHist_CutFlow_PosKaon = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", Const::Particle_Acronym[EParticle::PosKaon]).c_str(), hist_title.c_str(), x_nbins, x_min, x_max);
    fHist_CutFlow_PiMinus = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", Const::Particle_Acronym[EParticle::PiMinus]).c_str(), hist_title.c_str(), x_nbins, x_min, x_max);
    fHist_CutFlow_PiPlus = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", Const::Particle_Acronym[EParticle::PiPlus]).c_str(), hist_title.c_str(), x_nbins, x_min, x_max);

    switch (GetReactionChannel()) {
        case EReactionChannel::A:
            fHist_CutFlow_AntiLambda = std::make_unique<TH1D>(  //
                std::format("CutFlow_{}", Const::Particle_Acronym[EParticle::AntiLambda]).c_str(), hist_title.c_str(), x_nbins, x_min, x_max);
            fHist_CutFlow_Lambda = std::make_unique<TH1D>(  //
                std::format("CutFlow_{}", Const::Particle_Acronym[EParticle::Lambda]).c_str(), hist_title.c_str(), x_nbins, x_min, x_max);
            fHist_CutFlow_KaonZeroShort = std::make_unique<TH1D>(  //
                std::format("CutFlow_{}", Const::Particle_Acronym[EParticle::KaonZeroShort]).c_str(), hist_title.c_str(), x_nbins, x_min, x_max);
            break;
        case EReactionChannel::D:
        case EReactionChannel::E:
            fHist_CutFlow_AntiLambda = std::make_unique<TH1D>(  //
                std::format("CutFlow_{}", Const::Particle_Acronym[EParticle::AntiLambda]).c_str(), hist_title.c_str(), x_nbins, x_min, x_max);
            fHist_CutFlow_Lambda = std::make_unique<TH1D>(  //
                std::format("CutFlow_{}", Const::Particle_Acronym[EParticle::Lambda]).c_str(), hist_title.c_str(), x_nbins, x_min, x_max);
            break;
        default:
            break;
    }
}

// ## Event ZONE ## //

void Packager::ProcessEvent() {
    fOutput_Event = fInput_Event;
    if (IsMC()) fOutput_MC_PV = fInput_MC_PV;
    fHist_EventCounter->Fill(0.);
}

// ## MC/Injected ZONE ## //

void Packager::Injected_GetSecondaryVertex() {

    fVec_SV_X.resize(NumberInjected(), 0.);
    fVec_SV_Y.resize(NumberInjected(), 0.);
    fVec_SV_Z.resize(NumberInjected(), 0.);

    for (int idx_mc{0}; idx_mc < NumberMC(); ++idx_mc) {
        if (fInput_MC.MotherEntry->at(idx_mc) != -1) continue;
        if (fInput_MC.Generator->at(idx_mc) != 2) continue;

        int status{fInput_MC.Status->at(idx_mc)};
        if (status < 600 || status > 619) continue;

        if (fVec_SV_X[status - Const::ReactionID_Offset] != 0.) continue;
        fVec_SV_X[status - Const::ReactionID_Offset] = fInput_MC.X->at(idx_mc);
        fVec_SV_Y[status - Const::ReactionID_Offset] = fInput_MC.Y->at(idx_mc);
        fVec_SV_Z[status - Const::ReactionID_Offset] = fInput_MC.Z->at(idx_mc);
    }
}

void Packager::Injected_Store() {
    fOutput_Injected.ReactionID = fInput_Injected.ReactionID;
    fOutput_Injected.SV.X = &fVec_SV_X;
    fOutput_Injected.SV.Y = &fVec_SV_Y;
    fOutput_Injected.SV.Z = &fVec_SV_Z;
    fOutput_Injected.Px = fInput_Injected.Px;
    fOutput_Injected.Py = fInput_Injected.Py;
    fOutput_Injected.Pz = fInput_Injected.Pz;
    fOutput_Injected.Nucleon.Px = fInput_Injected.Nucleon.Px;
    fOutput_Injected.Nucleon.Py = fInput_Injected.Nucleon.Py;
    fOutput_Injected.Nucleon.Pz = fInput_Injected.Nucleon.Pz;
}

// ## Tracks ZONE ## //

// Filter and group tracks into indices vectors, according to their respective species.
void Packager::ProcessTracks() {

    fVec_AntiProtons.reserve(NumberTracks());
    fVec_NegKaons.reserve(NumberTracks());
    fVec_PiMinus.reserve(NumberTracks());
    fVec_Protons.reserve(NumberTracks());
    fVec_PosKaons.reserve(NumberTracks());
    fVec_PiPlus.reserve(NumberTracks());

    for (int esd_idx{0}; esd_idx < NumberTracks(); ++esd_idx) {
        // create track reference //
        View::Rec::Track track{&fInput_Tracks, esd_idx};
        // PID and pre-selection //
        if (track.Charge() < 0) {
            if (PassesCuts_Proton(track, fHist_CutFlow_AntiProton.get())) {
                fVec_AntiProtons.emplace_back(track, Const::Particle_Mass[EParticle::AntiProton]);
            }
            if (PassesCuts_Kaon(track, fHist_CutFlow_NegKaon.get())) {
                fVec_NegKaons.emplace_back(track, Const::Particle_Mass[EParticle::NegKaon]);
            }
            if (PassesCuts_Pion(track, fHist_CutFlow_PiMinus.get())) {
                fVec_PiMinus.emplace_back(track, Const::Particle_Mass[EParticle::PiMinus]);
            }
        }
        if (track.Charge() > 0) {
            if (PassesCuts_Proton(track, fHist_CutFlow_Proton.get())) {
                fVec_Protons.emplace_back(track, Const::Particle_Mass[EParticle::Proton]);
            }
            if (PassesCuts_Kaon(track, fHist_CutFlow_PosKaon.get())) {
                fVec_PosKaons.emplace_back(track, Const::Particle_Mass[EParticle::PosKaon]);
            }
            if (PassesCuts_Pion(track, fHist_CutFlow_PiPlus.get())) {
                fVec_PiPlus.emplace_back(track, Const::Particle_Mass[EParticle::PiPlus]);
            }
        }
    }  // end of loop over tracks

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "n_antiprotons = {}", fVec_AntiProtons.size());
    Logger::Debug(__FUNCTION__, "n_protons     = {}", fVec_Protons.size());
    Logger::Debug(__FUNCTION__, "n_negkaons    = {}", fVec_NegKaons.size());
    Logger::Debug(__FUNCTION__, "n_poskaons    = {}", fVec_PosKaons.size());
    Logger::Debug(__FUNCTION__, "n_piminus     = {}", fVec_PiMinus.size());
    Logger::Debug(__FUNCTION__, "n_piplus      = {}", fVec_PiPlus.size());
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
}

// NOTE: intended for light particles only, i.e., kaons and pions.
void Packager::PackTracks(EParticle pid) {

    // determine rules based on particle species
    std::vector<Fit::Track>* vec{nullptr};
    Storage::Vector::Tracks* out{nullptr};  // PENDING
    Storage::Vector::MC_Tracks* mc_out{nullptr};
    switch (pid) {
        case EParticle::NegKaon:
            vec = &fVec_NegKaons;
            out = &fOutput_NegKaons;
            if (IsMC()) mc_out = &fOutput_MC_NegKaons;
            break;
        case EParticle::PosKaon:
            vec = &fVec_PosKaons;
            out = &fOutput_PosKaons;
            if (IsMC()) mc_out = &fOutput_MC_PosKaons;
            break;
        case EParticle::PiMinus:
            vec = &fVec_PiMinus;
            out = &fOutput_PiMinus;
            if (IsMC()) mc_out = &fOutput_MC_PiMinus;
            break;
        case EParticle::PiPlus:
            vec = &fVec_PiPlus;
            out = &fOutput_PiPlus;
            if (IsMC()) mc_out = &fOutput_MC_PiPlus;
            break;
        default:
            return;
    }

    // loop over selected tracks
    for (const auto& kf_track : *vec) {

        // store
        Store(kf_track.View, *out);

        if (IsMC()) {
            View::MC::Particle mc{&fInput_MC, fInput_Tracks.McEntry->at(kf_track.View.Entry)};
            StoreMC(mc, *mc_out, pid);
        }
    }  // end of loop over selected tracks
}

void Packager::Store(const View::Rec::Track& view, Storage::Vector::Tracks& df) {
    // `Vector::Tracks`
    // -- `Vector::States_NoE`
    df.X->push_back(view.X());
    df.Y->push_back(view.Y());
    df.Z->push_back(view.Z());
    df.Px->push_back(view.Px());
    df.Py->push_back(view.Py());
    df.Pz->push_back(view.Pz());
    // -- `Vector::CovMatrices_NoE`
    df.SigmaX2->push_back(view.SigmaX2());
    df.SigmaXY->push_back(view.SigmaXY());
    df.SigmaY2->push_back(view.SigmaY2());
    df.SigmaXZ->push_back(view.SigmaXZ());
    df.SigmaYZ->push_back(view.SigmaYZ());
    df.SigmaZ2->push_back(view.SigmaZ2());
    df.SigmaXPx->push_back(view.SigmaXPx());
    df.SigmaYPx->push_back(view.SigmaYPx());
    df.SigmaZPx->push_back(view.SigmaZPx());
    df.SigmaPx2->push_back(view.SigmaPx2());
    df.SigmaXPy->push_back(view.SigmaXPy());
    df.SigmaYPy->push_back(view.SigmaYPy());
    df.SigmaZPy->push_back(view.SigmaZPy());
    df.SigmaPxPy->push_back(view.SigmaPxPy());
    df.SigmaPy2->push_back(view.SigmaPy2());
    df.SigmaXPz->push_back(view.SigmaXPz());
    df.SigmaYPz->push_back(view.SigmaYPz());
    df.SigmaZPz->push_back(view.SigmaZPz());
    df.SigmaPxPz->push_back(view.SigmaPxPz());
    df.SigmaPyPz->push_back(view.SigmaPyPz());
    df.SigmaPz2->push_back(view.SigmaPz2());
    // -- rest of variables
    df.Charge->push_back(view.Charge());
    df.DCAxy->push_back(view.DCAxy());
    df.DCAz->push_back(view.DCAz());
    df.TPCSignal->push_back(view.TPCSignal());
    df.NSigmaPion->push_back(view.NSigmaPion());
    df.NSigmaKaon->push_back(view.NSigmaKaon());
    df.NSigmaProton->push_back(view.NSigmaProton());
    // -- ESD index
    df.Index->push_back(view.Entry);  // NOTE: store current track entry as ESD index
}

void Packager::StoreMC(const View::MC::Particle& view, Storage::Vector::MC_Tracks& df, EParticle pid) {
    // -- `Vector::States`
    df.X->push_back(view.X());
    df.Y->push_back(view.Y());
    df.Z->push_back(view.Z());
    df.Px->push_back(view.Px());
    df.Py->push_back(view.Py());
    df.Pz->push_back(view.Pz());
    df.Energy->push_back(view.Energy());
    // -- Mother (`Vector::MC_Id`)
    View::MC::Particle mother{view.Source, view.MotherEntry()};
    df.Mother.McEntry->push_back(mother.Entry);
    df.Mother.PdgCode->push_back(mother.PdgCode());
    // -- `Vector::MC`
    df.McEntry->push_back(view.Entry);
    df.PdgCode->push_back(view.PdgCode());
    df.ReactionID->push_back(Truth::Particle::AsTrack_ReactionID(view, mother, pid));
    df.IsTrue->push_back(static_cast<char>(Truth::Particle::AsTrack_IsTrue(view, pid)));
    df.IsSignal->push_back(static_cast<char>(Truth::Particle::AsTrack_IsSignal(view, pid)));
    df.IsSecondary->push_back(static_cast<char>(Truth::Particle::AsTrack_IsSecondary(view, pid)));
    // -- Grandmother (`Vector::MC_Id`)
    View::MC::Particle grandmother{mother.Source, mother.MotherEntry()};
    df.GrandMother.McEntry->push_back(grandmother.Entry);
    df.GrandMother.PdgCode->push_back(grandmother.PdgCode());
}

bool Packager::PassesCuts_Proton(const View::Rec::Track& track, TH1D* cut_flow_hist) const {

    cut_flow_hist->Fill(0.);
    if (std::abs(track.NSigmaProton()) > Cuts::Proton::AbsMax_NSigmaProton) return false;
    cut_flow_hist->Fill(1.);
    // if (std::abs(track.DCAxy()) < Cuts::Proton::AbsMin_DCAxy) return false;  // TEMPORARY
    // cut_flow_hist->Fill(2.);                                                 // TEMPORARY

    return true;
}

bool Packager::PassesCuts_Kaon(const View::Rec::Track& track, TH1D* cut_flow_hist) const {

    cut_flow_hist->Fill(0.);
    if (std::abs(track.NSigmaKaon()) > Cuts::Kaon::AbsMax_NSigmaKaon) return false;
    cut_flow_hist->Fill(1.);
    // if (track.TPCSignal() < Cuts::Kaon::Min_TPCSignal) return false; // PENDING
    cut_flow_hist->Fill(2.);

    return true;
}

bool Packager::PassesCuts_Pion(const View::Rec::Track& track, TH1D* cut_flow_hist) const {

    cut_flow_hist->Fill(0.);
    if (std::abs(track.NSigmaPion()) > Cuts::Pion::AbsMax_NSigmaPion) return false;
    cut_flow_hist->Fill(1.);
    // if (std::abs(track.DCAxy()) < Cuts::Pion::AbsMin_DCAxy) return false;  // TEMPORARY
    // cut_flow_hist->Fill(2.);                                               // TEMPORARY

    return true;
}

// ## V0s ZONE ## //

void Packager::FindV0s(EParticle pid) {

    // determine rules based on V0 species
    Storage::Vector::V0s* out{nullptr};
    Storage::Vector::MC_V0s* mc_out{nullptr};
    const std::vector<Fit::Track>* vec_neg{&fVec_PiMinus};
    const std::vector<Fit::Track>* vec_pos{&fVec_PiPlus};
    switch (pid) {
        case EParticle::AntiLambda:
            out = &fOutput_AntiLambdas;
            if (IsMC()) mc_out = &fOutput_MC_AntiLambdas;
            vec_neg = &fVec_AntiProtons;
            break;
        case EParticle::Lambda:
            out = &fOutput_Lambdas;
            if (IsMC()) mc_out = &fOutput_MC_Lambdas;
            vec_pos = &fVec_Protons;
            break;
        case EParticle::KaonZeroShort:
            out = &fOutput_KaonsZeroShort;
            if (IsMC()) mc_out = &fOutput_MC_KaonsZeroShort;
            break;
        default:
            Logger::Error(__FUNCTION__, "Invalid PID {} for a V0.", Const::Particle_Acronym[pid]);
            return;
    }

    // loop over all possible pairs of tracks //
    int v0_entry{0};
    for (const auto& kf_neg : *vec_neg) {
        for (const auto& kf_pos : *vec_pos) {

            // sanity check //
            int esd_neg{kf_neg.View.Entry};
            int esd_pos{kf_pos.View.Entry};
            if (esd_neg == esd_pos) continue;

            // fit v0 //
            Fit::V0 v0{v0_entry, kf_neg, kf_pos};
            v0.DoFit(fInput_Event.MagneticField);

            // apply cuts //
            if (!PassesCuts(v0, pid)) continue;

#ifdef T2S_DEBUG
            Logger::Debug(__FUNCTION__, "idx,neg,pos={},{},{}", v0.View.Entry, esd_neg, esd_pos);
            Logger::Debug(__FUNCTION__, ";x,y,z={},{},{}", v0.X(), v0.Y(), v0.Z());
            Logger::Debug(__FUNCTION__, ";x,y,z(neg)={},{},{}", v0.Neg_PCA_XYZ()[0], v0.Neg_PCA_XYZ()[1], v0.Neg_PCA_XYZ()[2]);
            Logger::Debug(__FUNCTION__, ";x,y,z(pos)={},{},{}", v0.Pos_PCA_XYZ()[0], v0.Pos_PCA_XYZ()[1], v0.Pos_PCA_XYZ()[2]);
            Logger::Debug(__FUNCTION__, ";mass={}", v0.Mass());
            Logger::Debug(__FUNCTION__, ";dca_dau={}", v0.DCA_Daughters());
            Logger::Debug(__FUNCTION__, ";radius={}", v0.Radius2D());
            Logger::Debug(__FUNCTION__, ";dca_neg={}", v0.DCA_Neg_V0());
            Logger::Debug(__FUNCTION__, ";dca_pos={}", v0.DCA_Pos_V0());
            Logger::Debug(__FUNCTION__, ";pt={}", v0.Pt());
            Logger::Debug(__FUNCTION__, ";eta={}", v0.Eta());
            Logger::Debug(__FUNCTION__, ";qt={}", v0.ArmenterosQt());
            Logger::Debug(__FUNCTION__, ";alpha={}", v0.ArmenterosAlpha());
            Logger::Debug(__FUNCTION__, ";cpa_pv={}", v0.CPA_Point(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z));
            Logger::Debug(__FUNCTION__, ";dca_pv={}", v0.DCA_Point(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z));
            Logger::Debug(__FUNCTION__, "");
#endif
            // store //
            Store(v0, *out);

            if (IsMC()) {
                View::MC::V0 mc_v0{&fInput_MC, fInput_Tracks.McEntry->at(esd_neg), fInput_Tracks.McEntry->at(esd_pos)};
                StoreMC(mc_v0, *mc_out, pid);
            }
            ++v0_entry;
        }  // end of loop over pos
    }  // end of loop over neg
}

bool Packager::PassesCuts_Lambda(const Fit::V0& v0, TH1D* cut_flow_hist) const {

    cut_flow_hist->Fill(0.);

    double mass{v0.Mass()};  // cache
    if (mass < Cuts::Lambda::Min_Mass || mass > Cuts::Lambda::Max_Mass) return false;
    cut_flow_hist->Fill(1.);

    if (v0.DCA_Daughters() > Cuts::Lambda::Max_DCAbtwDau) return false;
    cut_flow_hist->Fill(2.);

    if (v0.Radius2D() < Cuts::Lambda::Min_Radius2D) return false;
    cut_flow_hist->Fill(3.);

    if (v0.DCA_Neg_V0() > Cuts::Lambda::Max_DCAnegV0) return false;
    cut_flow_hist->Fill(4.);

    if (v0.DCA_Pos_V0() > Cuts::Lambda::Max_DCAposV0) return false;
    cut_flow_hist->Fill(5.);

    // if (v0.Pt() < Cuts::Lambda::Min_Pt) return false; // TEMP
    // cut_flow_hist->Fill(6.); // TEMP

    if (std::abs(v0.Rapidity()) > Cuts::Lambda::AbsMax_Rapidity) return false;
    cut_flow_hist->Fill(7.);

    if (v0.AbsArmQtOverAlpha() > Cuts::Lambda::AbsMax_ArmQtOverAlpha) return false;
    cut_flow_hist->Fill(8.);

    double cpa_wrt_pv{v0.CPA_Point(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z)};  // cache
    if (cpa_wrt_pv < Cuts::Lambda::Min_CPAwrtPV || cpa_wrt_pv > Cuts::Lambda::Max_CPAwrtPV) return false;
    cut_flow_hist->Fill(9.);

    if (v0.DCA_Point(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) < Cuts::Lambda::Min_DCAwrtPV) return false;
    cut_flow_hist->Fill(10.);

    return true;
}

bool Packager::PassesCuts_KaonZeroShort(const Fit::V0& v0, TH1D* cut_flow_hist) const {

    cut_flow_hist->Fill(0.);

    if (v0.DCA_Daughters() > Cuts::KaonZeroShort::Max_DCAbtwDau) return false;
    cut_flow_hist->Fill(1.);

    // if (v0.Pt() < Cuts::KaonZeroShort::Min_Pt) return false; // TEMP
    // cut_flow_hist->Fill(2.); // TEMP

    double mass{v0.Mass()};  // cache
    if (mass < Cuts::KaonZeroShort::Min_Mass || mass > Cuts::KaonZeroShort::Max_Mass) return false;
    cut_flow_hist->Fill(3.);

    if (std::abs(v0.Rapidity()) > Cuts::KaonZeroShort::AbsMax_Rapidity) return false;
    cut_flow_hist->Fill(4.);

    if (v0.Radius2D() < Cuts::KaonZeroShort::Min_Radius2D) return false;
    cut_flow_hist->Fill(5.);

    if (v0.DCA_Neg_V0() > Cuts::KaonZeroShort::Max_DCAnegV0) return false;
    cut_flow_hist->Fill(6.);

    if (v0.DCA_Pos_V0() > Cuts::KaonZeroShort::Max_DCAposV0) return false;
    cut_flow_hist->Fill(7.);

    double cpa_wrt_pv{v0.CPA_Point(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z)};  // cache
    if (cpa_wrt_pv < Cuts::KaonZeroShort::Min_CPAwrtPV || cpa_wrt_pv > Cuts::KaonZeroShort::Max_CPAwrtPV) return false;
    cut_flow_hist->Fill(8.);

    if (v0.DCA_Point(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) < Cuts::KaonZeroShort::Min_DCAwrtPV) return false;
    cut_flow_hist->Fill(9.);

    return true;
}

void Packager::Store(const Fit::V0& v0, Storage::Vector::V0s& df) {

    // V0
    // -- `Vector::States`
    df.X->push_back(static_cast<float>(v0.X()));
    df.Y->push_back(static_cast<float>(v0.Y()));
    df.Z->push_back(static_cast<float>(v0.Z()));
    df.Px->push_back(static_cast<float>(v0.Px()));
    df.Py->push_back(static_cast<float>(v0.Py()));
    df.Pz->push_back(static_cast<float>(v0.Pz()));
    df.Energy->push_back(static_cast<float>(v0.E()));
    // -- `Vector::CovMatrices`
    df.SigmaX2->push_back(static_cast<float>(v0.GetCovariance(0)));
    df.SigmaXY->push_back(static_cast<float>(v0.GetCovariance(1)));
    df.SigmaY2->push_back(static_cast<float>(v0.GetCovariance(2)));
    df.SigmaXZ->push_back(static_cast<float>(v0.GetCovariance(3)));
    df.SigmaYZ->push_back(static_cast<float>(v0.GetCovariance(4)));
    df.SigmaZ2->push_back(static_cast<float>(v0.GetCovariance(5)));
    df.SigmaXPx->push_back(static_cast<float>(v0.GetCovariance(6)));
    df.SigmaYPx->push_back(static_cast<float>(v0.GetCovariance(7)));
    df.SigmaZPx->push_back(static_cast<float>(v0.GetCovariance(8)));
    df.SigmaPx2->push_back(static_cast<float>(v0.GetCovariance(9)));
    df.SigmaXPy->push_back(static_cast<float>(v0.GetCovariance(10)));
    df.SigmaYPy->push_back(static_cast<float>(v0.GetCovariance(11)));
    df.SigmaZPy->push_back(static_cast<float>(v0.GetCovariance(12)));
    df.SigmaPxPy->push_back(static_cast<float>(v0.GetCovariance(13)));
    df.SigmaPy2->push_back(static_cast<float>(v0.GetCovariance(14)));
    df.SigmaXPz->push_back(static_cast<float>(v0.GetCovariance(15)));
    df.SigmaYPz->push_back(static_cast<float>(v0.GetCovariance(16)));
    df.SigmaZPz->push_back(static_cast<float>(v0.GetCovariance(17)));
    df.SigmaPxPz->push_back(static_cast<float>(v0.GetCovariance(18)));
    df.SigmaPyPz->push_back(static_cast<float>(v0.GetCovariance(19)));
    df.SigmaPz2->push_back(static_cast<float>(v0.GetCovariance(20)));
    df.SigmaXE->push_back(static_cast<float>(v0.GetCovariance(21)));
    df.SigmaYE->push_back(static_cast<float>(v0.GetCovariance(22)));
    df.SigmaZE->push_back(static_cast<float>(v0.GetCovariance(23)));
    df.SigmaPxE->push_back(static_cast<float>(v0.GetCovariance(24)));
    df.SigmaPyE->push_back(static_cast<float>(v0.GetCovariance(25)));
    df.SigmaPzE->push_back(static_cast<float>(v0.GetCovariance(26)));
    df.SigmaE2->push_back(static_cast<float>(v0.GetCovariance(27)));
    // -- fit info
    df.Chi2NDF->push_back(static_cast<float>(v0.Chi2NDF()));

    // Neg Daughter (`Vector::Tracks`)
    // -- `States_NoE`
    df.Neg.X->push_back(v0.Neg.View.X());
    df.Neg.Y->push_back(v0.Neg.View.Y());
    df.Neg.Z->push_back(v0.Neg.View.Z());
    df.Neg.Px->push_back(v0.Neg.View.Px());
    df.Neg.Py->push_back(v0.Neg.View.Py());
    df.Neg.Pz->push_back(v0.Neg.View.Pz());
    // -- `CovMatrices_NoE`
    df.Neg.SigmaX2->push_back(v0.Neg.View.SigmaX2());
    df.Neg.SigmaXY->push_back(v0.Neg.View.SigmaXY());
    df.Neg.SigmaY2->push_back(v0.Neg.View.SigmaY2());
    df.Neg.SigmaXZ->push_back(v0.Neg.View.SigmaXZ());
    df.Neg.SigmaYZ->push_back(v0.Neg.View.SigmaYZ());
    df.Neg.SigmaZ2->push_back(v0.Neg.View.SigmaZ2());
    df.Neg.SigmaXPx->push_back(v0.Neg.View.SigmaXPx());
    df.Neg.SigmaYPx->push_back(v0.Neg.View.SigmaYPx());
    df.Neg.SigmaZPx->push_back(v0.Neg.View.SigmaZPx());
    df.Neg.SigmaPx2->push_back(v0.Neg.View.SigmaPx2());
    df.Neg.SigmaXPy->push_back(v0.Neg.View.SigmaXPy());
    df.Neg.SigmaYPy->push_back(v0.Neg.View.SigmaYPy());
    df.Neg.SigmaZPy->push_back(v0.Neg.View.SigmaZPy());
    df.Neg.SigmaPxPy->push_back(v0.Neg.View.SigmaPxPy());
    df.Neg.SigmaPy2->push_back(v0.Neg.View.SigmaPy2());
    df.Neg.SigmaXPz->push_back(v0.Neg.View.SigmaXPz());
    df.Neg.SigmaYPz->push_back(v0.Neg.View.SigmaYPz());
    df.Neg.SigmaZPz->push_back(v0.Neg.View.SigmaZPz());
    df.Neg.SigmaPxPz->push_back(v0.Neg.View.SigmaPxPz());
    df.Neg.SigmaPyPz->push_back(v0.Neg.View.SigmaPyPz());
    df.Neg.SigmaPz2->push_back(v0.Neg.View.SigmaPz2());
    // -- the rest
    df.Neg.Charge->push_back(v0.Neg.View.Charge());
    df.Neg.DCAxy->push_back(v0.Neg.View.DCAxy());
    df.Neg.DCAz->push_back(v0.Neg.View.DCAz());
    df.Neg.TPCSignal->push_back(v0.Neg.View.TPCSignal());
    df.Neg.NSigmaPion->push_back(v0.Neg.View.NSigmaPion());
    df.Neg.NSigmaKaon->push_back(v0.Neg.View.NSigmaKaon());
    df.Neg.NSigmaProton->push_back(v0.Neg.View.NSigmaProton());
    df.Neg.Index->push_back(v0.Neg.View.Entry);  // NOTE: store current track entry as ESD index
    // -- @ PCA w.r.t. V0 (`Storage::Vector::States_NoE`)
    df.Neg_atPCA.X->push_back(static_cast<float>(v0.Neg_PCA_XYZ()[0]));
    df.Neg_atPCA.Y->push_back(static_cast<float>(v0.Neg_PCA_XYZ()[1]));
    df.Neg_atPCA.Z->push_back(static_cast<float>(v0.Neg_PCA_XYZ()[2]));
    df.Neg_atPCA.Px->push_back(static_cast<float>(v0.Neg_PCA_PxPyPz()[0]));
    df.Neg_atPCA.Py->push_back(static_cast<float>(v0.Neg_PCA_PxPyPz()[1]));
    df.Neg_atPCA.Pz->push_back(static_cast<float>(v0.Neg_PCA_PxPyPz()[2]));

    // Pos Daughter (`Vector::Tracks`)
    // -- `States_NoE`
    df.Pos.X->push_back(v0.Pos.View.X());
    df.Pos.Y->push_back(v0.Pos.View.Y());
    df.Pos.Z->push_back(v0.Pos.View.Z());
    df.Pos.Px->push_back(v0.Pos.View.Px());
    df.Pos.Py->push_back(v0.Pos.View.Py());
    df.Pos.Pz->push_back(v0.Pos.View.Pz());
    // -- `CovMatrices_NoE`
    df.Pos.SigmaX2->push_back(v0.Pos.View.SigmaX2());
    df.Pos.SigmaXY->push_back(v0.Pos.View.SigmaXY());
    df.Pos.SigmaY2->push_back(v0.Pos.View.SigmaY2());
    df.Pos.SigmaXZ->push_back(v0.Pos.View.SigmaXZ());
    df.Pos.SigmaYZ->push_back(v0.Pos.View.SigmaYZ());
    df.Pos.SigmaZ2->push_back(v0.Pos.View.SigmaZ2());
    df.Pos.SigmaXPx->push_back(v0.Pos.View.SigmaXPx());
    df.Pos.SigmaYPx->push_back(v0.Pos.View.SigmaYPx());
    df.Pos.SigmaZPx->push_back(v0.Pos.View.SigmaZPx());
    df.Pos.SigmaPx2->push_back(v0.Pos.View.SigmaPx2());
    df.Pos.SigmaXPy->push_back(v0.Pos.View.SigmaXPy());
    df.Pos.SigmaYPy->push_back(v0.Pos.View.SigmaYPy());
    df.Pos.SigmaZPy->push_back(v0.Pos.View.SigmaZPy());
    df.Pos.SigmaPxPy->push_back(v0.Pos.View.SigmaPxPy());
    df.Pos.SigmaPy2->push_back(v0.Pos.View.SigmaPy2());
    df.Pos.SigmaXPz->push_back(v0.Pos.View.SigmaXPz());
    df.Pos.SigmaYPz->push_back(v0.Pos.View.SigmaYPz());
    df.Pos.SigmaZPz->push_back(v0.Pos.View.SigmaZPz());
    df.Pos.SigmaPxPz->push_back(v0.Pos.View.SigmaPxPz());
    df.Pos.SigmaPyPz->push_back(v0.Pos.View.SigmaPyPz());
    df.Pos.SigmaPz2->push_back(v0.Pos.View.SigmaPz2());
    // -- the rest
    df.Pos.Charge->push_back(v0.Pos.View.Charge());
    df.Pos.DCAxy->push_back(v0.Pos.View.DCAxy());
    df.Pos.DCAz->push_back(v0.Pos.View.DCAz());
    df.Pos.TPCSignal->push_back(v0.Pos.View.TPCSignal());
    df.Pos.NSigmaPion->push_back(v0.Pos.View.NSigmaPion());
    df.Pos.NSigmaKaon->push_back(v0.Pos.View.NSigmaKaon());
    df.Pos.NSigmaProton->push_back(v0.Pos.View.NSigmaProton());
    df.Pos.Index->push_back(v0.Pos.View.Entry);  // NOTE: store current track entry as ESD index
    // -- @ PCA w.r.t. V0 (`Storage::Vector::States_NoE`)
    df.Pos_atPCA.X->push_back(static_cast<float>(v0.Pos_PCA_XYZ()[0]));
    df.Pos_atPCA.Y->push_back(static_cast<float>(v0.Pos_PCA_XYZ()[1]));
    df.Pos_atPCA.Z->push_back(static_cast<float>(v0.Pos_PCA_XYZ()[2]));
    df.Pos_atPCA.Px->push_back(static_cast<float>(v0.Pos_PCA_PxPyPz()[0]));
    df.Pos_atPCA.Py->push_back(static_cast<float>(v0.Pos_PCA_PxPyPz()[1]));
    df.Pos_atPCA.Pz->push_back(static_cast<float>(v0.Pos_PCA_PxPyPz()[2]));
}

void Packager::StoreMC(const View::MC::V0& v0_view, Storage::Vector::MC_V0s& df, EParticle v0_pid) {
    // V0
    // -- `Vector::States`
    df.X->push_back(v0_view.X());
    df.Y->push_back(v0_view.Y());
    df.Z->push_back(v0_view.Z());
    df.Px->push_back(v0_view.Px());
    df.Py->push_back(v0_view.Py());
    df.Pz->push_back(v0_view.Pz());
    df.Energy->push_back(v0_view.Energy());
    // -- Mother (`Vector::MC_Id`)
    View::MC::Particle mother{v0_view.Source, v0_view.MotherEntry()};  // NOTE: no hypothesis
    df.Mother.McEntry->push_back(mother.Entry);
    df.Mother.PdgCode->push_back(mother.PdgCode());
    // -- `Vector::MC`
    df.McEntry->push_back(v0_view.Entry);
    df.PdgCode->push_back(v0_view.PdgCode());
    df.ReactionID->push_back(Truth::Particle::AsV0_ReactionID(v0_view, v0_pid));
    df.IsTrue->push_back(static_cast<char>(Truth::Particle::AsV0_IsTrue(v0_view, v0_pid)));
    df.IsSignal->push_back(static_cast<char>(Truth::Particle::AsV0_IsSignal(v0_view, v0_pid)));
    df.IsSecondary->push_back(static_cast<char>(Truth::Particle::AsV0_IsSecondary(v0_view, v0_pid)));
    // -- `Vector::Coordinates`
    df.AtDecay.X->push_back(Truth::Particle::AsV0_DecayX(v0_view));
    df.AtDecay.Y->push_back(Truth::Particle::AsV0_DecayY(v0_view));
    df.AtDecay.Z->push_back(Truth::Particle::AsV0_DecayZ(v0_view));
    // -- final flag
    df.IsHybrid->push_back(static_cast<char>(Truth::Particle::AsV0_IsHybrid(v0_view, v0_pid)));

    // Neg
    // -- `Vector::PxPyPz`
    df.Neg_Momentum.Px->push_back(v0_view.Neg.Px());
    df.Neg_Momentum.Py->push_back(v0_view.Neg.Py());
    df.Neg_Momentum.Pz->push_back(v0_view.Neg.Pz());
    // -- `Vector::MC`
    df.Neg.McEntry->push_back(v0_view.Neg.Entry);
    df.Neg.PdgCode->push_back(v0_view.Neg.PdgCode());
    df.Neg.ReactionID->push_back(Truth::Particle::AsTrack_ReactionID(v0_view.Neg, v0_view, Const::V0_NegativePID[v0_pid]));
    df.Neg.IsTrue->push_back(static_cast<char>(Truth::Particle::AsTrack_IsTrue(v0_view.Neg, Const::V0_NegativePID[v0_pid])));
    df.Neg.IsSignal->push_back(static_cast<char>(Truth::Particle::AsTrack_IsSignal(v0_view.Neg, Const::V0_NegativePID[v0_pid])));
    df.Neg.IsSecondary->push_back(static_cast<char>(Truth::Particle::AsTrack_IsSecondary(v0_view.Neg, Const::V0_NegativePID[v0_pid])));

    // Pos
    // -- `Vector::PxPyPz`
    df.Pos_Momentum.Px->push_back(v0_view.Pos.Px());
    df.Pos_Momentum.Py->push_back(v0_view.Pos.Py());
    df.Pos_Momentum.Pz->push_back(v0_view.Pos.Pz());
    // -- `Vector::MC`
    df.Pos.McEntry->push_back(v0_view.Pos.Entry);
    df.Pos.PdgCode->push_back(v0_view.Pos.PdgCode());
    df.Pos.ReactionID->push_back(Truth::Particle::AsTrack_ReactionID(v0_view.Pos, v0_view, Const::V0_PositivePID[v0_pid]));
    df.Pos.IsTrue->push_back(static_cast<char>(Truth::Particle::AsTrack_IsTrue(v0_view.Pos, Const::V0_PositivePID[v0_pid])));
    df.Pos.IsSignal->push_back(static_cast<char>(Truth::Particle::AsTrack_IsSignal(v0_view.Pos, Const::V0_PositivePID[v0_pid])));
    df.Pos.IsSecondary->push_back(static_cast<char>(Truth::Particle::AsTrack_IsSecondary(v0_view.Pos, Const::V0_PositivePID[v0_pid])));
}

// ## END OF CYCLES ## //

void Packager::EndOfEvent() {
    // fill tree
    fOutputTree->Fill();
    // clear temporary containers
    fVec_SV_X.clear();
    fVec_SV_Y.clear();
    fVec_SV_Z.clear();

    fVec_AntiProtons.clear();
    fVec_Protons.clear();
    fVec_NegKaons.clear();
    fVec_PosKaons.clear();
    fVec_PiMinus.clear();
    fVec_PiPlus.clear();
    // clear output vector branches
    if (IsMC()) fOutput_Injected.Clear_VectorInjected(true);
    // based on channel
    switch (GetReactionChannel()) {
        case EReactionChannel::A:
            fOutput_AntiLambdas.Clear_VectorV0s();
            fOutput_Lambdas.Clear_VectorV0s();
            fOutput_KaonsZeroShort.Clear_VectorV0s();
            if (IsMC()) {
                fOutput_MC_AntiLambdas.Clear_VectorMC_V0s();
                fOutput_MC_Lambdas.Clear_VectorMC_V0s();
                fOutput_MC_KaonsZeroShort.Clear_VectorMC_V0s();
            }
            break;
        case EReactionChannel::D:
            fOutput_AntiLambdas.Clear_VectorV0s();
            fOutput_Lambdas.Clear_VectorV0s();
            fOutput_NegKaons.Clear_VectorTracks(false, true);
            fOutput_PosKaons.Clear_VectorTracks(false, true);
            if (IsMC()) {
                fOutput_MC_AntiLambdas.Clear_VectorMC_V0s();
                fOutput_MC_Lambdas.Clear_VectorMC_V0s();
                fOutput_MC_NegKaons.Clear_VectorMC_Tracks();
                fOutput_MC_PosKaons.Clear_VectorMC_Tracks();
            }
            break;
        case EReactionChannel::E:
            fOutput_AntiLambdas.Clear_VectorV0s();
            fOutput_Lambdas.Clear_VectorV0s();
            fOutput_NegKaons.Clear_VectorTracks(false, true);
            fOutput_PosKaons.Clear_VectorTracks(false, true);
            fOutput_PiMinus.Clear_VectorTracks(false, true);
            fOutput_PiPlus.Clear_VectorTracks(false, true);
            if (IsMC()) {
                fOutput_MC_AntiLambdas.Clear_VectorMC_V0s();
                fOutput_MC_Lambdas.Clear_VectorMC_V0s();
                fOutput_MC_NegKaons.Clear_VectorMC_Tracks();
                fOutput_MC_PosKaons.Clear_VectorMC_Tracks();
                fOutput_MC_PiMinus.Clear_VectorMC_Tracks();
                fOutput_MC_PiPlus.Clear_VectorMC_Tracks();
            }
            break;
        case EReactionChannel::H:
            fOutput_NegKaons.Clear_VectorTracks(false, true);
            fOutput_PosKaons.Clear_VectorTracks(false, true);
            if (IsMC()) {
                fOutput_MC_NegKaons.Clear_VectorMC_Tracks();
                fOutput_MC_PosKaons.Clear_VectorMC_Tracks();
            }
            break;
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
        case EReactionChannel::H:
            break;
    }

    fInputChain_Events->ResetBranchAddresses();
    fOutputTree->ResetBranchAddresses();

    Logger::Info(__FUNCTION__, "All done.");
}

}  // namespace Tree2Secondaries
