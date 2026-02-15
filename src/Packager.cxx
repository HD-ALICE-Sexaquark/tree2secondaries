#include <array>
#include <filesystem>

#include "App/Logger.hxx"
#include "KalmanFitter/KalmanFitterV0.hxx"
#include "Math/Constants.hxx"
#include "Packager/Packager.hxx"
#include "Packager/PackagerCuts.hxx"
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
    if (fOutputTree == nullptr) {
        Logger::Error(__FUNCTION__, "Couldn't create TTree \"{}\"", Const::TreeName_PackedEvents);
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
            fOutput_AntiLambdas.CreateBranches_VectorV0s(fOutputTree.get(), Const::V0_Acronym[PID_V0::AntiLambda]);
            fOutput_Lambdas.CreateBranches_VectorV0s(fOutputTree.get(), Const::V0_Acronym[PID_V0::Lambda]);
            fOutput_KaonsZeroShort.CreateBranches_VectorV0s(fOutputTree.get(), Const::V0_Acronym[PID_V0::KaonZeroShort]);
            if (IsMC()) {
                fOutput_MC_AntiLambdas.CreateBranches_VectorMC_V0s(fOutputTree.get(), Const::V0_Acronym[PID_V0::AntiLambda]);
                fOutput_MC_Lambdas.CreateBranches_VectorMC_V0s(fOutputTree.get(), Const::V0_Acronym[PID_V0::Lambda]);
                fOutput_MC_KaonsZeroShort.CreateBranches_VectorMC_V0s(fOutputTree.get(), Const::V0_Acronym[PID_V0::KaonZeroShort]);
            }
            break;
        case EReactionChannel::D:
            fOutput_AntiLambdas.CreateBranches_VectorV0s(fOutputTree.get(), Const::V0_Acronym[PID_V0::AntiLambda]);
            fOutput_Lambdas.CreateBranches_VectorV0s(fOutputTree.get(), Const::V0_Acronym[PID_V0::Lambda]);
            fOutput_NegKaons.CreateBranches_VectorTracks(fOutputTree.get(), false, true, Const::Particle_Acronym[PID_StableParticle::NegKaon]);
            fOutput_PosKaons.CreateBranches_VectorTracks(fOutputTree.get(), false, true, Const::Particle_Acronym[PID_StableParticle::PosKaon]);
            if (IsMC()) {
                fOutput_MC_AntiLambdas.CreateBranches_VectorMC_V0s(fOutputTree.get(), Const::V0_Acronym[PID_V0::AntiLambda]);
                fOutput_MC_Lambdas.CreateBranches_VectorMC_V0s(fOutputTree.get(), Const::V0_Acronym[PID_V0::Lambda]);
                fOutput_MC_NegKaons.CreateBranches_VectorMC_Tracks(fOutputTree.get(), Const::Particle_Acronym[PID_StableParticle::NegKaon]);
                fOutput_MC_PosKaons.CreateBranches_VectorMC_Tracks(fOutputTree.get(), Const::Particle_Acronym[PID_StableParticle::PosKaon]);
            }
            break;
        case EReactionChannel::E:
            fOutput_AntiLambdas.CreateBranches_VectorV0s(fOutputTree.get(), Const::V0_Acronym[PID_V0::AntiLambda]);
            fOutput_Lambdas.CreateBranches_VectorV0s(fOutputTree.get(), Const::V0_Acronym[PID_V0::Lambda]);
            fOutput_NegKaons.CreateBranches_VectorTracks(fOutputTree.get(), false, true, Const::Particle_Acronym[PID_StableParticle::NegKaon]);
            fOutput_PosKaons.CreateBranches_VectorTracks(fOutputTree.get(), false, true, Const::Particle_Acronym[PID_StableParticle::PosKaon]);
            fOutput_PiMinus.CreateBranches_VectorTracks(fOutputTree.get(), false, true, Const::Particle_Acronym[PID_StableParticle::PiMinus]);
            fOutput_PiPlus.CreateBranches_VectorTracks(fOutputTree.get(), false, true, Const::Particle_Acronym[PID_StableParticle::PiPlus]);
            if (IsMC()) {
                fOutput_MC_AntiLambdas.CreateBranches_VectorMC_V0s(fOutputTree.get(), Const::V0_Acronym[PID_V0::AntiLambda]);
                fOutput_MC_Lambdas.CreateBranches_VectorMC_V0s(fOutputTree.get(), Const::V0_Acronym[PID_V0::Lambda]);
                fOutput_MC_NegKaons.CreateBranches_VectorMC_Tracks(fOutputTree.get(), Const::Particle_Acronym[PID_StableParticle::NegKaon]);
                fOutput_MC_PosKaons.CreateBranches_VectorMC_Tracks(fOutputTree.get(), Const::Particle_Acronym[PID_StableParticle::PosKaon]);
                fOutput_MC_PiMinus.CreateBranches_VectorMC_Tracks(fOutputTree.get(), Const::Particle_Acronym[PID_StableParticle::PiMinus]);
                fOutput_MC_PiPlus.CreateBranches_VectorMC_Tracks(fOutputTree.get(), Const::Particle_Acronym[PID_StableParticle::PiPlus]);
            }
            break;
        case EReactionChannel::H:
            fOutput_NegKaons.CreateBranches_VectorTracks(fOutputTree.get(), false, true, Const::Particle_Acronym[PID_StableParticle::NegKaon]);
            fOutput_PosKaons.CreateBranches_VectorTracks(fOutputTree.get(), false, true, Const::Particle_Acronym[PID_StableParticle::PosKaon]);
            if (IsMC()) {
                fOutput_MC_NegKaons.CreateBranches_VectorMC_Tracks(fOutputTree.get(), Const::Particle_Acronym[PID_StableParticle::NegKaon]);
                fOutput_MC_PosKaons.CreateBranches_VectorMC_Tracks(fOutputTree.get(), Const::Particle_Acronym[PID_StableParticle::PosKaon]);
            }
            break;
        default:
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
        std::format("CutFlow_{}", Const::Particle_Acronym[PID_StableParticle::AntiProton]).c_str(), hist_title.c_str(), x_nbins, x_min, x_max);
    fHist_CutFlow_Proton = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", Const::Particle_Acronym[PID_StableParticle::Proton]).c_str(), hist_title.c_str(), x_nbins, x_min, x_max);
    fHist_CutFlow_NegKaon = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", Const::Particle_Acronym[PID_StableParticle::NegKaon]).c_str(), hist_title.c_str(), x_nbins, x_min, x_max);
    fHist_CutFlow_PosKaon = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", Const::Particle_Acronym[PID_StableParticle::PosKaon]).c_str(), hist_title.c_str(), x_nbins, x_min, x_max);
    fHist_CutFlow_PiMinus = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", Const::Particle_Acronym[PID_StableParticle::PiMinus]).c_str(), hist_title.c_str(), x_nbins, x_min, x_max);
    fHist_CutFlow_PiPlus = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", Const::Particle_Acronym[PID_StableParticle::PiPlus]).c_str(), hist_title.c_str(), x_nbins, x_min, x_max);

    switch (GetReactionChannel()) {
        case EReactionChannel::A:
            fHist_CutFlow_AntiLambda = std::make_unique<TH1D>(  //
                std::format("CutFlow_{}", Const::V0_Acronym[PID_V0::AntiLambda]).c_str(), hist_title.c_str(), x_nbins, x_min, x_max);
            fHist_CutFlow_Lambda = std::make_unique<TH1D>(  //
                std::format("CutFlow_{}", Const::V0_Acronym[PID_V0::Lambda]).c_str(), hist_title.c_str(), x_nbins, x_min, x_max);
            fHist_CutFlow_KaonZeroShort = std::make_unique<TH1D>(  //
                std::format("CutFlow_{}", Const::V0_Acronym[PID_V0::KaonZeroShort]).c_str(), hist_title.c_str(), x_nbins, x_min, x_max);
            break;
        case EReactionChannel::D:
        case EReactionChannel::E:
            fHist_CutFlow_AntiLambda = std::make_unique<TH1D>(  //
                std::format("CutFlow_{}", Const::V0_Acronym[PID_V0::AntiLambda]).c_str(), hist_title.c_str(), x_nbins, x_min, x_max);
            fHist_CutFlow_Lambda = std::make_unique<TH1D>(  //
                std::format("CutFlow_{}", Const::V0_Acronym[PID_V0::Lambda]).c_str(), hist_title.c_str(), x_nbins, x_min, x_max);
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

    for (size_t mc_idx = 0; mc_idx < NumberMC(); ++mc_idx) {

        if ((*fInput_MC.MotherMcEntry)[mc_idx] != Const::DummyInt) continue;
        if ((*fInput_MC.Generator)[mc_idx] != Const::SignalGeneratorID) continue;

        int status = (*fInput_MC.Status)[mc_idx];
        if (status < reaction_id_lower || status > reaction_id_upper) continue;

        auto reaction_idx = static_cast<size_t>(status - Const::ReactionID_Offset);

        if (sv_found[reaction_idx]) continue;

        fVec_SV_X[reaction_idx] = (*fInput_MC.X)[mc_idx];
        fVec_SV_Y[reaction_idx] = (*fInput_MC.Y)[mc_idx];
        fVec_SV_Z[reaction_idx] = (*fInput_MC.Z)[mc_idx];
        sv_found[reaction_idx] = true;
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

    size_t n_reserve = NumberTracks() / 3;
    fIndices_AntiProtons.reserve(n_reserve);
    fIndices_NegKaons.reserve(n_reserve);
    fIndices_PiMinus.reserve(n_reserve);
    fIndices_Protons.reserve(n_reserve);
    fIndices_PosKaons.reserve(n_reserve);
    fIndices_PiPlus.reserve(n_reserve);

    for (size_t esd_idx = 0; esd_idx < NumberTracks(); ++esd_idx) {
        // create track reference //
        View::Rec::Track track{&fInput_Tracks, esd_idx};
        // PID and pre-selection //
        if (track.Charge() < 0) {
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
        if (track.Charge() > 0) {
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

    std::vector<size_t>* vec{nullptr};
    Storage::Vector::Tracks* out{nullptr};
    Storage::Vector::MC_Tracks* mc_out{nullptr};

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

void Packager::Store(const View::Rec::Track& track, Storage::Vector::Tracks& df) {

    // `Vector::Tracks` //

    // -- `Vector::States_NoE`
    df.X->push_back(track.X());
    df.Y->push_back(track.Y());
    df.Z->push_back(track.Z());
    df.Px->push_back(track.Px());
    df.Py->push_back(track.Py());
    df.Pz->push_back(track.Pz());
    // -- `Vector::CovMatrices_NoE`
    df.SigmaX2->push_back(track.SigmaX2());
    df.SigmaXY->push_back(track.SigmaXY());
    df.SigmaY2->push_back(track.SigmaY2());
    df.SigmaXZ->push_back(track.SigmaXZ());
    df.SigmaYZ->push_back(track.SigmaYZ());
    df.SigmaZ2->push_back(track.SigmaZ2());
    df.SigmaXPx->push_back(track.SigmaXPx());
    df.SigmaYPx->push_back(track.SigmaYPx());
    df.SigmaZPx->push_back(track.SigmaZPx());
    df.SigmaPx2->push_back(track.SigmaPx2());
    df.SigmaXPy->push_back(track.SigmaXPy());
    df.SigmaYPy->push_back(track.SigmaYPy());
    df.SigmaZPy->push_back(track.SigmaZPy());
    df.SigmaPxPy->push_back(track.SigmaPxPy());
    df.SigmaPy2->push_back(track.SigmaPy2());
    df.SigmaXPz->push_back(track.SigmaXPz());
    df.SigmaYPz->push_back(track.SigmaYPz());
    df.SigmaZPz->push_back(track.SigmaZPz());
    df.SigmaPxPz->push_back(track.SigmaPxPz());
    df.SigmaPyPz->push_back(track.SigmaPyPz());
    df.SigmaPz2->push_back(track.SigmaPz2());
    // -- rest of variables
    df.Charge->push_back(track.Charge());
    df.DCAxy->push_back(track.DCAxy());
    df.DCAz->push_back(track.DCAz());
    df.TPCSignal->push_back(track.TPCSignal());
    df.NSigmaPion->push_back(track.NSigmaPion());
    df.NSigmaKaon->push_back(track.NSigmaKaon());
    df.NSigmaProton->push_back(track.NSigmaProton());
    // -- ESD index
    df.Index->push_back(track.Entry);  // NOTE: store current track entry as ESD index
}

void Packager::StoreMC(const View::MC::Particle& mc, Storage::Vector::MC_Tracks& df, PID_StableParticle pid) {
    // -- `Vector::States`
    df.X->push_back(mc.X());
    df.Y->push_back(mc.Y());
    df.Z->push_back(mc.Z());
    df.Px->push_back(mc.Px());
    df.Py->push_back(mc.Py());
    df.Pz->push_back(mc.Pz());
    df.Energy->push_back(mc.Energy());
    // -- Mother (`Vector::MC_Id`)
    View::MC::Particle mother{mc.Source, mc.MotherMcEntry()};
    bool valid_mother = View::IsValid(mother);  // used again below for grandmother
    if (valid_mother) {
        df.Mother.McEntry->push_back(mother.Entry);
        df.Mother.PdgCode->push_back(mother.PdgCode());
    } else {
        df.Mother.McEntry->push_back(Const::DummyInt);
        df.Mother.PdgCode->push_back(Const::DummyInt);
    }
    // -- `Vector::MC`
    df.McEntry->push_back(mc.Entry);
    df.PdgCode->push_back(mc.PdgCode());
    df.ReactionID->push_back(Truth::Track::ReactionID(mc, mother, pid));
    df.IsTrue->push_back(static_cast<char>(Truth::Track::IsTrue(mc, pid)));
    df.IsSignal->push_back(static_cast<char>(Truth::Track::IsSignal(mc, pid)));
    df.IsSecondary->push_back(static_cast<char>(Truth::Track::IsSecondary(mc, pid)));
    // -- Grandmother (`Vector::MC_Id`)
    if (!valid_mother) return;  // early return
    View::MC::Particle grandmother{mother.Source, mother.MotherMcEntry()};
    if (View::IsValid(grandmother)) {
        df.GrandMother.McEntry->push_back(grandmother.Entry);
        df.GrandMother.PdgCode->push_back(grandmother.PdgCode());
    } else {
        df.GrandMother.McEntry->push_back(Const::DummyInt);
        df.GrandMother.PdgCode->push_back(Const::DummyInt);
    }
}

void Packager::StoreDummyMC(Storage::Vector::MC_Tracks& df) {
    // -- `Vector::States`
    df.X->push_back(Const::DummyFloat);
    df.Y->push_back(Const::DummyFloat);
    df.Z->push_back(Const::DummyFloat);
    df.Px->push_back(Const::DummyFloat);
    df.Py->push_back(Const::DummyFloat);
    df.Pz->push_back(Const::DummyFloat);
    df.Energy->push_back(Const::DummyFloat);
    // -- Mother (`Vector::MC_Id`)
    df.Mother.McEntry->push_back(Const::DummyInt);
    df.Mother.PdgCode->push_back(Const::DummyInt);
    // -- `Vector::MC`
    df.McEntry->push_back(Const::DummyInt);
    df.PdgCode->push_back(Const::DummyInt);
    df.ReactionID->push_back(Const::DummyInt);
    df.IsTrue->push_back(0);
    df.IsSignal->push_back(0);
    df.IsSecondary->push_back(0);
    // -- Grandmother (`Vector::MC_Id`)
    df.GrandMother.McEntry->push_back(Const::DummyInt);
    df.GrandMother.PdgCode->push_back(Const::DummyInt);
}

bool Packager::Cuts_Proton(const View::Rec::Track& track, TH1D* cut_flow_hist) const {

    cut_flow_hist->Fill(0.);
    if (std::abs(track.NSigmaProton()) > Cuts::Proton::AbsMax_NSigmaProton) return false;
    cut_flow_hist->Fill(1.);
    // if (std::abs(track.DCAxy()) < Cuts::Proton::AbsMin_DCAxy) return false;  // TEMPORARY
    // cut_flow_hist->Fill(2.);                                                 // TEMPORARY

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
    // if (std::abs(track.DCAxy()) < Cuts::Pion::AbsMin_DCAxy) return false;  // TEMPORARY
    // cut_flow_hist->Fill(2.);                                               // TEMPORARY

    return true;
}

// ## V0s ZONE ## //

void Packager::FindV0s(PID_V0 pid) {

    // determine rules based on V0 species //

    Storage::Vector::V0s* out{nullptr};
    Storage::Vector::MC_V0s* mc_out{nullptr};
    const std::vector<size_t>* vec_neg{&fIndices_PiMinus};
    const std::vector<size_t>* vec_pos{&fIndices_PiPlus};
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
            auto [seed_neg, seed_pos, deriv_neg, deriv_pos] = Legacy::HelixHelix::FullPCAs(neg, pos, mass_neg, mass_pos, fInput_Event.MagneticField);

            // fit vertex //
            auto l_fit = Legacy::Fit(neg, pos, mass_neg, mass_pos, fInput_Event.MagneticField);
            auto fit = KF::Particle::FromLegacy(l_fit);
#else
            // PCAs //
            auto [seed_neg, seed_pos, pca_cache] = Seeder::HelixHelix::FastCorrectPCAs(neg, pos, fInput_Event.MagneticField);

            // apply cuts (1) //
            if (!FastCuts(seed_neg, seed_pos, pid)) continue;

            // PCAs derivatives //
            auto [deriv_neg, deriv_pos] = Seeder::HelixHelix::ComputeDerivatives(seed_neg, seed_pos, pca_cache);

            // fit vertex //
            auto fit = KF::FitVertex(neg, pos, mass_neg, mass_pos, {seed_neg, deriv_neg}, {seed_pos, deriv_pos}, fInput_Event.MagneticField);
#endif

            // create composite particle //
            KF::V0 v0(fit, seed_neg.pca, seed_pos.pca, neg, pos);

#if T2S_DEBUG
            Logger::Debug(__FUNCTION__, "{}", fit);
            Logger::Debug(__FUNCTION__, "  neg,pos={},{}", esd_neg, esd_pos);
            Logger::Debug(__FUNCTION__, "  x,y,z(neg)={},{},{}", v0.Neg_at_PCA.X(), v0.Neg_at_PCA.Y(), v0.Neg_at_PCA.Z());
            Logger::Debug(__FUNCTION__, "  x,y,z(pos)={},{},{}", v0.Pos_at_PCA.X(), v0.Pos_at_PCA.Y(), v0.Pos_at_PCA.Z());
#if T2S_LEGACY_KF == OFF
            Logger::Debug(__FUNCTION__, "  dca_dau={}", v0.DCA_Daughters());
            Logger::Debug(__FUNCTION__, "  dca_neg={}", v0.DCA_Neg_V0());
            Logger::Debug(__FUNCTION__, "  dca_pos={}", v0.DCA_Pos_V0());
            Logger::Debug(__FUNCTION__, "  qt={}", v0.ArmenterosQt());
            Logger::Debug(__FUNCTION__, "  alpha={}", v0.ArmenterosAlpha());
            Logger::Debug(__FUNCTION__, "  cpa_pv={}", v0.CPA_Vertex(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z));
            Logger::Debug(__FUNCTION__, "  dca_pv={}", v0.DCA_Vertex(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z));
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

    double mass{v0.Mass().value_or(Const::DummyDouble)};  // cached
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

    double cpa_wrt_pv{v0.CPA_Vertex(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z)};  // cached
    if (cpa_wrt_pv < Cuts::Lambda::Min_CPAwrtPV || cpa_wrt_pv > Cuts::Lambda::Max_CPAwrtPV) return false;
    cut_flow_hist->Fill(9.);

    if (v0.DCA_Vertex(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) < Cuts::Lambda::Min_DCAwrtPV) return false;
    cut_flow_hist->Fill(10.);

    return true;
}

bool Packager::SlowCuts_KaonZeroShort(const KF::V0& v0, TH1D* cut_flow_hist) const {

    // if (v0.Pt() < Cuts::KaonZeroShort::Min_Pt) return false; // PENDING
    // cut_flow_hist->Fill(2.); // PENDING

    double mass{v0.Mass().value_or(Const::DummyDouble)};  // cached
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

    double cpa_wrt_pv{v0.CPA_Vertex(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z)};  // cached
    if (cpa_wrt_pv < Cuts::KaonZeroShort::Min_CPAwrtPV || cpa_wrt_pv > Cuts::KaonZeroShort::Max_CPAwrtPV) return false;
    cut_flow_hist->Fill(8.);

    if (v0.DCA_Vertex(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) < Cuts::KaonZeroShort::Min_DCAwrtPV) return false;
    cut_flow_hist->Fill(9.);

    return true;
}

void Packager::Store(const KF::V0& v0, Storage::Vector::V0s& df) {

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
    df.SigmaX2->push_back(static_cast<float>(v0.CovX2()));
    df.SigmaXY->push_back(static_cast<float>(v0.CovXY()));
    df.SigmaY2->push_back(static_cast<float>(v0.CovY2()));
    df.SigmaXZ->push_back(static_cast<float>(v0.CovXZ()));
    df.SigmaYZ->push_back(static_cast<float>(v0.CovYZ()));
    df.SigmaZ2->push_back(static_cast<float>(v0.CovZ2()));
    df.SigmaXPx->push_back(static_cast<float>(v0.CovXPx()));
    df.SigmaYPx->push_back(static_cast<float>(v0.CovYPx()));
    df.SigmaZPx->push_back(static_cast<float>(v0.CovZPx()));
    df.SigmaPx2->push_back(static_cast<float>(v0.CovPx2()));
    df.SigmaXPy->push_back(static_cast<float>(v0.CovXPy()));
    df.SigmaYPy->push_back(static_cast<float>(v0.CovYPy()));
    df.SigmaZPy->push_back(static_cast<float>(v0.CovZPy()));
    df.SigmaPxPy->push_back(static_cast<float>(v0.CovPxPy()));
    df.SigmaPy2->push_back(static_cast<float>(v0.CovPy2()));
    df.SigmaXPz->push_back(static_cast<float>(v0.CovXPz()));
    df.SigmaYPz->push_back(static_cast<float>(v0.CovYPz()));
    df.SigmaZPz->push_back(static_cast<float>(v0.CovZPz()));
    df.SigmaPxPz->push_back(static_cast<float>(v0.CovPxPz()));
    df.SigmaPyPz->push_back(static_cast<float>(v0.CovPyPz()));
    df.SigmaPz2->push_back(static_cast<float>(v0.CovPz2()));
    df.SigmaXE->push_back(static_cast<float>(v0.CovXE()));
    df.SigmaYE->push_back(static_cast<float>(v0.CovYE()));
    df.SigmaZE->push_back(static_cast<float>(v0.CovZE()));
    df.SigmaPxE->push_back(static_cast<float>(v0.CovPxE()));
    df.SigmaPyE->push_back(static_cast<float>(v0.CovPyE()));
    df.SigmaPzE->push_back(static_cast<float>(v0.CovPzE()));
    df.SigmaE2->push_back(static_cast<float>(v0.CovE2()));
    // -- fit info
    df.Chi2NDF->push_back(static_cast<float>(v0.Chi2NDF()));

    // Neg Daughter (`Vector::Tracks`)
    // -- `States_NoE`
    df.Neg.X->push_back(v0.Neg.X());
    df.Neg.Y->push_back(v0.Neg.Y());
    df.Neg.Z->push_back(v0.Neg.Z());
    df.Neg.Px->push_back(v0.Neg.Px());
    df.Neg.Py->push_back(v0.Neg.Py());
    df.Neg.Pz->push_back(v0.Neg.Pz());
    // -- `CovMatrices_NoE`
    df.Neg.SigmaX2->push_back(v0.Neg.SigmaX2());
    df.Neg.SigmaXY->push_back(v0.Neg.SigmaXY());
    df.Neg.SigmaY2->push_back(v0.Neg.SigmaY2());
    df.Neg.SigmaXZ->push_back(v0.Neg.SigmaXZ());
    df.Neg.SigmaYZ->push_back(v0.Neg.SigmaYZ());
    df.Neg.SigmaZ2->push_back(v0.Neg.SigmaZ2());
    df.Neg.SigmaXPx->push_back(v0.Neg.SigmaXPx());
    df.Neg.SigmaYPx->push_back(v0.Neg.SigmaYPx());
    df.Neg.SigmaZPx->push_back(v0.Neg.SigmaZPx());
    df.Neg.SigmaPx2->push_back(v0.Neg.SigmaPx2());
    df.Neg.SigmaXPy->push_back(v0.Neg.SigmaXPy());
    df.Neg.SigmaYPy->push_back(v0.Neg.SigmaYPy());
    df.Neg.SigmaZPy->push_back(v0.Neg.SigmaZPy());
    df.Neg.SigmaPxPy->push_back(v0.Neg.SigmaPxPy());
    df.Neg.SigmaPy2->push_back(v0.Neg.SigmaPy2());
    df.Neg.SigmaXPz->push_back(v0.Neg.SigmaXPz());
    df.Neg.SigmaYPz->push_back(v0.Neg.SigmaYPz());
    df.Neg.SigmaZPz->push_back(v0.Neg.SigmaZPz());
    df.Neg.SigmaPxPz->push_back(v0.Neg.SigmaPxPz());
    df.Neg.SigmaPyPz->push_back(v0.Neg.SigmaPyPz());
    df.Neg.SigmaPz2->push_back(v0.Neg.SigmaPz2());
    // -- the rest
    df.Neg.Charge->push_back(v0.Neg.Charge());
    df.Neg.DCAxy->push_back(v0.Neg.DCAxy());
    df.Neg.DCAz->push_back(v0.Neg.DCAz());
    df.Neg.TPCSignal->push_back(v0.Neg.TPCSignal());
    df.Neg.NSigmaPion->push_back(v0.Neg.NSigmaPion());
    df.Neg.NSigmaKaon->push_back(v0.Neg.NSigmaKaon());
    df.Neg.NSigmaProton->push_back(v0.Neg.NSigmaProton());
    df.Neg.Index->push_back(v0.Neg.Entry);  // NOTE: store current track entry as ESD index
    // -- @ PCA w.r.t. V0 (`Storage::Vector::States_NoE`)
    df.Neg_atPCA.X->push_back(static_cast<float>(v0.Neg_at_PCA.X()));
    df.Neg_atPCA.Y->push_back(static_cast<float>(v0.Neg_at_PCA.Y()));
    df.Neg_atPCA.Z->push_back(static_cast<float>(v0.Neg_at_PCA.Z()));
    df.Neg_atPCA.Px->push_back(static_cast<float>(v0.Neg_at_PCA.Px()));
    df.Neg_atPCA.Py->push_back(static_cast<float>(v0.Neg_at_PCA.Py()));
    df.Neg_atPCA.Pz->push_back(static_cast<float>(v0.Neg_at_PCA.Pz()));

    // Pos Daughter (`Vector::Tracks`)
    // -- `States_NoE`
    df.Pos.X->push_back(v0.Pos.X());
    df.Pos.Y->push_back(v0.Pos.Y());
    df.Pos.Z->push_back(v0.Pos.Z());
    df.Pos.Px->push_back(v0.Pos.Px());
    df.Pos.Py->push_back(v0.Pos.Py());
    df.Pos.Pz->push_back(v0.Pos.Pz());
    // -- `CovMatrices_NoE`
    df.Pos.SigmaX2->push_back(v0.Pos.SigmaX2());
    df.Pos.SigmaXY->push_back(v0.Pos.SigmaXY());
    df.Pos.SigmaY2->push_back(v0.Pos.SigmaY2());
    df.Pos.SigmaXZ->push_back(v0.Pos.SigmaXZ());
    df.Pos.SigmaYZ->push_back(v0.Pos.SigmaYZ());
    df.Pos.SigmaZ2->push_back(v0.Pos.SigmaZ2());
    df.Pos.SigmaXPx->push_back(v0.Pos.SigmaXPx());
    df.Pos.SigmaYPx->push_back(v0.Pos.SigmaYPx());
    df.Pos.SigmaZPx->push_back(v0.Pos.SigmaZPx());
    df.Pos.SigmaPx2->push_back(v0.Pos.SigmaPx2());
    df.Pos.SigmaXPy->push_back(v0.Pos.SigmaXPy());
    df.Pos.SigmaYPy->push_back(v0.Pos.SigmaYPy());
    df.Pos.SigmaZPy->push_back(v0.Pos.SigmaZPy());
    df.Pos.SigmaPxPy->push_back(v0.Pos.SigmaPxPy());
    df.Pos.SigmaPy2->push_back(v0.Pos.SigmaPy2());
    df.Pos.SigmaXPz->push_back(v0.Pos.SigmaXPz());
    df.Pos.SigmaYPz->push_back(v0.Pos.SigmaYPz());
    df.Pos.SigmaZPz->push_back(v0.Pos.SigmaZPz());
    df.Pos.SigmaPxPz->push_back(v0.Pos.SigmaPxPz());
    df.Pos.SigmaPyPz->push_back(v0.Pos.SigmaPyPz());
    df.Pos.SigmaPz2->push_back(v0.Pos.SigmaPz2());
    // -- the rest
    df.Pos.Charge->push_back(v0.Pos.Charge());
    df.Pos.DCAxy->push_back(v0.Pos.DCAxy());
    df.Pos.DCAz->push_back(v0.Pos.DCAz());
    df.Pos.TPCSignal->push_back(v0.Pos.TPCSignal());
    df.Pos.NSigmaPion->push_back(v0.Pos.NSigmaPion());
    df.Pos.NSigmaKaon->push_back(v0.Pos.NSigmaKaon());
    df.Pos.NSigmaProton->push_back(v0.Pos.NSigmaProton());
    df.Pos.Index->push_back(v0.Pos.Entry);  // NOTE: store current track entry as ESD index
    // -- @ PCA w.r.t. V0 (`Storage::Vector::States_NoE`)
    df.Pos_atPCA.X->push_back(static_cast<float>(v0.Pos_at_PCA.X()));
    df.Pos_atPCA.Y->push_back(static_cast<float>(v0.Pos_at_PCA.Y()));
    df.Pos_atPCA.Z->push_back(static_cast<float>(v0.Pos_at_PCA.Z()));
    df.Pos_atPCA.Px->push_back(static_cast<float>(v0.Pos_at_PCA.Px()));
    df.Pos_atPCA.Py->push_back(static_cast<float>(v0.Pos_at_PCA.Py()));
    df.Pos_atPCA.Pz->push_back(static_cast<float>(v0.Pos_at_PCA.Pz()));
}

void Packager::StoreMC(const View::MC::V0& v0, Storage::Vector::MC_V0s& df, PID_V0 pid) {

    // V0 //

    // -- `Vector::States`
    df.X->push_back(v0.X());
    df.Y->push_back(v0.Y());
    df.Z->push_back(v0.Z());
    df.Px->push_back(v0.Px());
    df.Py->push_back(v0.Py());
    df.Pz->push_back(v0.Pz());
    df.Energy->push_back(v0.Energy());
    // -- Mother (`Vector::MC_Id`)
    View::MC::Particle mother{v0.Source, v0.MotherMcEntry()};
    if (View::IsValid(mother)) {
        df.Mother.McEntry->push_back(mother.Entry);
        df.Mother.PdgCode->push_back(mother.PdgCode());
    } else {
        df.Mother.McEntry->push_back(Const::DummyInt);
        df.Mother.PdgCode->push_back(Const::DummyInt);
    }
    // -- `Vector::MC`
    df.McEntry->push_back(v0.Entry);
    df.PdgCode->push_back(v0.PdgCode());
    df.ReactionID->push_back(Truth::V0::ReactionID(v0, pid));
    df.IsTrue->push_back(static_cast<char>(Truth::V0::IsTrue(v0, pid)));
    df.IsSignal->push_back(static_cast<char>(Truth::V0::IsSignal(v0, pid)));
    df.IsSecondary->push_back(static_cast<char>(Truth::V0::IsSecondary(v0, pid)));
    // -- `Vector::Coordinates`
    df.AtDecay.X->push_back(Truth::V0::DecayX(v0));
    df.AtDecay.Y->push_back(Truth::V0::DecayY(v0));
    df.AtDecay.Z->push_back(Truth::V0::DecayZ(v0));
    // -- final flag
    df.IsHybrid->push_back(static_cast<char>(Truth::V0::IsHybrid(v0, pid)));

    // Negative Daughter //

    // -- `Vector::PxPyPz`
    df.Neg_Momentum.Px->push_back(v0.Neg.Px());
    df.Neg_Momentum.Py->push_back(v0.Neg.Py());
    df.Neg_Momentum.Pz->push_back(v0.Neg.Pz());
    // -- `Vector::MC`
    df.Neg.McEntry->push_back(v0.Neg.Entry);
    df.Neg.PdgCode->push_back(v0.Neg.PdgCode());
    df.Neg.ReactionID->push_back(Truth::Track::ReactionID(v0.Neg, v0, Const::V0_NegativePID[pid]));
    df.Neg.IsTrue->push_back(static_cast<char>(Truth::Track::IsTrue(v0.Neg, Const::V0_NegativePID[pid])));
    df.Neg.IsSignal->push_back(static_cast<char>(Truth::Track::IsSignal(v0.Neg, Const::V0_NegativePID[pid])));
    df.Neg.IsSecondary->push_back(static_cast<char>(Truth::Track::IsSecondary(v0.Neg, Const::V0_NegativePID[pid])));

    // Positive Daughter //

    // -- `Vector::PxPyPz`
    df.Pos_Momentum.Px->push_back(v0.Pos.Px());
    df.Pos_Momentum.Py->push_back(v0.Pos.Py());
    df.Pos_Momentum.Pz->push_back(v0.Pos.Pz());
    // -- `Vector::MC`
    df.Pos.McEntry->push_back(v0.Pos.Entry);
    df.Pos.PdgCode->push_back(v0.Pos.PdgCode());
    df.Pos.ReactionID->push_back(Truth::Track::ReactionID(v0.Pos, v0, Const::V0_PositivePID[pid]));
    df.Pos.IsTrue->push_back(static_cast<char>(Truth::Track::IsTrue(v0.Pos, Const::V0_PositivePID[pid])));
    df.Pos.IsSignal->push_back(static_cast<char>(Truth::Track::IsSignal(v0.Pos, Const::V0_PositivePID[pid])));
    df.Pos.IsSecondary->push_back(static_cast<char>(Truth::Track::IsSecondary(v0.Pos, Const::V0_PositivePID[pid])));
}

void Packager::StoreDummyMC(Storage::Vector::MC_V0s& df) {

    // V0 //

    // -- `Vector::States`
    df.X->push_back(Const::DummyFloat);
    df.Y->push_back(Const::DummyFloat);
    df.Z->push_back(Const::DummyFloat);
    df.Px->push_back(Const::DummyFloat);
    df.Py->push_back(Const::DummyFloat);
    df.Pz->push_back(Const::DummyFloat);
    df.Energy->push_back(Const::DummyFloat);
    // -- Mother (`Vector::MC_Id`)
    df.Mother.McEntry->push_back(Const::DummyInt);
    df.Mother.PdgCode->push_back(Const::DummyInt);
    // -- `Vector::MC`
    df.McEntry->push_back(Const::DummyInt);
    df.PdgCode->push_back(Const::DummyInt);
    df.ReactionID->push_back(Const::DummyInt);
    df.IsTrue->push_back(0);
    df.IsSignal->push_back(0);
    df.IsSecondary->push_back(0);
    // -- `Vector::Coordinates`
    df.AtDecay.X->push_back(Const::DummyFloat);
    df.AtDecay.Y->push_back(Const::DummyFloat);
    df.AtDecay.Z->push_back(Const::DummyFloat);
    // -- final flag
    df.IsHybrid->push_back(0);

    // Negative Daughter //

    // -- `Vector::PxPyPz`
    df.Neg_Momentum.Px->push_back(Const::DummyFloat);
    df.Neg_Momentum.Py->push_back(Const::DummyFloat);
    df.Neg_Momentum.Pz->push_back(Const::DummyFloat);
    // -- `Vector::MC`
    df.Neg.McEntry->push_back(Const::DummyInt);
    df.Neg.PdgCode->push_back(Const::DummyInt);
    df.Neg.ReactionID->push_back(Const::DummyInt);
    df.Neg.IsTrue->push_back(0);
    df.Neg.IsSignal->push_back(0);
    df.Neg.IsSecondary->push_back(0);

    // Positive Daughter //

    // -- `Vector::PxPyPz`
    df.Pos_Momentum.Px->push_back(Const::DummyFloat);
    df.Pos_Momentum.Py->push_back(Const::DummyFloat);
    df.Pos_Momentum.Pz->push_back(Const::DummyFloat);
    // -- `Vector::MC`
    df.Pos.McEntry->push_back(Const::DummyInt);
    df.Pos.PdgCode->push_back(Const::DummyInt);
    df.Pos.ReactionID->push_back(Const::DummyInt);
    df.Pos.IsTrue->push_back(0);
    df.Pos.IsSignal->push_back(0);
    df.Pos.IsSecondary->push_back(0);
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

    if (IsMC()) fOutput_Injected.Clear_VectorInjected(true);

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
        default:
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
        default:
            break;
    }

    fInputChain_Events->ResetBranchAddresses();
    fOutputTree->ResetBranchAddresses();

    Logger::Info(__FUNCTION__, "All done.");
}

}  // namespace Tree2Secondaries
