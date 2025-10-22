#include <filesystem>
#include <memory>

#include <KFParticle_Math.hxx>

#include "App/Logger.hxx"
#include "App/Utilities.hxx"
#include "KF/Utilities.hxx"
#include "Math/Constants.hxx"
#include "Packager/Cuts.hxx"
#include "Packager/Packager.hxx"

namespace Tree2Secondaries {

bool Packager::Initialize() {

    fInputChain_Events = std::make_unique<TChain>("Events");
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

    CreateCutFlowHistograms();

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
    ReadBranches_Events();
    if (IsMC()) {
        ReadBranches_MC();
        ReadBranches_Injected();
    }
    ReadBranches_Tracks();
}

void Packager::ReadBranches_Events() { fInput_Event.ReadBranches_Event(fInputChain_Events.get(), IsMC()); }

void Packager::ReadBranches_Injected() {
    // `DF::SOV::Injected`
    Utils::ReadBranch(fInputChain_Events.get(), "ReactionID", &fInput_Injected.ReactionID);
    Utils::ReadBranch(fInputChain_Events.get(), "Sexaquark_Px", &fInput_Injected.Px);
    Utils::ReadBranch(fInputChain_Events.get(), "Sexaquark_Py", &fInput_Injected.Py);
    Utils::ReadBranch(fInputChain_Events.get(), "Sexaquark_Pz", &fInput_Injected.Pz);
    Utils::ReadBranch(fInputChain_Events.get(), "Nucleon_Px", &fInput_Injected.Nucleon_Px);
    Utils::ReadBranch(fInputChain_Events.get(), "Nucleon_Py", &fInput_Injected.Nucleon_Py);
    Utils::ReadBranch(fInputChain_Events.get(), "Nucleon_Pz", &fInput_Injected.Nucleon_Pz);
}

void Packager::ReadBranches_MC() { fInput_MC.ReadBranches_MCParticle(fInputChain_Events.get()); }

void Packager::ReadBranches_Tracks() { fInput_Tracks.ReadBranches_MCParticle(fInputChain_Events.get(), IsMC()); }

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

    fOutputTree = std::make_unique<TTree>("PackedEvents", "Packed Events");
    if (!fOutputTree) {
        Logger::Error(__FUNCTION__, "Couldn't create TTree \"PackedEvents\"");
        return false;
    }

    return true;
}

void Packager::CreateOutputBranches() {

    CreateOutputBranches_Events();
    if (IsMC()) CreateOutputBranches_Injected();

    switch (GetReactionChannel()) {
        // standard channels //
        case EReactionChannel::A:
            CreateOutputBranches_V0s(EParticle::AntiLambda, fOutput_AntiLambdas);
            CreateOutputBranches_V0s(EParticle::KaonZeroShort, fOutput_KaonsZeroShort);
            if (IsMC()) {
                CreateOutputBranches_LinkedV0s(EParticle::AntiLambda, fOutput_Linked_AntiLambdas);
                CreateOutputBranches_LinkedV0s(EParticle::KaonZeroShort, fOutput_Linked_KaonsZeroShort);
            }
            break;
        case EReactionChannel::D:
            CreateOutputBranches_V0s(EParticle::AntiLambda, fOutput_AntiLambdas);
            CreateOutputBranches_Tracks(EParticle::PosKaon, fOutput_PosKaons);
            if (IsMC()) {
                CreateOutputBranches_LinkedV0s(EParticle::AntiLambda, fOutput_Linked_AntiLambdas);
                CreateOutputBranches_LinkedTracks(EParticle::PosKaon, fOutput_Linked_PosKaons);
            }
            break;
        case EReactionChannel::E:
            CreateOutputBranches_V0s(EParticle::AntiLambda, fOutput_AntiLambdas);
            CreateOutputBranches_Tracks(EParticle::PosKaon, fOutput_PosKaons);
            CreateOutputBranches_Tracks(EParticle::PiMinus, fOutput_PiMinus);
            CreateOutputBranches_Tracks(EParticle::PiPlus, fOutput_PiPlus);
            if (IsMC()) {
                CreateOutputBranches_LinkedV0s(EParticle::AntiLambda, fOutput_Linked_AntiLambdas);
                CreateOutputBranches_LinkedTracks(EParticle::PosKaon, fOutput_Linked_PosKaons);
                CreateOutputBranches_LinkedTracks(EParticle::PiMinus, fOutput_Linked_PiMinus);
                CreateOutputBranches_LinkedTracks(EParticle::PiPlus, fOutput_Linked_PiPlus);
            }
            break;
        case EReactionChannel::H:
            CreateOutputBranches_Tracks(EParticle::PosKaon, fOutput_PosKaons);
            if (IsMC()) CreateOutputBranches_LinkedTracks(EParticle::PosKaon, fOutput_Linked_PosKaons);
            break;
        // anti-channels //
        case EReactionChannel::AntiA:
            CreateOutputBranches_V0s(EParticle::Lambda, fOutput_Lambdas);
            CreateOutputBranches_V0s(EParticle::KaonZeroShort, fOutput_KaonsZeroShort);
            if (IsMC()) {
                CreateOutputBranches_V0s(EParticle::Lambda, fOutput_Lambdas);
                CreateOutputBranches_V0s(EParticle::KaonZeroShort, fOutput_KaonsZeroShort);
            }
            break;
        case EReactionChannel::AntiD:
            CreateOutputBranches_V0s(EParticle::Lambda, fOutput_Lambdas);
            CreateOutputBranches_Tracks(EParticle::NegKaon, fOutput_NegKaons);
            if (IsMC()) {
                CreateOutputBranches_LinkedV0s(EParticle::Lambda, fOutput_Linked_Lambdas);
                CreateOutputBranches_LinkedTracks(EParticle::NegKaon, fOutput_Linked_NegKaons);
            }
            break;
        case EReactionChannel::AntiE:
            CreateOutputBranches_V0s(EParticle::Lambda, fOutput_Lambdas);
            CreateOutputBranches_Tracks(EParticle::NegKaon, fOutput_NegKaons);
            CreateOutputBranches_Tracks(EParticle::PiMinus, fOutput_PiMinus);
            CreateOutputBranches_Tracks(EParticle::PiPlus, fOutput_PiPlus);
            if (IsMC()) {
                CreateOutputBranches_LinkedV0s(EParticle::Lambda, fOutput_Linked_Lambdas);
                CreateOutputBranches_LinkedTracks(EParticle::NegKaon, fOutput_Linked_NegKaons);
                CreateOutputBranches_LinkedTracks(EParticle::PiMinus, fOutput_Linked_PiMinus);
                CreateOutputBranches_LinkedTracks(EParticle::PiPlus, fOutput_Linked_PiPlus);
            }
            break;
        case EReactionChannel::AntiH:
            CreateOutputBranches_Tracks(EParticle::NegKaon, fOutput_NegKaons);
            if (IsMC()) CreateOutputBranches_LinkedTracks(EParticle::NegKaon, fOutput_Linked_NegKaons);
            break;
        // for data //
        case EReactionChannel::All:
            CreateOutputBranches_V0s(EParticle::AntiLambda, fOutput_AntiLambdas);
            CreateOutputBranches_V0s(EParticle::Lambda, fOutput_Lambdas);
            CreateOutputBranches_V0s(EParticle::KaonZeroShort, fOutput_KaonsZeroShort);
            CreateOutputBranches_Tracks(EParticle::NegKaon, fOutput_NegKaons);
            CreateOutputBranches_Tracks(EParticle::PosKaon, fOutput_PosKaons);
            CreateOutputBranches_Tracks(EParticle::PiMinus, fOutput_PiMinus);
            CreateOutputBranches_Tracks(EParticle::PiPlus, fOutput_PiPlus);
            if (IsMC()) {
                CreateOutputBranches_LinkedV0s(EParticle::AntiLambda, fOutput_Linked_AntiLambdas);
                CreateOutputBranches_LinkedV0s(EParticle::Lambda, fOutput_Linked_Lambdas);
                CreateOutputBranches_LinkedV0s(EParticle::KaonZeroShort, fOutput_Linked_KaonsZeroShort);
                CreateOutputBranches_LinkedTracks(EParticle::NegKaon, fOutput_Linked_NegKaons);
                CreateOutputBranches_LinkedTracks(EParticle::PosKaon, fOutput_Linked_PosKaons);
                CreateOutputBranches_LinkedTracks(EParticle::PiMinus, fOutput_Linked_PiMinus);
                CreateOutputBranches_LinkedTracks(EParticle::PiPlus, fOutput_Linked_PiPlus);
            }
            break;
    }  // end of switch statement
}

void Packager::CreateOutputBranches_Events() { fOutput_Event.CreateBranches_Event(fOutputTree.get(), IsMC()); }

void Packager::CreateOutputBranches_Injected() {
    // `DF::SOV::Injected`
    Utils::CreateBranch(fOutputTree.get(), "ReactionID", &fOutput_Injected.ReactionID);
    Utils::CreateBranch(fOutputTree.get(), "SV_X", &fOutput_Injected.X);
    Utils::CreateBranch(fOutputTree.get(), "SV_Y", &fOutput_Injected.Y);
    Utils::CreateBranch(fOutputTree.get(), "SV_Z", &fOutput_Injected.Z);
    Utils::CreateBranch(fOutputTree.get(), "Sexaquark_Px", &fOutput_Injected.Px);
    Utils::CreateBranch(fOutputTree.get(), "Sexaquark_Py", &fOutput_Injected.Py);
    Utils::CreateBranch(fOutputTree.get(), "Sexaquark_Pz", &fOutput_Injected.Pz);
    Utils::CreateBranch(fOutputTree.get(), "Nucleon_Px", &fOutput_Injected.Nucleon_Px);
    Utils::CreateBranch(fOutputTree.get(), "Nucleon_Py", &fOutput_Injected.Nucleon_Py);
    Utils::CreateBranch(fOutputTree.get(), "Nucleon_Pz", &fOutput_Injected.Nucleon_Pz);
}

void Packager::CreateOutputBranches_V0s(EParticle pid, DF::Packed::V0s& df) {
    df.CreateBranches_PackedV0s(fOutputTree.get(), Particle::Acronym[pid]);
}

void Packager::CreateOutputBranches_Tracks(EParticle pid, DF::Packed::Tracks& df) {
    df.CreateBranches_PackedTracks(fOutputTree.get(), Particle::Acronym[pid]);
}

void Packager::CreateOutputBranches_LinkedV0s(EParticle pid, DF::Packed::LinkedV0s& df) {
    df.CreateBranches_LinkedV0s(fOutputTree.get(), Particle::Acronym[pid]);
}

void Packager::CreateOutputBranches_LinkedTracks(EParticle pid, DF::Packed::LinkedTracks& df) {
    df.CreateBranches_LinkedTracks(fOutputTree.get(), Particle::Acronym[pid]);
}

void Packager::CreateCutFlowHistograms() {

    const int x_nbins{20};
    const float x_min{0.};
    const float x_max{20.};
    std::string hist_title{";Cut N;N Passed Cut"};

    switch (GetReactionChannel()) {
        case EReactionChannel::A:
            fCutFlowHist_AntiLambdas = std::make_unique<TH1D>("CutFlow_AL", hist_title.c_str(), x_nbins, x_min, x_max);
            fCutFlowHist_KaonsZeroShort = std::make_unique<TH1D>("CutFlow_K0S", hist_title.c_str(), x_nbins, x_min, x_max);
            break;
        case EReactionChannel::D:
        case EReactionChannel::E:
            fCutFlowHist_AntiLambdas = std::make_unique<TH1D>("CutFlow_AL", hist_title.c_str(), x_nbins, x_min, x_max);
            break;
        case EReactionChannel::AntiA:
            fCutFlowHist_Lambdas = std::make_unique<TH1D>("CutFlow_L", hist_title.c_str(), x_nbins, x_min, x_max);
            fCutFlowHist_KaonsZeroShort = std::make_unique<TH1D>("CutFlow_K0S", hist_title.c_str(), x_nbins, x_min, x_max);
            break;
        case EReactionChannel::AntiD:
        case EReactionChannel::AntiE:
            fCutFlowHist_Lambdas = std::make_unique<TH1D>("CutFlow_L", hist_title.c_str(), x_nbins, x_min, x_max);
            break;
        case EReactionChannel::H:
        case EReactionChannel::AntiH:
            break;
        case EReactionChannel::All:
            fCutFlowHist_AntiLambdas = std::make_unique<TH1D>("CutFlow_AL", hist_title.c_str(), x_nbins, x_min, x_max);
            fCutFlowHist_Lambdas = std::make_unique<TH1D>("CutFlow_L", hist_title.c_str(), x_nbins, x_min, x_max);
            fCutFlowHist_KaonsZeroShort = std::make_unique<TH1D>("CutFlow_K0S", hist_title.c_str(), x_nbins, x_min, x_max);
            break;
    }
}

// ## Event ZONE ## //

void Packager::ProcessEvent() { fOutput_Event = fInput_Event; }

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

        if (fVec_SV_X[status - 600] != 0.) continue;
        fVec_SV_X[status - 600] = fInput_MC.X->at(idx_mc);
        fVec_SV_Y[status - 600] = fInput_MC.Y->at(idx_mc);
        fVec_SV_Z[status - 600] = fInput_MC.Z->at(idx_mc);
    }
}

void Packager::Injected_Store() {
    fOutput_Injected.ReactionID = fInput_Injected.ReactionID;
    fOutput_Injected.X = &fVec_SV_X;
    fOutput_Injected.Y = &fVec_SV_Y;
    fOutput_Injected.Z = &fVec_SV_Z;
    fOutput_Injected.Px = fInput_Injected.Px;
    fOutput_Injected.Py = fInput_Injected.Py;
    fOutput_Injected.Pz = fInput_Injected.Pz;
    fOutput_Injected.Nucleon_Px = fInput_Injected.Nucleon_Px;
    fOutput_Injected.Nucleon_Py = fInput_Injected.Nucleon_Py;
    fOutput_Injected.Nucleon_Pz = fInput_Injected.Nucleon_Pz;
}

// ## Tracks ZONE ## //

// Store tracks' ESD indices into vectors, according to their respective track PID and charge.
void Packager::ProcessTracks() {

    for (int iESD{0}; iESD < NumberTracks(); ++iESD) {
        // get charge //
        int charge{fInput_Tracks.Charge->at(iESD)};
        // PID //
        if (std::abs(fInput_Tracks.NSigmaProton->at(iESD)) < Cuts::Track::AbsMax_PID_NSigma) {
            if (charge < 0) fVec_AntiProtons.push_back(iESD);
            if (charge > 0) fVec_Protons.push_back(iESD);
        }
        if (std::abs(fInput_Tracks.NSigmaKaon->at(iESD)) < Cuts::Track::AbsMax_PID_NSigma) {
            if (charge < 0) fVec_NegKaons.push_back(iESD);
            if (charge > 0) fVec_PosKaons.push_back(iESD);
        }
        if (std::abs(fInput_Tracks.NSigmaPion->at(iESD)) < Cuts::Track::AbsMax_PID_NSigma) {
            if (charge < 0) fVec_PiMinus.push_back(iESD);
            if (charge > 0) fVec_PiPlus.push_back(iESD);
        }
    }

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
    // determine rules based on particle id
    const std::vector<int>* vec{nullptr};
    DF::Packed::Tracks* out{nullptr};
    DF::Packed::LinkedTracks* mc_out{nullptr};
    double mass{Particle::Mass[pid]};
    switch (pid) {
        case EParticle::NegKaon:
            vec = &fVec_NegKaons;
            out = &fOutput_NegKaons;
            if (IsMC()) mc_out = &fOutput_Linked_NegKaons;
            break;
        case EParticle::PosKaon:
            vec = &fVec_PosKaons;
            out = &fOutput_PosKaons;
            if (IsMC()) mc_out = &fOutput_Linked_PosKaons;
            break;
        case EParticle::PiMinus:
            vec = &fVec_PiMinus;
            out = &fOutput_PiMinus;
            if (IsMC()) mc_out = &fOutput_Linked_PiMinus;
            break;
        case EParticle::PiPlus:
            vec = &fVec_PiPlus;
            out = &fOutput_PiPlus;
            if (IsMC()) mc_out = &fOutput_Linked_PiPlus;
            break;
        default:
            return;
    }

    // loop over selected tracks //
    for (auto esd_idx : *vec) {

        // prepare kf objects //
        ::KF::Vector<6> neg_params = KF::IntoKF_States_NoE(fInput_Tracks, esd_idx);
        ::KF::SymMatrix<6> neg_cov = KF::IntoKF_CovMatrices_NoE(fInput_Tracks, esd_idx);
        KF::Track kf_track(neg_params, neg_cov, fInput_Tracks.Charge->at(esd_idx), mass, esd_idx);

        // store //
        Store(kf_track, *out);
        if (IsMC()) {
            MC::Track mc_track{fInput_MC, fInput_Tracks.McEntry->at(esd_idx), pid};
            StoreMC(mc_track, *mc_out);
        }
    }  // end of loop over selected tracks
}

void Packager::Store(const KF::Track& track, DF::Packed::Tracks& df) {
    // `DF::SOV::States_NoE`
    df.X->push_back(static_cast<float>(track.GetParameter(0)));
    df.Y->push_back(static_cast<float>(track.GetParameter(1)));
    df.Z->push_back(static_cast<float>(track.GetParameter(2)));
    df.Px->push_back(static_cast<float>(track.GetParameter(3)));
    df.Py->push_back(static_cast<float>(track.GetParameter(4)));
    df.Pz->push_back(static_cast<float>(track.GetParameter(5)));
    // `DF::SOV::CovMatrices_NoE`
    df.SigmaX2->push_back(static_cast<float>(track.GetCovariance(0)));
    df.SigmaXY->push_back(static_cast<float>(track.GetCovariance(1)));
    df.SigmaY2->push_back(static_cast<float>(track.GetCovariance(2)));
    df.SigmaXZ->push_back(static_cast<float>(track.GetCovariance(3)));
    df.SigmaYZ->push_back(static_cast<float>(track.GetCovariance(4)));
    df.SigmaZ2->push_back(static_cast<float>(track.GetCovariance(5)));
    df.SigmaXPx->push_back(static_cast<float>(track.GetCovariance(6)));
    df.SigmaYPx->push_back(static_cast<float>(track.GetCovariance(7)));
    df.SigmaZPx->push_back(static_cast<float>(track.GetCovariance(8)));
    df.SigmaPx2->push_back(static_cast<float>(track.GetCovariance(9)));
    df.SigmaXPy->push_back(static_cast<float>(track.GetCovariance(10)));
    df.SigmaYPy->push_back(static_cast<float>(track.GetCovariance(11)));
    df.SigmaZPy->push_back(static_cast<float>(track.GetCovariance(12)));
    df.SigmaPxPy->push_back(static_cast<float>(track.GetCovariance(13)));
    df.SigmaPy2->push_back(static_cast<float>(track.GetCovariance(14)));
    df.SigmaXPz->push_back(static_cast<float>(track.GetCovariance(15)));
    df.SigmaYPz->push_back(static_cast<float>(track.GetCovariance(16)));
    df.SigmaZPz->push_back(static_cast<float>(track.GetCovariance(17)));
    df.SigmaPxPz->push_back(static_cast<float>(track.GetCovariance(18)));
    df.SigmaPyPz->push_back(static_cast<float>(track.GetCovariance(19)));
    df.SigmaPz2->push_back(static_cast<float>(track.GetCovariance(20)));
    // `DF::Packed::Tracks`
    df.Entry->push_back(track.idx);
}

void Packager::StoreMC(const MC::Track& mc_track, DF::Packed::LinkedTracks& df) {
    // `DF::SOV::States`
    df.X->push_back(mc_track.X);
    df.Y->push_back(mc_track.Y);
    df.Z->push_back(mc_track.Z);
    df.Px->push_back(mc_track.Px);
    df.Py->push_back(mc_track.Py);
    df.Pz->push_back(mc_track.Pz);
    df.Energy->push_back(mc_track.Energy);
    // `DF::Packed::LinkedTracks`
    df.Entry->push_back(mc_track.Entry);
    df.PdgCode->push_back(mc_track.PdgCode);
    df.Mother_Entry->push_back(mc_track.Mother_Entry);
    df.Mother_PdgCode->push_back(mc_track.Mother_PdgCode);
    df.GrandMother_Entry->push_back(mc_track.GrandMother_Entry);
    df.GrandMother_PdgCode->push_back(mc_track.GrandMother_PdgCode);
    df.ReactionID->push_back(mc_track.ReactionID);
    df.IsTrue->push_back(static_cast<char>(mc_track.IsTrue));
    df.IsSignal->push_back(static_cast<char>(mc_track.IsSignal));
    df.IsSecondary->push_back(static_cast<char>(mc_track.IsSecondary));
}

// ## V0s ZONE ## //

void Packager::FindV0s(EParticle pid) {
    // Determine rules based on V0 pdg code
    const std::vector<int>* vec_neg{nullptr};
    const std::vector<int>* vec_pos{nullptr};
    DF::Packed::V0s* out{nullptr};
    DF::Packed::LinkedV0s* mc_out{nullptr};
    EParticle pid_neg;
    EParticle pid_pos;
    switch (pid) {
        case EParticle::AntiLambda:
            vec_neg = &fVec_AntiProtons;
            vec_pos = &fVec_PiPlus;
            out = &fOutput_AntiLambdas;
            if (IsMC()) mc_out = &fOutput_Linked_AntiLambdas;
            pid_neg = EParticle::AntiProton;
            pid_pos = EParticle::PiPlus;
            break;
        case EParticle::Lambda:
            vec_neg = &fVec_PiMinus;
            vec_pos = &fVec_Protons;
            out = &fOutput_Lambdas;
            if (IsMC()) mc_out = &fOutput_Linked_Lambdas;
            pid_neg = EParticle::PiMinus;
            pid_pos = EParticle::Proton;
            break;
        case EParticle::KaonZeroShort:
            vec_neg = &fVec_PiMinus;
            vec_pos = &fVec_PiPlus;
            out = &fOutput_KaonsZeroShort;
            if (IsMC()) mc_out = &fOutput_Linked_KaonsZeroShort;
            pid_neg = EParticle::PiMinus;
            pid_pos = EParticle::PiPlus;
            break;
        default:
            Logger::Error(__FUNCTION__, "Invalid PID {} for a V0.", Particle::Acronym[pid]);
            return;
    }
    double mass_neg{Particle::Mass[pid_neg]};
    double mass_pos{Particle::Mass[pid_pos]};

    // loop over all possible pairs of tracks //
    int v0_entry{0};
    for (auto esd_neg : *vec_neg) {
        for (auto esd_pos : *vec_pos) {

            // sanity check //
            if (esd_neg == esd_pos) continue;

            // prepare neg //
            ::KF::Vector<6> neg_params = KF::IntoKF_States_NoE(fInput_Tracks, esd_neg);
            ::KF::SymMatrix<6> neg_cov = KF::IntoKF_CovMatrices_NoE(fInput_Tracks, esd_neg);
            KF::Track neg(neg_params, neg_cov, fInput_Tracks.Charge->at(esd_neg), mass_neg, esd_neg);

            // prepare pos //
            ::KF::Vector<6> pos_params = KF::IntoKF_States_NoE(fInput_Tracks, esd_pos);
            ::KF::SymMatrix<6> pos_cov = KF::IntoKF_CovMatrices_NoE(fInput_Tracks, esd_pos);
            KF::Track pos(pos_params, pos_cov, fInput_Tracks.Charge->at(esd_pos), mass_pos, esd_pos);

            // fit v0 //
            KF::V0 v0{v0_entry, pid, neg, pos};
            v0.DoFit(fInput_Event.MagneticField);
            // apply cuts //
            if (!PassesCuts(v0, pid)) continue;
#ifdef T2S_DEBUG
            Logger::Debug(__FUNCTION__, "idx,neg,pos={},{},{}", v0.idx, neg.idx, pos.idx);
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
                MC::V0 mc_v0{fInput_MC, fInput_Tracks.McEntry->at(esd_neg), fInput_Tracks.McEntry->at(esd_pos), pid, pid_neg, pid_pos};
                StoreMC(mc_v0, *mc_out);
            }
            ++v0_entry;
        }  // end of loop over pos
    }  // end of loop over neg
}

bool Packager::PassesCuts_Lambda(const KF::V0& v0) const {

    const auto& CutFlowHist{v0.hypothesis == EParticle::AntiLambda ? fCutFlowHist_AntiLambdas : fCutFlowHist_Lambdas};

    CutFlowHist->Fill(0.);
    if (v0.Mass() < Cuts::Lambda::Min_Mass || v0.Mass() > Cuts::Lambda::Max_Mass) return false;
    CutFlowHist->Fill(1.);
    if (v0.DCA_Daughters() > Cuts::Lambda::Max_DCAbtwDau) return false;
    CutFlowHist->Fill(2.);
    if (v0.AbsZ() > Cuts::Lambda::AbsMax_Zv) return false;
    CutFlowHist->Fill(3.);
    if (v0.Radius2D() < Cuts::Lambda::Min_Radius2D || v0.Radius2D() > Cuts::Lambda::Max_Radius2D) return false;
    CutFlowHist->Fill(4.);
    if (v0.DCA_Neg_V0() > Cuts::Lambda::Max_DCAnegV0) return false;
    CutFlowHist->Fill(5.);
    if (v0.DCA_Pos_V0() > Cuts::Lambda::Max_DCAposV0) return false;
    CutFlowHist->Fill(6.);
    if (v0.Pt() < Cuts::Lambda::Min_Pt) return false;
    CutFlowHist->Fill(7.);
    if (v0.AbsEta() > Cuts::Lambda::AbsMax_Eta) return false;
    CutFlowHist->Fill(8.);
    if (v0.AbsArmQtOverAlpha() > Cuts::Lambda::AbsMax_ArmQtOverAlpha) return false;
    CutFlowHist->Fill(9.);
    if (v0.CPA_Point(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) < Cuts::Lambda::Min_CPAwrtPV ||
        v0.CPA_Point(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) > Cuts::Lambda::Max_CPAwrtPV) {
        return false;
    }
    CutFlowHist->Fill(10.);
    if (v0.DCA_Point(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) < Cuts::Lambda::Min_DCAwrtPV) return false;
    CutFlowHist->Fill(11.);

    return true;
}

bool Packager::PassesCuts_KaonZeroShort(const KF::V0& v0) const {

    fCutFlowHist_KaonsZeroShort->Fill(0.);
    if (v0.DCA_Daughters() > Cuts::KaonZeroShort::Max_DCAbtwDau) return false;
    fCutFlowHist_KaonsZeroShort->Fill(1.);
    if (v0.Pt() < Cuts::KaonZeroShort::Min_Pt) return false;
    fCutFlowHist_KaonsZeroShort->Fill(2.);
    if (v0.Mass() < Cuts::KaonZeroShort::Min_Mass || v0.Mass() > Cuts::KaonZeroShort::Max_Mass) return false;
    fCutFlowHist_KaonsZeroShort->Fill(3.);
    if (v0.AbsEta() > Cuts::KaonZeroShort::AbsMax_Eta) return false;
    fCutFlowHist_KaonsZeroShort->Fill(4.);
    if (v0.AbsZ() > Cuts::KaonZeroShort::AbsMax_Zv) return false;
    fCutFlowHist_KaonsZeroShort->Fill(5.);
    if (v0.Radius2D() < Cuts::KaonZeroShort::Min_Radius2D || v0.Radius2D() > Cuts::KaonZeroShort::Max_Radius2D) return false;
    fCutFlowHist_KaonsZeroShort->Fill(6.);
    if (v0.DCA_Neg_V0() > Cuts::KaonZeroShort::Max_DCAnegV0) return false;
    fCutFlowHist_KaonsZeroShort->Fill(7.);
    if (v0.DCA_Pos_V0() > Cuts::KaonZeroShort::Max_DCAposV0) return false;
    fCutFlowHist_KaonsZeroShort->Fill(8.);
    if (v0.CPA_Point(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) < Cuts::KaonZeroShort::Min_CPAwrtPV ||
        v0.CPA_Point(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) > Cuts::KaonZeroShort::Max_CPAwrtPV) {
        return false;
    }
    fCutFlowHist_KaonsZeroShort->Fill(9.);
    if (v0.DCA_Point(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) < Cuts::KaonZeroShort::Min_DCAwrtPV) return false;
    fCutFlowHist_KaonsZeroShort->Fill(10.);

    return true;
}

void Packager::Store(const KF::V0& v0, DF::Packed::V0s& df) {
    // `DF::SOV::States`
    df.X->push_back(static_cast<float>(v0.X()));
    df.Y->push_back(static_cast<float>(v0.Y()));
    df.Z->push_back(static_cast<float>(v0.Z()));
    df.Px->push_back(static_cast<float>(v0.Px()));
    df.Py->push_back(static_cast<float>(v0.Py()));
    df.Pz->push_back(static_cast<float>(v0.Pz()));
    df.Energy->push_back(static_cast<float>(v0.E()));
    // `DF::SOV::CovMatrices`
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

    // Neg Daughter
    // -- `DF::SOV::States_NoE`
    df.Neg.X->push_back(static_cast<float>(v0.Neg.X()));
    df.Neg.Y->push_back(static_cast<float>(v0.Neg.Y()));
    df.Neg.Z->push_back(static_cast<float>(v0.Neg.Z()));
    df.Neg.Px->push_back(static_cast<float>(v0.Neg.Px()));
    df.Neg.Py->push_back(static_cast<float>(v0.Neg.Py()));
    df.Neg.Pz->push_back(static_cast<float>(v0.Neg.Pz()));
    // -- `DF::SOV::CovMatrices_NoE`
    df.Neg.SigmaX2->push_back(static_cast<float>(v0.Neg.GetCovariance(0)));
    df.Neg.SigmaXY->push_back(static_cast<float>(v0.Neg.GetCovariance(1)));
    df.Neg.SigmaY2->push_back(static_cast<float>(v0.Neg.GetCovariance(2)));
    df.Neg.SigmaXZ->push_back(static_cast<float>(v0.Neg.GetCovariance(3)));
    df.Neg.SigmaYZ->push_back(static_cast<float>(v0.Neg.GetCovariance(4)));
    df.Neg.SigmaZ2->push_back(static_cast<float>(v0.Neg.GetCovariance(5)));
    df.Neg.SigmaXPx->push_back(static_cast<float>(v0.Neg.GetCovariance(6)));
    df.Neg.SigmaYPx->push_back(static_cast<float>(v0.Neg.GetCovariance(7)));
    df.Neg.SigmaZPx->push_back(static_cast<float>(v0.Neg.GetCovariance(8)));
    df.Neg.SigmaPx2->push_back(static_cast<float>(v0.Neg.GetCovariance(9)));
    df.Neg.SigmaXPy->push_back(static_cast<float>(v0.Neg.GetCovariance(10)));
    df.Neg.SigmaYPy->push_back(static_cast<float>(v0.Neg.GetCovariance(11)));
    df.Neg.SigmaZPy->push_back(static_cast<float>(v0.Neg.GetCovariance(12)));
    df.Neg.SigmaPxPy->push_back(static_cast<float>(v0.Neg.GetCovariance(13)));
    df.Neg.SigmaPy2->push_back(static_cast<float>(v0.Neg.GetCovariance(14)));
    df.Neg.SigmaXPz->push_back(static_cast<float>(v0.Neg.GetCovariance(15)));
    df.Neg.SigmaYPz->push_back(static_cast<float>(v0.Neg.GetCovariance(16)));
    df.Neg.SigmaZPz->push_back(static_cast<float>(v0.Neg.GetCovariance(17)));
    df.Neg.SigmaPxPz->push_back(static_cast<float>(v0.Neg.GetCovariance(18)));
    df.Neg.SigmaPyPz->push_back(static_cast<float>(v0.Neg.GetCovariance(19)));
    df.Neg.SigmaPz2->push_back(static_cast<float>(v0.Neg.GetCovariance(20)));
    // -- `DF::Packed::Tracks`
    df.Neg.Entry->push_back(v0.Neg.idx);

    // Neg Daughter @ PCA w.r.t. V0
    // -- `DF::SOV::States_NoE`
    df.Neg_atPCA.X->push_back(Const::DummyFloat);   // PENDING!
    df.Neg_atPCA.Y->push_back(Const::DummyFloat);   // PENDING!
    df.Neg_atPCA.Z->push_back(Const::DummyFloat);   // PENDING!
    df.Neg_atPCA.Px->push_back(Const::DummyFloat);  // PENDING!
    df.Neg_atPCA.Py->push_back(Const::DummyFloat);  // PENDING!
    df.Neg_atPCA.Pz->push_back(Const::DummyFloat);  // PENDING!

    // Pos Daughter
    // -- `DF::SOV::States_NoE`
    df.Pos.X->push_back(static_cast<float>(v0.Pos.X()));
    df.Pos.Y->push_back(static_cast<float>(v0.Pos.Y()));
    df.Pos.Z->push_back(static_cast<float>(v0.Pos.Z()));
    df.Pos.Px->push_back(static_cast<float>(v0.Pos.Px()));
    df.Pos.Py->push_back(static_cast<float>(v0.Pos.Py()));
    df.Pos.Pz->push_back(static_cast<float>(v0.Pos.Pz()));
    // -- `DF::SOV::CovMatrices_NoE`
    df.Pos.SigmaX2->push_back(static_cast<float>(v0.Pos.GetCovariance(0)));
    df.Pos.SigmaXY->push_back(static_cast<float>(v0.Pos.GetCovariance(1)));
    df.Pos.SigmaY2->push_back(static_cast<float>(v0.Pos.GetCovariance(2)));
    df.Pos.SigmaXZ->push_back(static_cast<float>(v0.Pos.GetCovariance(3)));
    df.Pos.SigmaYZ->push_back(static_cast<float>(v0.Pos.GetCovariance(4)));
    df.Pos.SigmaZ2->push_back(static_cast<float>(v0.Pos.GetCovariance(5)));
    df.Pos.SigmaXPx->push_back(static_cast<float>(v0.Pos.GetCovariance(6)));
    df.Pos.SigmaYPx->push_back(static_cast<float>(v0.Pos.GetCovariance(7)));
    df.Pos.SigmaZPx->push_back(static_cast<float>(v0.Pos.GetCovariance(8)));
    df.Pos.SigmaPx2->push_back(static_cast<float>(v0.Pos.GetCovariance(9)));
    df.Pos.SigmaXPy->push_back(static_cast<float>(v0.Pos.GetCovariance(10)));
    df.Pos.SigmaYPy->push_back(static_cast<float>(v0.Pos.GetCovariance(11)));
    df.Pos.SigmaZPy->push_back(static_cast<float>(v0.Pos.GetCovariance(12)));
    df.Pos.SigmaPxPy->push_back(static_cast<float>(v0.Pos.GetCovariance(13)));
    df.Pos.SigmaPy2->push_back(static_cast<float>(v0.Pos.GetCovariance(14)));
    df.Pos.SigmaXPz->push_back(static_cast<float>(v0.Pos.GetCovariance(15)));
    df.Pos.SigmaYPz->push_back(static_cast<float>(v0.Pos.GetCovariance(16)));
    df.Pos.SigmaZPz->push_back(static_cast<float>(v0.Pos.GetCovariance(17)));
    df.Pos.SigmaPxPz->push_back(static_cast<float>(v0.Pos.GetCovariance(18)));
    df.Pos.SigmaPyPz->push_back(static_cast<float>(v0.Pos.GetCovariance(19)));
    df.Pos.SigmaPz2->push_back(static_cast<float>(v0.Pos.GetCovariance(20)));
    // -- `DF::Packed::Tracks`
    df.Pos.Entry->push_back(v0.Pos.idx);

    // Pos Daughter @ PCA w.r.t. V0
    // -- `DF::SOV::States_NoE`
    df.Pos_atPCA.X->push_back(Const::DummyFloat);   // PENDING!
    df.Pos_atPCA.Y->push_back(Const::DummyFloat);   // PENDING!
    df.Pos_atPCA.Z->push_back(Const::DummyFloat);   // PENDING!
    df.Pos_atPCA.Px->push_back(Const::DummyFloat);  // PENDING!
    df.Pos_atPCA.Py->push_back(Const::DummyFloat);  // PENDING!
    df.Pos_atPCA.Pz->push_back(Const::DummyFloat);  // PENDING!

    // `DF::Packed::V0s`
    df.Entry->push_back(v0.idx);
    df.Chi2NDF->push_back(Const::DummyFloat);  // PENDING!
}

void Packager::StoreMC(const MC::V0& mc_v0, DF::Packed::LinkedV0s& df) {
    // `DF::SOV::States`
    df.X->push_back(mc_v0.X);
    df.Y->push_back(mc_v0.Y);
    df.Z->push_back(mc_v0.Z);
    df.Px->push_back(mc_v0.Px);
    df.Py->push_back(mc_v0.Py);
    df.Pz->push_back(mc_v0.Pz);
    df.Energy->push_back(mc_v0.Energy);
    // `DF::SOV::MCInfo`
    df.Entry->push_back(mc_v0.Entry);
    df.PdgCode->push_back(mc_v0.PdgCode);
    df.Mother_Entry->push_back(mc_v0.Mother_Entry);
    df.Mother_PdgCode->push_back(mc_v0.Mother_PdgCode);
    df.ReactionID->push_back(mc_v0.ReactionID);
    df.IsTrue->push_back(static_cast<char>(mc_v0.IsTrue));
    df.IsSignal->push_back(static_cast<char>(mc_v0.IsSignal));
    df.IsSecondary->push_back(static_cast<char>(mc_v0.IsSecondary));
    // `DF::Packed::LinkedV0s`
    df.DecayX->push_back(mc_v0.DecayX());
    df.DecayY->push_back(mc_v0.DecayY());
    df.DecayZ->push_back(mc_v0.DecayZ());
    df.IsHybrid->push_back(static_cast<char>(mc_v0.IsHybrid));
    // -- neg. daughter (`DF::SOV::MCInfo_Reduced`)
    df.Neg.Px->push_back(mc_v0.neg.Px);
    df.Neg.Py->push_back(mc_v0.neg.Py);
    df.Neg.Pz->push_back(mc_v0.neg.Pz);
    df.Neg.Entry->push_back(mc_v0.neg.Entry);
    df.Neg.PdgCode->push_back(mc_v0.neg.PdgCode);
    df.Neg.ReactionID->push_back(mc_v0.neg.ReactionID);
    df.Neg.IsTrue->push_back(static_cast<char>(mc_v0.neg.IsTrue));
    df.Neg.IsSignal->push_back(static_cast<char>(mc_v0.neg.IsSignal));
    df.Neg.IsSecondary->push_back(static_cast<char>(mc_v0.neg.IsSecondary));
    // -- pos. daughter (`DF::SOV::MCInfo_Reduced`)
    df.Pos.Px->push_back(mc_v0.pos.Px);
    df.Pos.Py->push_back(mc_v0.pos.Py);
    df.Pos.Pz->push_back(mc_v0.pos.Pz);
    df.Pos.Entry->push_back(mc_v0.pos.Entry);
    df.Pos.PdgCode->push_back(mc_v0.pos.PdgCode);
    df.Pos.ReactionID->push_back(mc_v0.pos.ReactionID);
    df.Pos.IsTrue->push_back(static_cast<char>(mc_v0.pos.IsTrue));
    df.Pos.IsSignal->push_back(static_cast<char>(mc_v0.pos.IsSignal));
    df.Pos.IsSecondary->push_back(static_cast<char>(mc_v0.pos.IsSecondary));
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
    // clear output branches
    if (IsMC()) fOutput_Injected.Clear_Injected();
    switch (GetReactionChannel()) {
        // standard channels
        case EReactionChannel::A:
            fOutput_AntiLambdas.Clear_PackedV0s();
            fOutput_KaonsZeroShort.Clear_PackedV0s();
            if (IsMC()) {
                fOutput_Linked_AntiLambdas.Clear_LinkedV0s();
                fOutput_Linked_KaonsZeroShort.Clear_LinkedV0s();
            }
            break;
        case EReactionChannel::D:
            fOutput_AntiLambdas.Clear_PackedV0s();
            fOutput_PosKaons.Clear_PackedTracks();
            if (IsMC()) {
                fOutput_Linked_AntiLambdas.Clear_LinkedV0s();
                fOutput_Linked_PosKaons.Clear_LinkedTracks();
            }
            break;
        case EReactionChannel::E:
            fOutput_AntiLambdas.Clear_PackedV0s();
            fOutput_PosKaons.Clear_PackedTracks();
            fOutput_PiMinus.Clear_PackedTracks();
            fOutput_PiPlus.Clear_PackedTracks();
            if (IsMC()) {
                fOutput_Linked_AntiLambdas.Clear_LinkedV0s();
                fOutput_Linked_PosKaons.Clear_LinkedTracks();
                fOutput_Linked_PiMinus.Clear_LinkedTracks();
                fOutput_Linked_PiPlus.Clear_LinkedTracks();
            }
            break;
        case EReactionChannel::H:
            fOutput_PosKaons.Clear_PackedTracks();
            if (IsMC()) fOutput_Linked_PosKaons.Clear_LinkedTracks();
            break;
        // anti-channels
        case EReactionChannel::AntiA:
            fOutput_Lambdas.Clear_PackedV0s();
            fOutput_KaonsZeroShort.Clear_PackedV0s();
            if (IsMC()) {
                fOutput_Linked_Lambdas.Clear_LinkedV0s();
                fOutput_Linked_KaonsZeroShort.Clear_LinkedV0s();
            }
            break;
        case EReactionChannel::AntiD:
            fOutput_Lambdas.Clear_PackedV0s();
            fOutput_NegKaons.Clear_PackedTracks();
            if (IsMC()) {
                fOutput_Linked_Lambdas.Clear_LinkedV0s();
                fOutput_Linked_NegKaons.Clear_LinkedTracks();
            }
            break;
        case EReactionChannel::AntiE:
            fOutput_Lambdas.Clear_PackedV0s();
            fOutput_NegKaons.Clear_PackedTracks();
            fOutput_PiMinus.Clear_PackedTracks();
            fOutput_PiPlus.Clear_PackedTracks();
            if (IsMC()) {
                fOutput_Linked_Lambdas.Clear_LinkedV0s();
                fOutput_Linked_NegKaons.Clear_LinkedTracks();
                fOutput_Linked_PiMinus.Clear_LinkedTracks();
                fOutput_Linked_PiPlus.Clear_LinkedTracks();
            }
            break;
        case EReactionChannel::AntiH:
            fOutput_NegKaons.Clear_PackedTracks();
            if (IsMC()) fOutput_Linked_NegKaons.Clear_LinkedTracks();
            break;
        // for data
        case EReactionChannel::All:
            fOutput_AntiLambdas.Clear_PackedV0s();
            fOutput_Lambdas.Clear_PackedV0s();
            fOutput_KaonsZeroShort.Clear_PackedV0s();
            fOutput_NegKaons.Clear_PackedTracks();
            fOutput_PosKaons.Clear_PackedTracks();
            fOutput_PiMinus.Clear_PackedTracks();
            fOutput_PiPlus.Clear_PackedTracks();
            if (IsMC()) {
                fOutput_Linked_AntiLambdas.Clear_LinkedV0s();
                fOutput_Linked_Lambdas.Clear_LinkedV0s();
                fOutput_Linked_KaonsZeroShort.Clear_LinkedV0s();
                fOutput_Linked_NegKaons.Clear_LinkedTracks();
                fOutput_Linked_PosKaons.Clear_LinkedTracks();
                fOutput_Linked_PiMinus.Clear_LinkedTracks();
                fOutput_Linked_PiPlus.Clear_LinkedTracks();
            }
            break;
    }
}

void Packager::EndOfAnalysis() {

    fOutputTree->Write();
    Logger::Info(__FUNCTION__, "TTree \"{}\" has been written into TFile {}", fOutputTree->GetName(), fSettings.PathOutputFile);

    switch (GetReactionChannel()) {
        case EReactionChannel::A:
            fCutFlowHist_AntiLambdas->Write();
            fCutFlowHist_KaonsZeroShort->Write();
            break;
        case EReactionChannel::D:
        case EReactionChannel::E:
            fCutFlowHist_AntiLambdas->Write();
            break;
        case EReactionChannel::AntiA:
            fCutFlowHist_Lambdas->Write();
            fCutFlowHist_KaonsZeroShort->Write();
            break;
        case EReactionChannel::AntiD:
        case EReactionChannel::AntiE:
            fCutFlowHist_Lambdas->Write();
            break;
        case EReactionChannel::H:
        case EReactionChannel::AntiH:
            break;
        case EReactionChannel::All:
            fCutFlowHist_AntiLambdas->Write();
            fCutFlowHist_Lambdas->Write();
            fCutFlowHist_KaonsZeroShort->Write();
            break;
    }

    fInputChain_Events->ResetBranchAddresses();
    fOutputTree->ResetBranchAddresses();

    Logger::Info(__FUNCTION__, "All done.");
}

}  // namespace Tree2Secondaries
