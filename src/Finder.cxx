#include <filesystem>
#include <memory>
#include <set>

#include "App/Logger.hxx"
#include "App/Utilities.hxx"
#include "DataFormats/PackedEvents.hxx"
#include "Finder/Cuts.hxx"
#include "Finder/Finder.hxx"
#include "KF/Utilities.hxx"
#include "MC/Particles.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries {

bool Finder::Initialize() {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    fTree_PackedEvents = std::make_unique<TChain>("PackedEvents");
    for (const auto& path : fSettings.PathInputFiles) {
        if (fTree_PackedEvents->Add(path.c_str()) == 0) {
            Logger::Error(__FUNCTION__, "Couldn't add TFile {}");
        }
    }
    if (!fTree_PackedEvents->GetEntries()) {
        Logger::Error(__FUNCTION__, "Couldn't manage to read any entry.");
        return false;
    }
    Logger::Info(__FUNCTION__, "TChain \"{}\" loaded successfully with {} trees and {} total entries.", fTree_PackedEvents->GetName(),
                 fTree_PackedEvents->GetNtrees(), fTree_PackedEvents->GetEntries());

    ConnectInputBranches();

    if (!PrepareOutputFile()) return false;

    CreateCutFlowHistogram();

    if (!PrepareOutputTree()) return false;
    CreateOutputBranches();

    if (IsMC()) {
        if (!Injected_PrepareOutputTree()) return false;
        Injected_CreateOutputBranches();
    }

    Logger::Info(__FUNCTION__, "Finder initialized successfully.");

    return true;
}

// ## INPUT ZONE ## //

void Finder::ConnectInputBranches() {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    fTree_PackedEvents->SetBranchStatus("*", false);

    ConnectBranches_Events();
    if (IsMC()) ConnectBranches_Injected();

    switch (GetReactionChannel()) {
        // standard channels //
        case EReactionChannel::A:
            ConnectBranches_V0s(EParticle::AntiLambda, fPacked_AntiLambdas);
            ConnectBranches_V0s(EParticle::KaonZeroShort, fPacked_KaonsZeroShort);
            if (IsMC()) {
                ConnectBranches_MC_V0s(EParticle::AntiLambda, fPacked_MC_AntiLambdas);
                ConnectBranches_MC_V0s(EParticle::KaonZeroShort, fPacked_MC_KaonsZeroShort);
            }
            break;
        case EReactionChannel::D:
            ConnectBranches_V0s(EParticle::AntiLambda, fPacked_AntiLambdas);
            ConnectBranches_Tracks(EParticle::PosKaon, fPacked_PosKaons);
            if (IsMC()) {
                ConnectBranches_MC_V0s(EParticle::AntiLambda, fPacked_MC_AntiLambdas);
                ConnectBranches_MC_Tracks(EParticle::PosKaon, fPacked_MC_PosKaons);
            }
            break;
        case EReactionChannel::E:
            ConnectBranches_V0s(EParticle::AntiLambda, fPacked_AntiLambdas);
            ConnectBranches_Tracks(EParticle::PosKaon, fPacked_PosKaons);
            ConnectBranches_Tracks(EParticle::PiMinus, fPacked_PiMinus);
            ConnectBranches_Tracks(EParticle::PiPlus, fPacked_PiPlus);
            if (IsMC()) {
                ConnectBranches_MC_V0s(EParticle::AntiLambda, fPacked_MC_AntiLambdas);
                ConnectBranches_MC_Tracks(EParticle::PosKaon, fPacked_MC_PosKaons);
                ConnectBranches_MC_Tracks(EParticle::PiMinus, fPacked_MC_PiMinus);
                ConnectBranches_MC_Tracks(EParticle::PiPlus, fPacked_MC_PiPlus);
            }
            break;
        case EReactionChannel::H:
            ConnectBranches_Tracks(EParticle::PosKaon, fPacked_PosKaons);
            if (IsMC()) ConnectBranches_MC_Tracks(EParticle::PosKaon, fPacked_MC_PosKaons);
            break;
        // anti-channels //
        case EReactionChannel::AntiA:
            ConnectBranches_V0s(EParticle::Lambda, fPacked_Lambdas);
            ConnectBranches_V0s(EParticle::KaonZeroShort, fPacked_KaonsZeroShort);
            if (IsMC()) {
                ConnectBranches_MC_V0s(EParticle::Lambda, fPacked_MC_Lambdas);
                ConnectBranches_MC_V0s(EParticle::KaonZeroShort, fPacked_MC_KaonsZeroShort);
            }
            break;
        case EReactionChannel::AntiD:
            ConnectBranches_V0s(EParticle::Lambda, fPacked_Lambdas);
            ConnectBranches_Tracks(EParticle::NegKaon, fPacked_NegKaons);
            if (IsMC()) {
                ConnectBranches_MC_V0s(EParticle::Lambda, fPacked_MC_Lambdas);
                ConnectBranches_MC_Tracks(EParticle::NegKaon, fPacked_MC_NegKaons);
            }
            break;
        case EReactionChannel::AntiE:
            ConnectBranches_V0s(EParticle::Lambda, fPacked_Lambdas);
            ConnectBranches_Tracks(EParticle::NegKaon, fPacked_NegKaons);
            ConnectBranches_Tracks(EParticle::PiMinus, fPacked_PiMinus);
            ConnectBranches_Tracks(EParticle::PiPlus, fPacked_PiPlus);
            if (IsMC()) {
                ConnectBranches_MC_V0s(EParticle::Lambda, fPacked_MC_Lambdas);
                ConnectBranches_MC_Tracks(EParticle::NegKaon, fPacked_MC_NegKaons);
                ConnectBranches_MC_Tracks(EParticle::PiMinus, fPacked_MC_PiMinus);
                ConnectBranches_MC_Tracks(EParticle::PiPlus, fPacked_MC_PiPlus);
            }
            break;
        case EReactionChannel::AntiH:
            ConnectBranches_Tracks(EParticle::NegKaon, fPacked_NegKaons);
            if (IsMC()) ConnectBranches_MC_Tracks(EParticle::NegKaon, fPacked_MC_NegKaons);
            break;
        default:
            break;
    }  // end of switch statement

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
}

void Finder::ConnectBranches_Events() {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    Utils::ConnectBranch(fTree_PackedEvents.get(), "RunNumber", &fInput_Event.RunNumber);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "DirNumber", &fInput_Event.DirNumber);
    if (!IsMC()) Utils::ConnectBranch(fTree_PackedEvents.get(), "DirNumberB", &fInput_Event.DirNumberB);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "EventNumber", &fInput_Event.EventNumber);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "Centrality", &fInput_Event.Centrality);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "MagneticField", &fInput_Event.MagneticField);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "PV_Xv", &fInput_Event.PV_Xv);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "PV_Yv", &fInput_Event.PV_Yv);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "PV_Zv", &fInput_Event.PV_Zv);

    if (IsMC()) {
        Utils::ConnectBranch(fTree_PackedEvents.get(), "MC_PV_Xv", &fInput_Event.MC_PV_Xv);
        Utils::ConnectBranch(fTree_PackedEvents.get(), "MC_PV_Yv", &fInput_Event.MC_PV_Yv);
        Utils::ConnectBranch(fTree_PackedEvents.get(), "MC_PV_Zv", &fInput_Event.MC_PV_Zv);
    }

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
}

void Finder::ConnectBranches_Injected() {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    Utils::ConnectBranch(fTree_PackedEvents.get(), "ReactionID", &fInput_Injected.ReactionID);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "SV_X", &fInput_Injected.X);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "SV_Y", &fInput_Injected.Y);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "SV_Z", &fInput_Injected.Z);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "Sexaquark_Px", &fInput_Injected.Px);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "Sexaquark_Py", &fInput_Injected.Py);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "Sexaquark_Pz", &fInput_Injected.Pz);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "Nucleon_Px", &fInput_Injected.Nucleon_Px);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "Nucleon_Py", &fInput_Injected.Nucleon_Py);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "Nucleon_Pz", &fInput_Injected.Nucleon_Pz);

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
}

void Finder::ConnectBranches_V0s(EParticle pid, PackedEvents::V0s& vec_v0s) {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Entry", Particle::Acronym[pid]), &vec_v0s.Entry);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_X", Particle::Acronym[pid]), &vec_v0s.X);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Y", Particle::Acronym[pid]), &vec_v0s.Y);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Z", Particle::Acronym[pid]), &vec_v0s.Z);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Px", Particle::Acronym[pid]), &vec_v0s.Px);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Py", Particle::Acronym[pid]), &vec_v0s.Py);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Pz", Particle::Acronym[pid]), &vec_v0s.Pz);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_E", Particle::Acronym[pid]), &vec_v0s.E);

    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaX2", Particle::Acronym[pid]), &vec_v0s.Sigma.X2);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaXY", Particle::Acronym[pid]), &vec_v0s.Sigma.XY);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaY2", Particle::Acronym[pid]), &vec_v0s.Sigma.Y2);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaXZ", Particle::Acronym[pid]), &vec_v0s.Sigma.XZ);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaYZ", Particle::Acronym[pid]), &vec_v0s.Sigma.YZ);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaZ2", Particle::Acronym[pid]), &vec_v0s.Sigma.Z2);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaXPx", Particle::Acronym[pid]), &vec_v0s.Sigma.XPx);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaYPx", Particle::Acronym[pid]), &vec_v0s.Sigma.YPx);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaZPx", Particle::Acronym[pid]), &vec_v0s.Sigma.ZPx);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaPx2", Particle::Acronym[pid]), &vec_v0s.Sigma.Px2);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaXPy", Particle::Acronym[pid]), &vec_v0s.Sigma.XPy);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaYPy", Particle::Acronym[pid]), &vec_v0s.Sigma.YPy);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaZPy", Particle::Acronym[pid]), &vec_v0s.Sigma.ZPy);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaPxPy", Particle::Acronym[pid]), &vec_v0s.Sigma.PxPy);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaPy2", Particle::Acronym[pid]), &vec_v0s.Sigma.Py2);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaXPz", Particle::Acronym[pid]), &vec_v0s.Sigma.XPz);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaYPz", Particle::Acronym[pid]), &vec_v0s.Sigma.YPz);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaZPz", Particle::Acronym[pid]), &vec_v0s.Sigma.ZPz);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaPxPz", Particle::Acronym[pid]), &vec_v0s.Sigma.PxPz);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaPyPz", Particle::Acronym[pid]), &vec_v0s.Sigma.PyPz);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaPz2", Particle::Acronym[pid]), &vec_v0s.Sigma.Pz2);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaXE", Particle::Acronym[pid]), &vec_v0s.Sigma.XE);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaYE", Particle::Acronym[pid]), &vec_v0s.Sigma.YE);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaZE", Particle::Acronym[pid]), &vec_v0s.Sigma.ZE);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaPxE", Particle::Acronym[pid]), &vec_v0s.Sigma.PxE);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaPyE", Particle::Acronym[pid]), &vec_v0s.Sigma.PyE);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaPzE", Particle::Acronym[pid]), &vec_v0s.Sigma.PzE);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaE2", Particle::Acronym[pid]), &vec_v0s.Sigma.E2);

    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Neg_Entry", Particle::Acronym[pid]), &vec_v0s.Neg.Entry);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Neg_X", Particle::Acronym[pid]), &vec_v0s.Neg.X);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Neg_Y", Particle::Acronym[pid]), &vec_v0s.Neg.Y);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Neg_Z", Particle::Acronym[pid]), &vec_v0s.Neg.Z);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Neg_Px", Particle::Acronym[pid]), &vec_v0s.Neg.Px);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Neg_Py", Particle::Acronym[pid]), &vec_v0s.Neg.Py);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Neg_Pz", Particle::Acronym[pid]), &vec_v0s.Neg.Pz);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Neg_E", Particle::Acronym[pid]), &vec_v0s.Neg.E);

    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Neg_X_AtPCA", Particle::Acronym[pid]), &vec_v0s.Neg_X_AtPCA);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Neg_Y_AtPCA", Particle::Acronym[pid]), &vec_v0s.Neg_Y_AtPCA);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Neg_Z_AtPCA", Particle::Acronym[pid]), &vec_v0s.Neg_Z_AtPCA);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Neg_Px_AtPCA", Particle::Acronym[pid]), &vec_v0s.Neg_Px_AtPCA);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Neg_Py_AtPCA", Particle::Acronym[pid]), &vec_v0s.Neg_Py_AtPCA);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Neg_Pz_AtPCA", Particle::Acronym[pid]), &vec_v0s.Neg_Pz_AtPCA);

    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Pos_Entry", Particle::Acronym[pid]), &vec_v0s.Pos.Entry);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Pos_X", Particle::Acronym[pid]), &vec_v0s.Pos.X);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Pos_Y", Particle::Acronym[pid]), &vec_v0s.Pos.Y);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Pos_Z", Particle::Acronym[pid]), &vec_v0s.Pos.Z);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Pos_Px", Particle::Acronym[pid]), &vec_v0s.Pos.Px);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Pos_Py", Particle::Acronym[pid]), &vec_v0s.Pos.Py);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Pos_Pz", Particle::Acronym[pid]), &vec_v0s.Pos.Pz);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Pos_E", Particle::Acronym[pid]), &vec_v0s.Pos.E);

    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Pos_X_AtPCA", Particle::Acronym[pid]), &vec_v0s.Pos_X_AtPCA);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Pos_Y_AtPCA", Particle::Acronym[pid]), &vec_v0s.Pos_Y_AtPCA);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Pos_Z_AtPCA", Particle::Acronym[pid]), &vec_v0s.Pos_Z_AtPCA);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Pos_Px_AtPCA", Particle::Acronym[pid]), &vec_v0s.Pos_Px_AtPCA);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Pos_Py_AtPCA", Particle::Acronym[pid]), &vec_v0s.Pos_Py_AtPCA);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Pos_Pz_AtPCA", Particle::Acronym[pid]), &vec_v0s.Pos_Pz_AtPCA);

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
}

void Finder::ConnectBranches_Tracks(EParticle pid, PackedEvents::Tracks& vec_tracks) {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Entry", Particle::Acronym[pid]), &vec_tracks.Entry);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_X", Particle::Acronym[pid]), &vec_tracks.X);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Y", Particle::Acronym[pid]), &vec_tracks.Y);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Z", Particle::Acronym[pid]), &vec_tracks.Z);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Px", Particle::Acronym[pid]), &vec_tracks.Px);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Py", Particle::Acronym[pid]), &vec_tracks.Py);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_Pz", Particle::Acronym[pid]), &vec_tracks.Pz);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_E", Particle::Acronym[pid]), &vec_tracks.E);

    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaX2", Particle::Acronym[pid]), &vec_tracks.Sigma.X2);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaXY", Particle::Acronym[pid]), &vec_tracks.Sigma.XY);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaY2", Particle::Acronym[pid]), &vec_tracks.Sigma.Y2);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaXZ", Particle::Acronym[pid]), &vec_tracks.Sigma.XZ);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaYZ", Particle::Acronym[pid]), &vec_tracks.Sigma.YZ);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaZ2", Particle::Acronym[pid]), &vec_tracks.Sigma.Z2);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaXPx", Particle::Acronym[pid]), &vec_tracks.Sigma.XPx);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaYPx", Particle::Acronym[pid]), &vec_tracks.Sigma.YPx);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaZPx", Particle::Acronym[pid]), &vec_tracks.Sigma.ZPx);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaPx2", Particle::Acronym[pid]), &vec_tracks.Sigma.Px2);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaXPy", Particle::Acronym[pid]), &vec_tracks.Sigma.XPy);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaYPy", Particle::Acronym[pid]), &vec_tracks.Sigma.YPy);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaZPy", Particle::Acronym[pid]), &vec_tracks.Sigma.ZPy);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaPxPy", Particle::Acronym[pid]), &vec_tracks.Sigma.PxPy);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaPy2", Particle::Acronym[pid]), &vec_tracks.Sigma.Py2);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaXPz", Particle::Acronym[pid]), &vec_tracks.Sigma.XPz);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaYPz", Particle::Acronym[pid]), &vec_tracks.Sigma.YPz);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaZPz", Particle::Acronym[pid]), &vec_tracks.Sigma.ZPz);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaPxPz", Particle::Acronym[pid]), &vec_tracks.Sigma.PxPz);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaPyPz", Particle::Acronym[pid]), &vec_tracks.Sigma.PyPz);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaPz2", Particle::Acronym[pid]), &vec_tracks.Sigma.Pz2);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaXE", Particle::Acronym[pid]), &vec_tracks.Sigma.XE);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaYE", Particle::Acronym[pid]), &vec_tracks.Sigma.YE);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaZE", Particle::Acronym[pid]), &vec_tracks.Sigma.ZE);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaPxE", Particle::Acronym[pid]), &vec_tracks.Sigma.PxE);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaPyE", Particle::Acronym[pid]), &vec_tracks.Sigma.PyE);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaPzE", Particle::Acronym[pid]), &vec_tracks.Sigma.PzE);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_SigmaE2", Particle::Acronym[pid]), &vec_tracks.Sigma.E2);

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
}

void Finder::ConnectBranches_MC_V0s(EParticle pid, PackedEvents::MC_V0s& sov) {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Entry", Particle::Acronym[pid]), &sov.Entry);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_X", Particle::Acronym[pid]), &sov.X);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Y", Particle::Acronym[pid]), &sov.Y);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Z", Particle::Acronym[pid]), &sov.Z);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Px", Particle::Acronym[pid]), &sov.Px);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Py", Particle::Acronym[pid]), &sov.Py);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Pz", Particle::Acronym[pid]), &sov.Pz);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_E", Particle::Acronym[pid]), &sov.E);

    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_DecayX", Particle::Acronym[pid]), &sov.DecayX);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_DecayY", Particle::Acronym[pid]), &sov.DecayY);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_DecayZ", Particle::Acronym[pid]), &sov.DecayZ);

    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_PdgCode", Particle::Acronym[pid]), &sov.PdgCode);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Mother_Entry", Particle::Acronym[pid]), &sov.Mother_Entry);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Mother_PdgCode", Particle::Acronym[pid]), &sov.Mother_PdgCode);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_IsTrue", Particle::Acronym[pid]), &sov.IsTrue);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_IsSignal", Particle::Acronym[pid]), &sov.IsSignal);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_IsSecondary", Particle::Acronym[pid]), &sov.IsSecondary);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_ReactionID", Particle::Acronym[pid]), &sov.ReactionID);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_IsHybrid", Particle::Acronym[pid]), &sov.IsHybrid);

    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Neg_Entry", Particle::Acronym[pid]), &sov.Neg_Entry);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Neg_Px", Particle::Acronym[pid]), &sov.Neg_Px);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Neg_Py", Particle::Acronym[pid]), &sov.Neg_Py);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Neg_Pz", Particle::Acronym[pid]), &sov.Neg_Pz);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Neg_PdgCode", Particle::Acronym[pid]), &sov.Neg_PdgCode);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Neg_IsTrue", Particle::Acronym[pid]), &sov.Neg_IsTrue);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Neg_IsSignal", Particle::Acronym[pid]), &sov.Neg_IsSignal);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Neg_IsSecondary", Particle::Acronym[pid]), &sov.Neg_IsSecondary);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Neg_ReactionID", Particle::Acronym[pid]), &sov.Neg_ReactionID);

    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Pos_Entry", Particle::Acronym[pid]), &sov.Pos_Entry);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Pos_Px", Particle::Acronym[pid]), &sov.Pos_Px);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Pos_Py", Particle::Acronym[pid]), &sov.Pos_Py);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Pos_Pz", Particle::Acronym[pid]), &sov.Pos_Pz);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Pos_PdgCode", Particle::Acronym[pid]), &sov.Pos_PdgCode);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Pos_IsTrue", Particle::Acronym[pid]), &sov.Pos_IsTrue);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Pos_IsSignal", Particle::Acronym[pid]), &sov.Pos_IsSignal);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Pos_IsSecondary", Particle::Acronym[pid]), &sov.Pos_IsSecondary);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Pos_ReactionID", Particle::Acronym[pid]), &sov.Pos_ReactionID);

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
}

void Finder::ConnectBranches_MC_Tracks(EParticle pid, PackedEvents::MC_Tracks& sov) {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Entry", Particle::Acronym[pid]), &sov.Entry);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_X", Particle::Acronym[pid]), &sov.X);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Y", Particle::Acronym[pid]), &sov.Y);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Z", Particle::Acronym[pid]), &sov.Z);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Px", Particle::Acronym[pid]), &sov.Px);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Py", Particle::Acronym[pid]), &sov.Py);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Pz", Particle::Acronym[pid]), &sov.Pz);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_E", Particle::Acronym[pid]), &sov.E);

    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_PdgCode", Particle::Acronym[pid]), &sov.PdgCode);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Mother_Entry", Particle::Acronym[pid]), &sov.Mother_Entry);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_Mother_PdgCode", Particle::Acronym[pid]), &sov.Mother_PdgCode);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_GrandMother_Entry", Particle::Acronym[pid]), &sov.GrandMother_Entry);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_GrandMother_PdgCode", Particle::Acronym[pid]), &sov.GrandMother_PdgCode);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_IsTrue", Particle::Acronym[pid]), &sov.IsTrue);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_IsSignal", Particle::Acronym[pid]), &sov.IsSignal);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_IsSecondary", Particle::Acronym[pid]), &sov.IsSecondary);
    Utils::ConnectBranch(fTree_PackedEvents.get(), fmt::format("{}_MC_ReactionID", Particle::Acronym[pid]), &sov.ReactionID);

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
}

// ## OUTPUT ZONE ## //

bool Finder::PrepareOutputFile() {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    const std::filesystem::path output_path(fSettings.PathOutputFile);
    if (output_path.has_parent_path()) std::filesystem::create_directories(output_path.parent_path());

    fOutputFile = std::unique_ptr<TFile>(TFile::Open(fSettings.PathOutputFile.c_str(), "RECREATE"));
    if (!fOutputFile) {
        Logger::Error(__FUNCTION__, "Couldn't create TFile {}", fSettings.PathOutputFile);
        return false;
    }

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
    return true;
}

bool Finder::PrepareOutputTree() {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    std::string tree_name{fmt::format("Candidates_{}", Name::ReactionChannel[GetReactionChannel()])};

    fOutputTree = std::make_unique<TTree>(tree_name.c_str(), "");
    if (!fOutputTree) {
        Logger::Error(__FUNCTION__, "Couldn't create TTree \"{}\"");
        return false;
    }

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
    return true;
}

void Finder::CreateOutputBranches(Found::ChannelA& out_branches) {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    Utils::CreateBranch(fOutputTree.get(), "RunNumber", &out_branches.RunNumber);
    Utils::CreateBranch(fOutputTree.get(), "DirNumber", &out_branches.DirNumber);
    if (!IsMC()) Utils::CreateBranch(fOutputTree.get(), "DirNumberB", &out_branches.DirNumberB);
    Utils::CreateBranch(fOutputTree.get(), "EventNumber", &out_branches.EventNumber);
    Utils::CreateBranch(fOutputTree.get(), "MagneticField", &out_branches.MagneticField);
    Utils::CreateBranch(fOutputTree.get(), "PV_Xv", &out_branches.PV_Xv);
    Utils::CreateBranch(fOutputTree.get(), "PV_Yv", &out_branches.PV_Yv);
    Utils::CreateBranch(fOutputTree.get(), "PV_Zv", &out_branches.PV_Zv);
    if (IsMC()) {
        Utils::CreateBranch(fOutputTree.get(), "MC_PV_Xv", &out_branches.MC_PV_Xv);
        Utils::CreateBranch(fOutputTree.get(), "MC_PV_Yv", &out_branches.MC_PV_Yv);
        Utils::CreateBranch(fOutputTree.get(), "MC_PV_Zv", &out_branches.MC_PV_Zv);
    }

    Utils::CreateBranch(fOutputTree.get(), "X", &out_branches.X);
    Utils::CreateBranch(fOutputTree.get(), "Y", &out_branches.Y);
    Utils::CreateBranch(fOutputTree.get(), "Z", &out_branches.Z);
    Utils::CreateBranch(fOutputTree.get(), "Px", &out_branches.Px);
    Utils::CreateBranch(fOutputTree.get(), "Py", &out_branches.Py);
    Utils::CreateBranch(fOutputTree.get(), "Pz", &out_branches.Pz);
    Utils::CreateBranch(fOutputTree.get(), "E", &out_branches.E);
    Utils::CreateBranch(fOutputTree.get(), "E_MinusNucleon", &out_branches.E_MinusNucleon);

    Utils::CreateBranch(fOutputTree.get(), "V0A_X_AtPCA", &out_branches.V0A_X_AtPCA);
    Utils::CreateBranch(fOutputTree.get(), "V0A_Y_AtPCA", &out_branches.V0A_Y_AtPCA);
    Utils::CreateBranch(fOutputTree.get(), "V0A_Z_AtPCA", &out_branches.V0A_Z_AtPCA);

    Utils::CreateBranch(fOutputTree.get(), "V0B_X_AtPCA", &out_branches.V0B_X_AtPCA);
    Utils::CreateBranch(fOutputTree.get(), "V0B_Y_AtPCA", &out_branches.V0B_Y_AtPCA);
    Utils::CreateBranch(fOutputTree.get(), "V0B_Z_AtPCA", &out_branches.V0B_Z_AtPCA);

    Utils::CreateBranch(fOutputTree.get(), "V0A_Entry", &out_branches.V0A_Entry);
    Utils::CreateBranch(fOutputTree.get(), "V0A_X_AtDecay", &out_branches.V0A.X);
    Utils::CreateBranch(fOutputTree.get(), "V0A_Y_AtDecay", &out_branches.V0A.Y);
    Utils::CreateBranch(fOutputTree.get(), "V0A_Z_AtDecay", &out_branches.V0A.Z);
    Utils::CreateBranch(fOutputTree.get(), "V0A_Px", &out_branches.V0A.Px);
    Utils::CreateBranch(fOutputTree.get(), "V0A_Py", &out_branches.V0A.Py);
    Utils::CreateBranch(fOutputTree.get(), "V0A_Pz", &out_branches.V0A.Pz);
    Utils::CreateBranch(fOutputTree.get(), "V0A_E", &out_branches.V0A.E);

    Utils::CreateBranch(fOutputTree.get(), "V0A_Neg_Entry", &out_branches.V0A_Neg_Entry);
    Utils::CreateBranch(fOutputTree.get(), "V0A_Pos_Entry", &out_branches.V0A_Pos_Entry);

    Utils::CreateBranch(fOutputTree.get(), "V0B_Entry", &out_branches.V0B_Entry);
    Utils::CreateBranch(fOutputTree.get(), "V0B_X_AtDecay", &out_branches.V0B.X);
    Utils::CreateBranch(fOutputTree.get(), "V0B_Y_AtDecay", &out_branches.V0B.Y);
    Utils::CreateBranch(fOutputTree.get(), "V0B_Z_AtDecay", &out_branches.V0B.Z);
    Utils::CreateBranch(fOutputTree.get(), "V0B_Px", &out_branches.V0B.Px);
    Utils::CreateBranch(fOutputTree.get(), "V0B_Py", &out_branches.V0B.Py);
    Utils::CreateBranch(fOutputTree.get(), "V0B_Pz", &out_branches.V0B.Pz);
    Utils::CreateBranch(fOutputTree.get(), "V0B_E", &out_branches.V0B.E);

    Utils::CreateBranch(fOutputTree.get(), "V0B_Neg_Entry", &out_branches.V0B_Neg_Entry);
    Utils::CreateBranch(fOutputTree.get(), "V0B_Pos_Entry", &out_branches.V0B_Pos_Entry);

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
}

void Finder::CreateOutputBranches(Found::ChannelD& out_branches) {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    Utils::CreateBranch(fOutputTree.get(), "RunNumber", &out_branches.RunNumber);
    Utils::CreateBranch(fOutputTree.get(), "DirNumber", &out_branches.DirNumber);
    if (!IsMC()) Utils::CreateBranch(fOutputTree.get(), "DirNumberB", &out_branches.DirNumberB);
    Utils::CreateBranch(fOutputTree.get(), "EventNumber", &out_branches.EventNumber);
    Utils::CreateBranch(fOutputTree.get(), "MagneticField", &out_branches.MagneticField);
    Utils::CreateBranch(fOutputTree.get(), "PV_Xv", &out_branches.PV_Xv);
    Utils::CreateBranch(fOutputTree.get(), "PV_Yv", &out_branches.PV_Yv);
    Utils::CreateBranch(fOutputTree.get(), "PV_Zv", &out_branches.PV_Zv);
    if (IsMC()) {
        Utils::CreateBranch(fOutputTree.get(), "MC_PV_Xv", &out_branches.MC_PV_Xv);
        Utils::CreateBranch(fOutputTree.get(), "MC_PV_Yv", &out_branches.MC_PV_Yv);
        Utils::CreateBranch(fOutputTree.get(), "MC_PV_Zv", &out_branches.MC_PV_Zv);
    }

    Utils::CreateBranch(fOutputTree.get(), "X", &out_branches.X);
    Utils::CreateBranch(fOutputTree.get(), "Y", &out_branches.Y);
    Utils::CreateBranch(fOutputTree.get(), "Z", &out_branches.Z);
    Utils::CreateBranch(fOutputTree.get(), "Px", &out_branches.Px);
    Utils::CreateBranch(fOutputTree.get(), "Py", &out_branches.Py);
    Utils::CreateBranch(fOutputTree.get(), "Pz", &out_branches.Pz);
    Utils::CreateBranch(fOutputTree.get(), "E", &out_branches.E);
    Utils::CreateBranch(fOutputTree.get(), "E_MinusNucleon", &out_branches.E_MinusNucleon);

    Utils::CreateBranch(fOutputTree.get(), "V0_X_AtPCA", &out_branches.V0_X_AtPCA);
    Utils::CreateBranch(fOutputTree.get(), "V0_Y_AtPCA", &out_branches.V0_Y_AtPCA);
    Utils::CreateBranch(fOutputTree.get(), "V0_Z_AtPCA", &out_branches.V0_Z_AtPCA);

    Utils::CreateBranch(fOutputTree.get(), "Kaon_X_AtPCA", &out_branches.Kaon_X_AtPCA);
    Utils::CreateBranch(fOutputTree.get(), "Kaon_Y_AtPCA", &out_branches.Kaon_Y_AtPCA);
    Utils::CreateBranch(fOutputTree.get(), "Kaon_Z_AtPCA", &out_branches.Kaon_Z_AtPCA);
    Utils::CreateBranch(fOutputTree.get(), "Kaon_Px_AtPCA", &out_branches.Kaon_Px_AtPCA);
    Utils::CreateBranch(fOutputTree.get(), "Kaon_Py_AtPCA", &out_branches.Kaon_Py_AtPCA);
    Utils::CreateBranch(fOutputTree.get(), "Kaon_Pz_AtPCA", &out_branches.Kaon_Pz_AtPCA);

    Utils::CreateBranch(fOutputTree.get(), "V0_Entry", &out_branches.V0_Entry);
    Utils::CreateBranch(fOutputTree.get(), "V0_X_AtDecay", &out_branches.V0.X);
    Utils::CreateBranch(fOutputTree.get(), "V0_Y_AtDecay", &out_branches.V0.Y);
    Utils::CreateBranch(fOutputTree.get(), "V0_Z_AtDecay", &out_branches.V0.Z);
    Utils::CreateBranch(fOutputTree.get(), "V0_Px", &out_branches.V0.Px);
    Utils::CreateBranch(fOutputTree.get(), "V0_Py", &out_branches.V0.Py);
    Utils::CreateBranch(fOutputTree.get(), "V0_Pz", &out_branches.V0.Pz);
    Utils::CreateBranch(fOutputTree.get(), "V0_E", &out_branches.V0.E);

    Utils::CreateBranch(fOutputTree.get(), "V0_Neg_Entry", &out_branches.V0_Neg_Entry);
    Utils::CreateBranch(fOutputTree.get(), "V0_Pos_Entry", &out_branches.V0_Pos_Entry);

    Utils::CreateBranch(fOutputTree.get(), "Kaon_Entry", &out_branches.Kaon_Entry);
    Utils::CreateBranch(fOutputTree.get(), "Kaon_X", &out_branches.Kaon.X);
    Utils::CreateBranch(fOutputTree.get(), "Kaon_Y", &out_branches.Kaon.Y);
    Utils::CreateBranch(fOutputTree.get(), "Kaon_Z", &out_branches.Kaon.Z);
    Utils::CreateBranch(fOutputTree.get(), "Kaon_Px", &out_branches.Kaon.Px);
    Utils::CreateBranch(fOutputTree.get(), "Kaon_Py", &out_branches.Kaon.Py);
    Utils::CreateBranch(fOutputTree.get(), "Kaon_Pz", &out_branches.Kaon.Pz);
    Utils::CreateBranch(fOutputTree.get(), "Kaon_E", &out_branches.Kaon.E);

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
}

void Finder::CreateOutputBranches(Found::ChannelE& out_branches) {
    // PENDING
}

void Finder::CreateOutputBranches(Found::ChannelH& out_branches) {
    // PENDING
}

void Finder::CreateOutputBranches(Found::MC_ChannelA& sov) {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    Utils::CreateBranch(fOutputTree.get(), "MC_X", &sov.X);
    Utils::CreateBranch(fOutputTree.get(), "MC_Y", &sov.Y);
    Utils::CreateBranch(fOutputTree.get(), "MC_Z", &sov.Z);
    Utils::CreateBranch(fOutputTree.get(), "MC_Px_Before", &sov.BeforePx);
    Utils::CreateBranch(fOutputTree.get(), "MC_Py_Before", &sov.BeforePy);
    Utils::CreateBranch(fOutputTree.get(), "MC_Pz_Before", &sov.BeforePz);
    Utils::CreateBranch(fOutputTree.get(), "MC_E_Before", &sov.BeforeE);

    Utils::CreateBranch(fOutputTree.get(), "MC_Px_After", &sov.AfterPx);
    Utils::CreateBranch(fOutputTree.get(), "MC_Py_After", &sov.AfterPy);
    Utils::CreateBranch(fOutputTree.get(), "MC_Pz_After", &sov.AfterPz);
    Utils::CreateBranch(fOutputTree.get(), "MC_E_After", &sov.AfterE);

    Utils::CreateBranch(fOutputTree.get(), "MC_IsSignal", &sov.IsSignal);
    Utils::CreateBranch(fOutputTree.get(), "MC_ReactionID", &sov.ReactionID);
    Utils::CreateBranch(fOutputTree.get(), "MC_IsHybrid", &sov.IsHybrid);

    Utils::CreateBranch(fOutputTree.get(), "MC_Px_Nucleon", &sov.NucleonPx);
    Utils::CreateBranch(fOutputTree.get(), "MC_Py_Nucleon", &sov.NucleonPy);
    Utils::CreateBranch(fOutputTree.get(), "MC_Pz_Nucleon", &sov.NucleonPz);
    Utils::CreateBranch(fOutputTree.get(), "MC_E_Nucleon", &sov.NucleonPz);

    Utils::CreateBranch(fOutputTree.get(), "MC_V0A_Entry", &sov.V0A_Entry);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0A_PdgCode", &sov.V0A_PdgCode);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0A_Mother_Entry", &sov.V0A_Mother_Entry);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0A_Mother_PdgCode", &sov.V0A_Mother_PdgCode);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0A_Neg_Entry", &sov.V0A_Neg_Entry);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0A_Pos_Entry", &sov.V0A_Pos_Entry);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0A_Px", &sov.V0A_Px);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0A_Py", &sov.V0A_Py);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0A_Pz", &sov.V0A_Pz);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0A_E", &sov.V0A_E);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0A_IsTrue", &sov.V0A_IsTrue);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0A_IsSignal", &sov.V0A_IsSignal);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0A_IsSecondary", &sov.V0A_IsSecondary);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0A_ReactionID", &sov.V0A_ReactionID);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0A_IsHybrid", &sov.V0A_IsHybrid);

    Utils::CreateBranch(fOutputTree.get(), "MC_V0B_Entry", &sov.V0B_Entry);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0B_PdgCode", &sov.V0B_PdgCode);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0B_Mother_Entry", &sov.V0B_Mother_Entry);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0B_Mother_PdgCode", &sov.V0B_Mother_PdgCode);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0B_Neg_Entry", &sov.V0B_Neg_Entry);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0B_Pos_Entry", &sov.V0B_Pos_Entry);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0B_Px", &sov.V0B_Px);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0B_Py", &sov.V0B_Py);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0B_Pz", &sov.V0B_Pz);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0B_E", &sov.V0B_E);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0B_IsTrue", &sov.V0B_IsTrue);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0B_IsSignal", &sov.V0B_IsSignal);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0B_IsSecondary", &sov.V0B_IsSecondary);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0B_ReactionID", &sov.V0B_ReactionID);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0B_IsHybrid", &sov.V0B_IsHybrid);

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
}

void Finder::CreateOutputBranches(Found::MC_ChannelD& sov) {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    Utils::CreateBranch(fOutputTree.get(), "MC_X", &sov.X);
    Utils::CreateBranch(fOutputTree.get(), "MC_Y", &sov.Y);
    Utils::CreateBranch(fOutputTree.get(), "MC_Z", &sov.Z);
    Utils::CreateBranch(fOutputTree.get(), "MC_Px_Before", &sov.BeforePx);
    Utils::CreateBranch(fOutputTree.get(), "MC_Py_Before", &sov.BeforePy);
    Utils::CreateBranch(fOutputTree.get(), "MC_Pz_Before", &sov.BeforePz);
    Utils::CreateBranch(fOutputTree.get(), "MC_E_Before", &sov.BeforeE);

    Utils::CreateBranch(fOutputTree.get(), "MC_Px_After", &sov.AfterPx);
    Utils::CreateBranch(fOutputTree.get(), "MC_Py_After", &sov.AfterPy);
    Utils::CreateBranch(fOutputTree.get(), "MC_Pz_After", &sov.AfterPz);
    Utils::CreateBranch(fOutputTree.get(), "MC_E_After", &sov.AfterE);

    Utils::CreateBranch(fOutputTree.get(), "MC_IsSignal", &sov.IsSignal);
    Utils::CreateBranch(fOutputTree.get(), "MC_ReactionID", &sov.ReactionID);
    Utils::CreateBranch(fOutputTree.get(), "MC_IsHybrid", &sov.IsHybrid);

    Utils::CreateBranch(fOutputTree.get(), "MC_Nucleon_Px", &sov.NucleonPx);
    Utils::CreateBranch(fOutputTree.get(), "MC_Nucleon_Py", &sov.NucleonPy);
    Utils::CreateBranch(fOutputTree.get(), "MC_Nucleon_Pz", &sov.NucleonPz);
    Utils::CreateBranch(fOutputTree.get(), "MC_Nucleon_E", &sov.NucleonPz);

    Utils::CreateBranch(fOutputTree.get(), "MC_V0_Entry", &sov.V0_Entry);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0_PdgCode", &sov.V0_PdgCode);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0_Mother_Entry", &sov.V0_Mother_Entry);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0_Mother_PdgCode", &sov.V0_Mother_PdgCode);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0_Neg_Entry", &sov.V0_Neg_Entry);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0_Pos_Entry", &sov.V0_Pos_Entry);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0_Px", &sov.V0_Px);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0_Py", &sov.V0_Py);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0_Pz", &sov.V0_Pz);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0_E", &sov.V0_E);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0_IsTrue", &sov.V0_IsTrue);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0_IsSignal", &sov.V0_IsSignal);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0_IsSecondary", &sov.V0_IsSecondary);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0_ReactionID", &sov.V0_ReactionID);
    Utils::CreateBranch(fOutputTree.get(), "MC_V0_IsHybrid", &sov.V0_IsHybrid);

    Utils::CreateBranch(fOutputTree.get(), "MC_Kaon_Entry", &sov.Kaon_Entry);
    Utils::CreateBranch(fOutputTree.get(), "MC_Kaon_PdgCode", &sov.Kaon_PdgCode);
    Utils::CreateBranch(fOutputTree.get(), "MC_Kaon_Mother_Entry", &sov.Kaon_Mother_Entry);
    Utils::CreateBranch(fOutputTree.get(), "MC_Kaon_Mother_PdgCode", &sov.Kaon_Mother_PdgCode);
    Utils::CreateBranch(fOutputTree.get(), "MC_Kaon_GrandMother_Entry", &sov.Kaon_GrandMother_Entry);
    Utils::CreateBranch(fOutputTree.get(), "MC_Kaon_GrandMother_PdgCode", &sov.Kaon_GrandMother_PdgCode);
    Utils::CreateBranch(fOutputTree.get(), "MC_Kaon_Px", &sov.Kaon_Px);
    Utils::CreateBranch(fOutputTree.get(), "MC_Kaon_Py", &sov.Kaon_Py);
    Utils::CreateBranch(fOutputTree.get(), "MC_Kaon_Pz", &sov.Kaon_Pz);
    Utils::CreateBranch(fOutputTree.get(), "MC_Kaon_E", &sov.Kaon_E);
    Utils::CreateBranch(fOutputTree.get(), "MC_Kaon_IsTrue", &sov.Kaon_IsTrue);
    Utils::CreateBranch(fOutputTree.get(), "MC_Kaon_IsSignal", &sov.Kaon_IsSignal);
    Utils::CreateBranch(fOutputTree.get(), "MC_Kaon_IsSecondary", &sov.Kaon_IsSecondary);
    Utils::CreateBranch(fOutputTree.get(), "MC_Kaon_ReactionID", &sov.Kaon_ReactionID);

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
}

void Finder::CreateCutFlowHistogram() {
    const int x_nbins{20};
    const float x_min{0.};
    const float x_max{20.};
    std::string hist_title{";Cut N;N Passed Cut"};
    fCutFlowHist = std::make_unique<TH1D>("CutFlow", hist_title.c_str(), x_nbins, x_min, x_max);
}

// ## Injected ZONE ## //

bool Finder::Injected_PrepareOutputTree() {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    std::string tree_name{"Injected"};

    fOutputTree_Injected = std::make_unique<TTree>(tree_name.c_str(), "");
    if (!fOutputTree_Injected) {
        Logger::Error(__FUNCTION__, "Couldn't create TTree \"{}\"");
        return false;
    }

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
    return true;
}

void Finder::Injected_CreateOutputBranches() {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    fOutputTree_Injected->Branch("RunNumber", &fOutput_Injected.RunNumber);
    fOutputTree_Injected->Branch("DirNumber", &fOutput_Injected.DirNumber);
    fOutputTree_Injected->Branch("EventNumber", &fOutput_Injected.EventNumber);
    fOutputTree_Injected->Branch("ReactionID", &fOutput_Injected.ReactionID);
    fOutputTree_Injected->Branch("X", &fOutput_Injected.X);
    fOutputTree_Injected->Branch("Y", &fOutput_Injected.Y);
    fOutputTree_Injected->Branch("Z", &fOutput_Injected.Z);
    fOutputTree_Injected->Branch("Px", &fOutput_Injected.Px);
    fOutputTree_Injected->Branch("Py", &fOutput_Injected.Py);
    fOutputTree_Injected->Branch("Pz", &fOutput_Injected.Pz);
    fOutputTree_Injected->Branch("E", &fOutput_Injected.E);

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
};

void Finder::Injected_FlattenAndStore() {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    auto n_injected = static_cast<int>(fInput_Injected.ReactionID->size());
    for (int idx_inj{0}; idx_inj < n_injected; ++idx_inj) {
        fOutput_Injected.RunNumber = fInput_Event.RunNumber;
        fOutput_Injected.DirNumber = fInput_Event.DirNumber;
        fOutput_Injected.EventNumber = fInput_Event.EventNumber;
        fOutput_Injected.ReactionID = fInput_Injected.ReactionID->at(idx_inj);
        fOutput_Injected.X = fInput_Injected.X->at(idx_inj);
        fOutput_Injected.Y = fInput_Injected.Y->at(idx_inj);
        fOutput_Injected.Z = fInput_Injected.Z->at(idx_inj);
        fOutput_Injected.Px = fInput_Injected.Px->at(idx_inj);
        fOutput_Injected.Py = fInput_Injected.Py->at(idx_inj);
        fOutput_Injected.Pz = fInput_Injected.Pz->at(idx_inj);
        fOutput_Injected.E = static_cast<float>(std::sqrt(fSettings.SexaquarkMass * fSettings.SexaquarkMass +
                                                          static_cast<double>(fOutput_Injected.Px) * static_cast<double>(fOutput_Injected.Px) +
                                                          static_cast<double>(fOutput_Injected.Py) * static_cast<double>(fOutput_Injected.Py) +
                                                          static_cast<double>(fOutput_Injected.Pz) * static_cast<double>(fOutput_Injected.Pz)));
        fOutputTree_Injected->Fill();
    }

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
}

// ## Channel A ZONE ## //

void Finder::FindSexaquarks_ChannelA(bool anti_channel) {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    const PackedEvents::V0s* Packed_Lambdas{anti_channel ? &fPacked_Lambdas : &fPacked_AntiLambdas};
    const PackedEvents::MC_V0s* MC_Lambdas{anti_channel ? &fPacked_MC_Lambdas : &fPacked_MC_AntiLambdas};
    EParticle PID_Lambda{anti_channel ? EParticle::Lambda : EParticle::AntiLambda};

    // loop over all possible pairs of (anti)lambda + K0S //
    auto n_lambdas = static_cast<int>(Packed_Lambdas->Entry->size());
    auto n_k0s = static_cast<int>(fPacked_KaonsZeroShort.Entry->size());
    for (int idx_v0a{0}; idx_v0a < n_lambdas; ++idx_v0a) {

        // unpack (anti)lambda //
        auto v0a = KF::UnpackV0(*Packed_Lambdas, idx_v0a, PID_Lambda);

        for (int idx_v0b{0}; idx_v0b < n_k0s; ++idx_v0b) {

            // unpack K0S //
            auto v0b = KF::UnpackV0(fPacked_KaonsZeroShort, idx_v0b, EParticle::KaonZeroShort);

            // sanity check //
            std::set<int> unique_track_entries{v0a.Neg.idx, v0a.Pos.idx, v0b.Neg.idx, v0b.Pos.idx};
            if (unique_track_entries.size() < 4) continue;

            // fit //
            KF::ChannelA sexa{v0a, v0b};
            sexa.DoFit(fInput_Event.MagneticField);

#ifdef T2S_DEBUG
            Logger::Debug(__FUNCTION__, "idx(v0a,neg,pos)={},{},{}", v0a.idx, v0a.Neg.idx, v0a.Pos.idx);
            Logger::Debug(__FUNCTION__, ";x,y,z={},{},{}", v0a.X(), v0a.Y(), v0a.Z());
            Logger::Debug(__FUNCTION__, ";px,py,pz={},{},{}", v0a.Px(), v0a.Py(), v0a.Pz());
            Logger::Debug(__FUNCTION__, ";mass={}", v0a.Mass());

            Logger::Debug(__FUNCTION__, "idx(v0b,neg,pos)={},{},{}", v0b.idx, v0b.Neg.idx, v0b.Pos.idx);
            Logger::Debug(__FUNCTION__, ";x,y,z={},{},{}", v0b.X(), v0b.Y(), v0b.Z());
            Logger::Debug(__FUNCTION__, ";px,py,pz={},{},{}", v0b.Px(), v0b.Py(), v0b.Pz());
            Logger::Debug(__FUNCTION__, ";mass={}", v0b.Mass());

            Logger::Debug(__FUNCTION__, "x,y,z={},{},{}", sexa.X(), sexa.Y(), sexa.Z());
            Logger::Debug(__FUNCTION__, ";x,y,z(v0a)={},{},{}", sexa.V0A_PCA_XYZ()[0], sexa.V0A_PCA_XYZ()[1], sexa.V0A_PCA_XYZ()[2]);
            Logger::Debug(__FUNCTION__, ";x,y,z(v0b)={},{},{}", sexa.V0B_PCA_XYZ()[0], sexa.V0B_PCA_XYZ()[1], sexa.V0B_PCA_XYZ()[2]);
            Logger::Debug(__FUNCTION__, ";mass={}", sexa.Mass());
            Logger::Debug(__FUNCTION__, ";mass_minus_n={}", sexa.Mass_MinusNucleon());
            Logger::Debug(__FUNCTION__, ";dca_btw_v0s={}", sexa.DCA_btw_V0s());
            Logger::Debug(__FUNCTION__, ";radius={}", sexa.Radius2D());
            Logger::Debug(__FUNCTION__, ";dca_v0a={}", sexa.DCA_V0A_wrt_SV());
            Logger::Debug(__FUNCTION__, ";dca_v0b={}", sexa.DCA_V0B_wrt_SV());
            Logger::Debug(__FUNCTION__, ";dca_v0a_neg={}", sexa.DCA_V0ANeg_wrt_SV(fInput_Event.MagneticField));
            Logger::Debug(__FUNCTION__, ";dca_v0a_pos={}", sexa.DCA_V0APos_wrt_SV(fInput_Event.MagneticField));
            Logger::Debug(__FUNCTION__, ";dca_v0b_neg={}", sexa.DCA_V0BNeg_wrt_SV(fInput_Event.MagneticField));
            Logger::Debug(__FUNCTION__, ";dca_v0b_pos={}", sexa.DCA_V0BPos_wrt_SV(fInput_Event.MagneticField));
            Logger::Debug(__FUNCTION__, ";pt={}", sexa.Pt());
            Logger::Debug(__FUNCTION__, ";eta={}", sexa.Eta());
            Logger::Debug(__FUNCTION__, ";decay_length(v0a,v0b)={},{}", sexa.DecayLength_V0A(), sexa.DecayLength_V0B());
            Logger::Debug(__FUNCTION__, ";cpa_pv={}", sexa.CPA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv));

            MC::V0 mc_v0a_debug{*MC_Lambdas, v0a.idx};
            MC::V0 mc_v0b_debug{fPacked_MC_KaonsZeroShort, v0b.idx};
            MC::ChannelA mc_sexa_debug{fInput_Injected, fSettings.SexaquarkMass, mc_v0a_debug, mc_v0b_debug};
            if (IsMC()) {
                Logger::Debug(__FUNCTION__, "sexa(is_signal,reaction_id,is_hybrid)={},{},{}", mc_sexa_debug.IsSignal, mc_sexa_debug.ReactionID,
                              mc_sexa_debug.IsHybrid);
                Logger::Debug(__FUNCTION__, ";v0a(entry,pdg_code)={},{}", mc_sexa_debug.V0A.Entry, mc_sexa_debug.V0A.PdgCode);
                Logger::Debug(__FUNCTION__, ";v0a(is_signal,reaction_id,is_hybrid)={},{},{}", mc_sexa_debug.V0A.IsSignal,
                              mc_sexa_debug.V0A.ReactionID, mc_sexa_debug.V0A.IsHybrid);
                Logger::Debug(__FUNCTION__, ";v0b(entry,pdg_code)={},{}", mc_sexa_debug.V0B.Entry, mc_sexa_debug.V0B.PdgCode);
                Logger::Debug(__FUNCTION__, ";v0b(is_signal,reaction_id,is_hybrid)={},{},{}", mc_sexa_debug.V0B.IsSignal,
                              mc_sexa_debug.V0B.ReactionID, mc_sexa_debug.V0B.IsHybrid);
            }
#endif

            // apply cuts //
            if (!PassesCuts(sexa)) continue;

            // store //
            Store(sexa);
            if (IsMC()) {
                MC::V0 mc_v0a{*MC_Lambdas, v0a.idx};
                MC::V0 mc_v0b{fPacked_MC_KaonsZeroShort, v0b.idx};
                MC::ChannelA mc_sexa{fInput_Injected, fSettings.SexaquarkMass, mc_v0a, mc_v0b};
                StoreMC(mc_sexa);
            }
            fOutputTree->Fill();
        }
    }

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
}

bool Finder::PassesCuts(const KF::ChannelA& sexa) const {

    fCutFlowHist->Fill(0.);
    if (sexa.Radius2D() < Cuts::ChannelA::Min_Radius2D || sexa.Radius2D() > Cuts::ChannelA::Max_Radius2D) return false;
    fCutFlowHist->Fill(1.);
    if (sexa.DecayLength_V0A() > Cuts::ChannelA::Max_DecayLengthLa) return false;
    fCutFlowHist->Fill(2.);
    if (sexa.DecayLength_V0B() > Cuts::ChannelA::Max_DecayLengthK0) return false;
    fCutFlowHist->Fill(3.);
    if (sexa.AbsRapidity_MinusNucleon() > Cuts::ChannelA::AbsMax_Rapidity) return false;  // PENDING: kinematic cut, affected by Fermi motion
    fCutFlowHist->Fill(4.);
    if (sexa.Mass_MinusNucleon() < Cuts::ChannelA::Min_MassMinusNucleon || sexa.Mass_MinusNucleon() > Cuts::ChannelA::Max_MassMinusNucleon) {
        return false;  // PENDING: kinematic cut, affected by Fermi motion
    }
    fCutFlowHist->Fill(5.);
    if (sexa.CPA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv) < Cuts::ChannelA::Min_CPAwrtPV ||
        sexa.CPA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv) > Cuts::ChannelA::Max_CPAwrtPV) {
        return false;  // PENDING: kinematic cut, affected by Fermi motion
    }
    fCutFlowHist->Fill(6.);
    if (sexa.DCA_V0ANeg_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelA::Max_DCALaNegSV) return false;
    fCutFlowHist->Fill(7.);
    if (sexa.DCA_V0APos_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelA::Max_DCALaPosSV) return false;
    fCutFlowHist->Fill(8.);
    if (sexa.DCA_V0BNeg_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelA::Max_DCAK0NegSV) return false;
    fCutFlowHist->Fill(9.);
    if (sexa.DCA_V0BPos_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelA::Max_DCAK0PosSV) return false;
    fCutFlowHist->Fill(10.);
    if (sexa.DCA_V0A_wrt_SV() > Cuts::ChannelA::Max_DCALaSV) return false;
    fCutFlowHist->Fill(11.);
    if (sexa.DCA_V0B_wrt_SV() > Cuts::ChannelA::Max_DCAK0SV) return false;
    fCutFlowHist->Fill(12.);
    if (sexa.DCA_btw_V0s() > Cuts::ChannelA::Max_DCAbtwV0s) return false;
    fCutFlowHist->Fill(13.);

    return true;
}

void Finder::Store(const KF::ChannelA& sexa) {

    fOutput_ChannelA.RunNumber = fInput_Event.RunNumber;
    fOutput_ChannelA.DirNumber = fInput_Event.DirNumber;
    if (!IsMC()) fOutput_ChannelA.DirNumberB = fInput_Event.DirNumberB;
    fOutput_ChannelA.EventNumber = fInput_Event.EventNumber;
    fOutput_ChannelA.MagneticField = fInput_Event.MagneticField;
    fOutput_ChannelA.PV_Xv = fInput_Event.PV_Xv;
    fOutput_ChannelA.PV_Yv = fInput_Event.PV_Yv;
    fOutput_ChannelA.PV_Zv = fInput_Event.PV_Zv;
    if (IsMC()) {
        fOutput_ChannelA.MC_PV_Xv = fInput_Event.MC_PV_Xv;
        fOutput_ChannelA.MC_PV_Yv = fInput_Event.MC_PV_Yv;
        fOutput_ChannelA.MC_PV_Zv = fInput_Event.MC_PV_Zv;
    }

    fOutput_ChannelA.X = static_cast<float>(sexa.X());
    fOutput_ChannelA.Y = static_cast<float>(sexa.Y());
    fOutput_ChannelA.Z = static_cast<float>(sexa.Z());
    fOutput_ChannelA.Px = static_cast<float>(sexa.Px());
    fOutput_ChannelA.Py = static_cast<float>(sexa.Py());
    fOutput_ChannelA.Pz = static_cast<float>(sexa.Pz());
    fOutput_ChannelA.E = static_cast<float>(sexa.E());
    fOutput_ChannelA.E_MinusNucleon = static_cast<float>(sexa.E_MinusNucleon());

    fOutput_ChannelA.V0A_X_AtPCA = static_cast<float>(sexa.V0A_PCA_XYZ()[0]);
    fOutput_ChannelA.V0A_Y_AtPCA = static_cast<float>(sexa.V0A_PCA_XYZ()[1]);
    fOutput_ChannelA.V0A_Z_AtPCA = static_cast<float>(sexa.V0A_PCA_XYZ()[2]);

    fOutput_ChannelA.V0B_X_AtPCA = static_cast<float>(sexa.V0B_PCA_XYZ()[0]);
    fOutput_ChannelA.V0B_Y_AtPCA = static_cast<float>(sexa.V0B_PCA_XYZ()[1]);
    fOutput_ChannelA.V0B_Z_AtPCA = static_cast<float>(sexa.V0B_PCA_XYZ()[2]);

    fOutput_ChannelA.V0A_Entry = sexa.V0A.idx;
    fOutput_ChannelA.V0A.X = static_cast<float>(sexa.V0A.X());
    fOutput_ChannelA.V0A.Y = static_cast<float>(sexa.V0A.Y());
    fOutput_ChannelA.V0A.Z = static_cast<float>(sexa.V0A.Z());
    fOutput_ChannelA.V0A.Px = static_cast<float>(sexa.V0A.Px());
    fOutput_ChannelA.V0A.Py = static_cast<float>(sexa.V0A.Py());
    fOutput_ChannelA.V0A.Pz = static_cast<float>(sexa.V0A.Pz());
    fOutput_ChannelA.V0A.E = static_cast<float>(sexa.V0A.E());

    fOutput_ChannelA.V0A_Neg_Entry = sexa.V0A.Neg.idx;
    fOutput_ChannelA.V0A_Pos_Entry = sexa.V0A.Pos.idx;

    fOutput_ChannelA.V0B_Entry = sexa.V0B.idx;
    fOutput_ChannelA.V0B.X = static_cast<float>(sexa.V0B.X());
    fOutput_ChannelA.V0B.Y = static_cast<float>(sexa.V0B.Y());
    fOutput_ChannelA.V0B.Z = static_cast<float>(sexa.V0B.Z());
    fOutput_ChannelA.V0B.Px = static_cast<float>(sexa.V0B.Px());
    fOutput_ChannelA.V0B.Py = static_cast<float>(sexa.V0B.Py());
    fOutput_ChannelA.V0B.Pz = static_cast<float>(sexa.V0B.Pz());
    fOutput_ChannelA.V0B.E = static_cast<float>(sexa.V0B.E());

    fOutput_ChannelA.V0B_Neg_Entry = sexa.V0B.Neg.idx;
    fOutput_ChannelA.V0B_Pos_Entry = sexa.V0B.Pos.idx;
}

void Finder::StoreMC(const MC::ChannelA& sexa) {

    fOutput_MC_ChannelA.X = static_cast<float>(sexa.X);
    fOutput_MC_ChannelA.Y = static_cast<float>(sexa.Y);
    fOutput_MC_ChannelA.Z = static_cast<float>(sexa.Z);
    fOutput_MC_ChannelA.BeforePx = static_cast<float>(sexa.BeforePx);
    fOutput_MC_ChannelA.BeforePy = static_cast<float>(sexa.BeforePy);
    fOutput_MC_ChannelA.BeforePz = static_cast<float>(sexa.BeforePz);
    fOutput_MC_ChannelA.BeforeE = static_cast<float>(sexa.BeforeE);

    fOutput_MC_ChannelA.NucleonPx = static_cast<float>(sexa.NucleonPx);
    fOutput_MC_ChannelA.NucleonPy = static_cast<float>(sexa.NucleonPy);
    fOutput_MC_ChannelA.NucleonPz = static_cast<float>(sexa.NucleonPz);
    fOutput_MC_ChannelA.NucleonE = static_cast<float>(sexa.NucleonE);

    fOutput_MC_ChannelA.AfterPx = static_cast<float>(sexa.AfterPx);
    fOutput_MC_ChannelA.AfterPy = static_cast<float>(sexa.AfterPy);
    fOutput_MC_ChannelA.AfterPz = static_cast<float>(sexa.AfterPz);
    fOutput_MC_ChannelA.AfterE = static_cast<float>(sexa.AfterE);

    fOutput_MC_ChannelA.IsSignal = sexa.IsSignal;
    fOutput_MC_ChannelA.ReactionID = sexa.ReactionID;
    fOutput_MC_ChannelA.IsHybrid = sexa.IsHybrid;

    fOutput_MC_ChannelA.V0A_Entry = sexa.V0A.Entry;
    fOutput_MC_ChannelA.V0A_PdgCode = sexa.V0A.PdgCode;
    fOutput_MC_ChannelA.V0A_Mother_Entry = sexa.V0A.Mother_Entry;
    fOutput_MC_ChannelA.V0A_Mother_PdgCode = sexa.V0A.Mother_PdgCode;
    fOutput_MC_ChannelA.V0A_Neg_Entry = sexa.V0A.neg.Entry;
    fOutput_MC_ChannelA.V0A_Pos_Entry = sexa.V0A.pos.Entry;
    fOutput_MC_ChannelA.V0A_Px = sexa.V0A.Px;
    fOutput_MC_ChannelA.V0A_Py = sexa.V0A.Py;
    fOutput_MC_ChannelA.V0A_Pz = sexa.V0A.Pz;
    fOutput_MC_ChannelA.V0A_E = sexa.V0A.Energy;
    fOutput_MC_ChannelA.V0A_IsTrue = sexa.V0A.IsTrue;
    fOutput_MC_ChannelA.V0A_IsSignal = sexa.V0A.IsSignal;
    fOutput_MC_ChannelA.V0A_IsSecondary = sexa.V0A.IsSecondary;
    fOutput_MC_ChannelA.V0A_ReactionID = sexa.V0A.ReactionID;
    fOutput_MC_ChannelA.V0A_IsHybrid = sexa.V0A.IsHybrid;

    fOutput_MC_ChannelA.V0B_Entry = sexa.V0B.Entry;
    fOutput_MC_ChannelA.V0B_PdgCode = sexa.V0B.PdgCode;
    fOutput_MC_ChannelA.V0B_Mother_Entry = sexa.V0B.Mother_Entry;
    fOutput_MC_ChannelA.V0B_Mother_PdgCode = sexa.V0B.Mother_PdgCode;
    fOutput_MC_ChannelA.V0B_Neg_Entry = sexa.V0B.neg.Entry;
    fOutput_MC_ChannelA.V0B_Pos_Entry = sexa.V0B.pos.Entry;
    fOutput_MC_ChannelA.V0B_Px = sexa.V0B.Px;
    fOutput_MC_ChannelA.V0B_Py = sexa.V0B.Py;
    fOutput_MC_ChannelA.V0B_Pz = sexa.V0B.Pz;
    fOutput_MC_ChannelA.V0B_E = sexa.V0B.Energy;
    fOutput_MC_ChannelA.V0B_IsTrue = sexa.V0B.IsTrue;
    fOutput_MC_ChannelA.V0B_IsSignal = sexa.V0B.IsSignal;
    fOutput_MC_ChannelA.V0B_IsSecondary = sexa.V0B.IsSecondary;
    fOutput_MC_ChannelA.V0B_ReactionID = sexa.V0B.ReactionID;
    fOutput_MC_ChannelA.V0B_IsHybrid = sexa.V0B.IsHybrid;
}

// ## Channel D ZONE ## //

void Finder::FindSexaquarks_ChannelD(bool anti_channel) {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    const PackedEvents::V0s* Packed_Lambdas{anti_channel ? &fPacked_Lambdas : &fPacked_AntiLambdas};
    const PackedEvents::Tracks* Packed_Kaons{anti_channel ? &fPacked_NegKaons : &fPacked_PosKaons};
    const PackedEvents::MC_V0s* MC_Lambdas{anti_channel ? &fPacked_MC_Lambdas : &fPacked_MC_AntiLambdas};
    const PackedEvents::MC_Tracks* MC_Kaons{anti_channel ? &fPacked_MC_NegKaons : &fPacked_MC_PosKaons};
    EParticle PID_Lambda{anti_channel ? EParticle::Lambda : EParticle::AntiLambda};
    int Charge_Kaon{anti_channel ? -1 : +1};

    // loop over all possible pairs of (anti)lambda + (pos/neg)kaon //
    auto n_lambdas = static_cast<int>(Packed_Lambdas->Entry->size());
    auto n_kaons = static_cast<int>(Packed_Kaons->Entry->size());
    for (int idx_v0{0}; idx_v0 < n_lambdas; ++idx_v0) {

        // unpack (anti)lambda //
        auto v0 = KF::UnpackV0(*Packed_Lambdas, idx_v0, PID_Lambda);

        for (int idx_kaon{0}; idx_kaon < n_kaons; ++idx_kaon) {

            // unpack kaon //
            auto kaon = KF::UnpackTrack(*Packed_Kaons, idx_kaon, Charge_Kaon);

            // sanity check //
            std::set<int> unique_track_entries{v0.Neg.idx, v0.Pos.idx, kaon.idx};
            if (unique_track_entries.size() < 3) continue;

            // fit //
            KF::ChannelD sexa{v0, kaon};
            sexa.DoFit(fInput_Event.MagneticField);

            // apply cuts //
            if (!PassesCuts(sexa)) continue;
#ifdef T2S_DEBUG
            Logger::Debug(__FUNCTION__, "idx(sexa,v0,neg,pos,kaon)={},{},{},{}", sexa.V0.idx, sexa.V0.Neg.idx, sexa.V0.Pos.idx, sexa.Kaon.idx);
            Logger::Debug(__FUNCTION__, ";x,y,z={},{},{}", sexa.X(), sexa.Y(), sexa.Z());
            Logger::Debug(__FUNCTION__, ";x,y,z(v0)={},{},{}", sexa.V0_PCA_XYZ()[0], sexa.V0_PCA_XYZ()[1], sexa.V0_PCA_XYZ()[2]);
            Logger::Debug(__FUNCTION__, ";x,y,z(kaon)={},{},{}", sexa.Kaon_PCA_XYZ()[0], sexa.Kaon_PCA_XYZ()[1], sexa.Kaon_PCA_XYZ()[2]);
            Logger::Debug(__FUNCTION__, ";mass={}", sexa.Mass());
            Logger::Debug(__FUNCTION__, ";dca_v0_kaon={}", sexa.DCA_btw_V0_Kaon());
            Logger::Debug(__FUNCTION__, ";radius={}", sexa.Radius2D());
            Logger::Debug(__FUNCTION__, ";dca_v0={}", sexa.DCA_V0_wrt_SV());
            Logger::Debug(__FUNCTION__, ";dca_kaon={}", sexa.DCA_Kaon_wrt_SV());
            Logger::Debug(__FUNCTION__, ";pt={}", sexa.Pt());
            Logger::Debug(__FUNCTION__, ";eta={}", sexa.Eta());
            Logger::Debug(__FUNCTION__, ";cpa_pv={}", sexa.CPA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv));
#endif

            // store //
            Store(sexa);
            if (IsMC()) {
                MC::V0 mc_v0{*MC_Lambdas, v0.idx};
                MC::Track mc_kaon{*MC_Kaons, kaon.idx};
                MC::ChannelD mc_sexa{fInput_Injected, fSettings.SexaquarkMass, mc_v0, mc_kaon};
                StoreMC(mc_sexa);
            }
            fOutputTree->Fill();
        }
    }

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
}

bool Finder::PassesCuts(const KF::ChannelD& sexa) const {

    fCutFlowHist->Fill(0.);
    if (sexa.Radius2D() < Cuts::ChannelD::Min_Radius2D || sexa.Radius2D() > Cuts::ChannelD::Max_Radius2D) return false;
    fCutFlowHist->Fill(1.);
    if (sexa.AbsRapidity_MinusNucleon() > Cuts::ChannelD::AbsMax_Rapidity) return false;  // PENDING: kinematics, affected by Fermi motion
    fCutFlowHist->Fill(2.);
    if (sexa.CPA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv) < Cuts::ChannelD::Min_CPAwrtPV ||
        sexa.CPA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv) > Cuts::ChannelD::Max_CPAwrtPV) {
        return false;  // PENDING: kinematics, affected by Fermi motion
    }
    fCutFlowHist->Fill(3.);
    if (sexa.DCA_V0_wrt_SV() > Cuts::ChannelD::Max_DCALaSV) return false;
    fCutFlowHist->Fill(4.);
    if (sexa.DCA_Kaon_wrt_SV() > Cuts::ChannelD::Max_DCAKaSV) return false;
    fCutFlowHist->Fill(5.);
    if (sexa.DCA_V0Neg_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelD::Max_DCALaNegSV) return false;
    fCutFlowHist->Fill(6.);
    if (sexa.DCA_V0Pos_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelD::Max_DCALaPosSV) return false;
    fCutFlowHist->Fill(7.);
    if (sexa.DCA_btw_V0_Kaon() > Cuts::ChannelD::Max_DCAKaLa) return false;
    fCutFlowHist->Fill(8.);
    if (sexa.DCA_btw_V0Neg_Kaon(fInput_Event.MagneticField) > Cuts::ChannelD::Max_DCALaNegKa) return false;
    fCutFlowHist->Fill(9.);
    if (sexa.DCA_btw_V0Pos_Kaon(fInput_Event.MagneticField) > Cuts::ChannelD::Max_DCALaPosKa) return false;
    fCutFlowHist->Fill(10.);

    return true;
}

void Finder::Store(const KF::ChannelD& sexa) {

    fOutput_ChannelD.RunNumber = fInput_Event.RunNumber;
    fOutput_ChannelD.DirNumber = fInput_Event.DirNumber;
    if (!IsMC()) fOutput_ChannelD.DirNumberB = fInput_Event.DirNumberB;
    fOutput_ChannelD.EventNumber = fInput_Event.EventNumber;
    fOutput_ChannelD.MagneticField = fInput_Event.MagneticField;
    fOutput_ChannelD.PV_Xv = fInput_Event.PV_Xv;
    fOutput_ChannelD.PV_Yv = fInput_Event.PV_Yv;
    fOutput_ChannelD.PV_Zv = fInput_Event.PV_Zv;
    if (IsMC()) {
        fOutput_ChannelD.MC_PV_Xv = fInput_Event.MC_PV_Xv;
        fOutput_ChannelD.MC_PV_Yv = fInput_Event.MC_PV_Yv;
        fOutput_ChannelD.MC_PV_Zv = fInput_Event.MC_PV_Zv;
    }

    fOutput_ChannelD.X = static_cast<float>(sexa.X());
    fOutput_ChannelD.Y = static_cast<float>(sexa.Y());
    fOutput_ChannelD.Z = static_cast<float>(sexa.Z());
    fOutput_ChannelD.Px = static_cast<float>(sexa.Px());
    fOutput_ChannelD.Py = static_cast<float>(sexa.Py());
    fOutput_ChannelD.Pz = static_cast<float>(sexa.Pz());
    fOutput_ChannelD.E = static_cast<float>(sexa.E());
    fOutput_ChannelD.E_MinusNucleon = static_cast<float>(sexa.E_MinusNucleon());

    fOutput_ChannelD.V0_X_AtPCA = static_cast<float>(sexa.V0_PCA_XYZ()[0]);
    fOutput_ChannelD.V0_Y_AtPCA = static_cast<float>(sexa.V0_PCA_XYZ()[1]);
    fOutput_ChannelD.V0_Z_AtPCA = static_cast<float>(sexa.V0_PCA_XYZ()[2]);

    fOutput_ChannelD.Kaon_X_AtPCA = static_cast<float>(sexa.Kaon_PCA_XYZ()[0]);
    fOutput_ChannelD.Kaon_Y_AtPCA = static_cast<float>(sexa.Kaon_PCA_XYZ()[1]);
    fOutput_ChannelD.Kaon_Z_AtPCA = static_cast<float>(sexa.Kaon_PCA_XYZ()[2]);
    fOutput_ChannelD.Kaon_Px_AtPCA = static_cast<float>(sexa.Kaon_PCA_PxPyPz()[0]);
    fOutput_ChannelD.Kaon_Py_AtPCA = static_cast<float>(sexa.Kaon_PCA_PxPyPz()[1]);
    fOutput_ChannelD.Kaon_Pz_AtPCA = static_cast<float>(sexa.Kaon_PCA_PxPyPz()[2]);

    fOutput_ChannelD.V0_Entry = sexa.V0.idx;
    fOutput_ChannelD.V0.X = static_cast<float>(sexa.V0.X());
    fOutput_ChannelD.V0.Y = static_cast<float>(sexa.V0.Y());
    fOutput_ChannelD.V0.Z = static_cast<float>(sexa.V0.Z());
    fOutput_ChannelD.V0.Px = static_cast<float>(sexa.V0.Px());
    fOutput_ChannelD.V0.Py = static_cast<float>(sexa.V0.Py());
    fOutput_ChannelD.V0.Pz = static_cast<float>(sexa.V0.Pz());
    fOutput_ChannelD.V0.E = static_cast<float>(sexa.V0.E());

    fOutput_ChannelD.V0_Neg_Entry = sexa.V0.Neg.idx;
    fOutput_ChannelD.V0_Pos_Entry = sexa.V0.Pos.idx;

    fOutput_ChannelD.Kaon_Entry = sexa.Kaon.idx;
    fOutput_ChannelD.Kaon.X = static_cast<float>(sexa.Kaon.X());
    fOutput_ChannelD.Kaon.Y = static_cast<float>(sexa.Kaon.Y());
    fOutput_ChannelD.Kaon.Z = static_cast<float>(sexa.Kaon.Z());
    fOutput_ChannelD.Kaon.Px = static_cast<float>(sexa.Kaon.Px());
    fOutput_ChannelD.Kaon.Py = static_cast<float>(sexa.Kaon.Py());
    fOutput_ChannelD.Kaon.Pz = static_cast<float>(sexa.Kaon.Pz());
    fOutput_ChannelD.Kaon.E = static_cast<float>(sexa.Kaon.E());
}

void Finder::StoreMC(const MC::ChannelD& sexa) {

    fOutput_MC_ChannelD.X = static_cast<float>(sexa.X);
    fOutput_MC_ChannelD.Y = static_cast<float>(sexa.Y);
    fOutput_MC_ChannelD.Z = static_cast<float>(sexa.Z);
    fOutput_MC_ChannelD.BeforePx = static_cast<float>(sexa.BeforePx);
    fOutput_MC_ChannelD.BeforePy = static_cast<float>(sexa.BeforePy);
    fOutput_MC_ChannelD.BeforePz = static_cast<float>(sexa.BeforePz);
    fOutput_MC_ChannelD.BeforeE = static_cast<float>(sexa.BeforeE);

    fOutput_MC_ChannelD.NucleonPx = static_cast<float>(sexa.NucleonPx);
    fOutput_MC_ChannelD.NucleonPy = static_cast<float>(sexa.NucleonPy);
    fOutput_MC_ChannelD.NucleonPz = static_cast<float>(sexa.NucleonPz);
    fOutput_MC_ChannelD.NucleonE = static_cast<float>(sexa.NucleonE);

    fOutput_MC_ChannelD.AfterPx = static_cast<float>(sexa.AfterPx);
    fOutput_MC_ChannelD.AfterPy = static_cast<float>(sexa.AfterPy);
    fOutput_MC_ChannelD.AfterPz = static_cast<float>(sexa.AfterPz);
    fOutput_MC_ChannelD.AfterE = static_cast<float>(sexa.AfterE);

    fOutput_MC_ChannelD.IsSignal = sexa.IsSignal;
    fOutput_MC_ChannelD.ReactionID = sexa.ReactionID;
    fOutput_MC_ChannelD.IsHybrid = sexa.IsHybrid;

    fOutput_MC_ChannelD.V0_Entry = sexa.V0.Entry;
    fOutput_MC_ChannelD.V0_PdgCode = sexa.V0.PdgCode;
    fOutput_MC_ChannelD.V0_Mother_Entry = sexa.V0.Mother_Entry;
    fOutput_MC_ChannelD.V0_Mother_PdgCode = sexa.V0.Mother_PdgCode;
    fOutput_MC_ChannelD.V0_Neg_Entry = sexa.V0.neg.Entry;
    fOutput_MC_ChannelD.V0_Pos_Entry = sexa.V0.pos.Entry;
    fOutput_MC_ChannelD.V0_Px = sexa.V0.Px;
    fOutput_MC_ChannelD.V0_Py = sexa.V0.Py;
    fOutput_MC_ChannelD.V0_Pz = sexa.V0.Pz;
    fOutput_MC_ChannelD.V0_E = sexa.V0.Energy;
    fOutput_MC_ChannelD.V0_IsTrue = sexa.V0.IsTrue;
    fOutput_MC_ChannelD.V0_IsSignal = sexa.V0.IsSignal;
    fOutput_MC_ChannelD.V0_IsSecondary = sexa.V0.IsSecondary;
    fOutput_MC_ChannelD.V0_ReactionID = sexa.V0.ReactionID;
    fOutput_MC_ChannelD.V0_IsHybrid = sexa.V0.IsHybrid;

    fOutput_MC_ChannelD.Kaon_Entry = sexa.Kaon.Entry;
    fOutput_MC_ChannelD.Kaon_PdgCode = sexa.Kaon.PdgCode;
    fOutput_MC_ChannelD.Kaon_Mother_Entry = sexa.Kaon.Mother_Entry;
    fOutput_MC_ChannelD.Kaon_Mother_PdgCode = sexa.Kaon.Mother_PdgCode;
    fOutput_MC_ChannelD.Kaon_GrandMother_Entry = sexa.Kaon.GrandMother_Entry;
    fOutput_MC_ChannelD.Kaon_GrandMother_PdgCode = sexa.Kaon.GrandMother_PdgCode;
    fOutput_MC_ChannelD.Kaon_Px = sexa.Kaon.Px;
    fOutput_MC_ChannelD.Kaon_Py = sexa.Kaon.Py;
    fOutput_MC_ChannelD.Kaon_Pz = sexa.Kaon.Pz;
    fOutput_MC_ChannelD.Kaon_E = sexa.Kaon.Energy;
    fOutput_MC_ChannelD.Kaon_IsTrue = sexa.Kaon.IsTrue;
    fOutput_MC_ChannelD.Kaon_IsSignal = sexa.Kaon.IsSignal;
    fOutput_MC_ChannelD.Kaon_IsSecondary = sexa.Kaon.IsSecondary;
    fOutput_MC_ChannelD.Kaon_ReactionID = sexa.Kaon.ReactionID;
}

// ## Channel E ZONE ## //

void Finder::FindSexaquarks_ChannelE(bool anti_channel) {
    // PENDING
}

bool Finder::PassesCuts(const KF::ChannelE& sexa) const {
    // PENDING
    return true;
}

void Finder::Store(const KF::ChannelE& sexa) {
    // PENDING
}

// ## Channel H ZONE ## //

void Finder::FindSexaquarks_ChannelH(bool anti_channel) {
    // PENDING
}

bool Finder::PassesCuts(const KF::ChannelH& sexa) const {
    // PENDING
    return true;
}

void Finder::Store(const KF::ChannelH& sexa) {
    // PENDING
}

// ## END OF CYCLES ## //

void Finder::EndOfAnalysis() {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    if (IsMC()) {
        fOutputTree_Injected->Write();
        Logger::Info(__FUNCTION__, "TTree \"{}\" has been written onto TFile {}", fOutputTree_Injected->GetName(), fSettings.PathOutputFile);
    }

    fOutputTree->Write();
    Logger::Info(__FUNCTION__, "TTree \"{}\" has been written onto TFile {}", fOutputTree->GetName(), fSettings.PathOutputFile);

    fCutFlowHist->Write();

    fTree_PackedEvents->ResetBranchAddresses();
    fOutputTree->ResetBranchAddresses();

    Logger::Info(__FUNCTION__, "All done.");
}

}  // namespace Tree2Secondaries
