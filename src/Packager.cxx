#include <filesystem>
#include <memory>

#include "ALICE/ESD.hxx"
#include "ALICE/Utilities.hxx"
#include "App/Logger.hxx"
#include "App/Utilities.hxx"
#include "KF/Utilities.hxx"
#include "Math/Constants.hxx"
#include "Packager/Cuts.hxx"
#include "Packager/Packager.hxx"

#ifdef T2S_USE_ALICE
#include "ALICE/Vertexer.hxx"
#endif

namespace Tree2Secondaries {

bool Packager::Initialize() {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    fEventsTree = std::make_unique<TChain>("Events");
    for (const auto& path : fSettings.PathInputFiles) {
        if (fEventsTree->Add(path.c_str()) == 0) {
            Logger::Error(__FUNCTION__, "Couldn't add TFile {}", path);
        }
    }
    if (!fEventsTree->GetEntries()) {
        Logger::Error(__FUNCTION__, "Couldn't read any entry.");
        return false;
    }
    Logger::Info(__FUNCTION__, "TChain \"{}\" loaded successfully with {} trees and {} total entries.", fEventsTree->GetName(),
                 fEventsTree->GetNtrees(), fEventsTree->GetEntries());

    ConnectInputBranches();

    if (!PrepareOutputFile()) return false;

    CreateCutFlowHistograms();

    if (!PrepareOutputTree()) return false;
    CreateOutputBranches();

    Logger::Info(__FUNCTION__, "Packager initialized successfully.");

    return true;
}

// ## INPUT ZONE ## //

void Packager::ConnectInputBranches() {
    fEventsTree->SetBranchStatus("*", false);

    ConnectBranches_Events();
    if (IsMC()) {
        ConnectBranches_MC();
        if (IsMC()) ConnectBranches_Injected();
    }
    ConnectBranches_Tracks();
}

void Packager::ConnectBranches_Events() {

    Utils::ConnectBranch(fEventsTree.get(), "RunNumber", &fInput_Event.RunNumber);
    Utils::ConnectBranch(fEventsTree.get(), "DirNumber", &fInput_Event.DirNumber);
    if (!IsMC()) Utils::ConnectBranch(fEventsTree.get(), "DirNumberB", &fInput_Event.DirNumberB);
    Utils::ConnectBranch(fEventsTree.get(), "EventNumber", &fInput_Event.EventNumber);
    Utils::ConnectBranch(fEventsTree.get(), "Centrality", &fInput_Event.Centrality);
    Utils::ConnectBranch(fEventsTree.get(), "MagneticField", &fInput_Event.MagneticField);
    Utils::ConnectBranch(fEventsTree.get(), "PV_Xv", &fInput_Event.PV_Xv);
    Utils::ConnectBranch(fEventsTree.get(), "PV_Yv", &fInput_Event.PV_Yv);
    Utils::ConnectBranch(fEventsTree.get(), "PV_Zv", &fInput_Event.PV_Zv);

    if (IsMC()) {
        Utils::ConnectBranch(fEventsTree.get(), "MC_PV_Xv", &fInput_Event.MC_PV_Xv);
        Utils::ConnectBranch(fEventsTree.get(), "MC_PV_Yv", &fInput_Event.MC_PV_Yv);
        Utils::ConnectBranch(fEventsTree.get(), "MC_PV_Zv", &fInput_Event.MC_PV_Zv);
    }
}

void Packager::ConnectBranches_Injected() {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    Utils::ConnectBranch(fEventsTree.get(), "ReactionID", &fInput_Injected.ReactionID);
    Utils::ConnectBranch(fEventsTree.get(), "Sexaquark_Px", &fInput_Injected.Px);
    Utils::ConnectBranch(fEventsTree.get(), "Sexaquark_Py", &fInput_Injected.Py);
    Utils::ConnectBranch(fEventsTree.get(), "Sexaquark_Pz", &fInput_Injected.Pz);
    Utils::ConnectBranch(fEventsTree.get(), "Nucleon_Px", &fInput_Injected.Nucleon_Px);
    Utils::ConnectBranch(fEventsTree.get(), "Nucleon_Py", &fInput_Injected.Nucleon_Py);
    Utils::ConnectBranch(fEventsTree.get(), "Nucleon_Pz", &fInput_Injected.Nucleon_Pz);

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
}

void Packager::ConnectBranches_MC() {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    Utils::ConnectBranch(fEventsTree.get(), "MC_X", &fInput_MC.X);
    Utils::ConnectBranch(fEventsTree.get(), "MC_Y", &fInput_MC.Y);
    Utils::ConnectBranch(fEventsTree.get(), "MC_Z", &fInput_MC.Z);
    Utils::ConnectBranch(fEventsTree.get(), "MC_Px", &fInput_MC.Px);
    Utils::ConnectBranch(fEventsTree.get(), "MC_Py", &fInput_MC.Py);
    Utils::ConnectBranch(fEventsTree.get(), "MC_Pz", &fInput_MC.Pz);
    Utils::ConnectBranch(fEventsTree.get(), "MC_E", &fInput_MC.E);

    Utils::ConnectBranch(fEventsTree.get(), "MC_PdgCode", &fInput_MC.PdgCode);
    Utils::ConnectBranch(fEventsTree.get(), "MC_Mother_McEntry", &fInput_MC.MotherEntry);
    Utils::ConnectBranch(fEventsTree.get(), "MC_Status", &fInput_MC.Status);
    Utils::ConnectBranch(fEventsTree.get(), "MC_Generator", &fInput_MC.Generator);
    Utils::ConnectBranch(fEventsTree.get(), "MC_IsPrimary", &fInput_MC.IsPrimary);
    Utils::ConnectBranch(fEventsTree.get(), "MC_IsSecFromMat", &fInput_MC.IsSecFromMat);
    Utils::ConnectBranch(fEventsTree.get(), "MC_IsSecFromWeak", &fInput_MC.IsSecFromWeak);

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
}

void Packager::ConnectBranches_Tracks() {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    Utils::ConnectBranch(fEventsTree.get(), "Track_X", &fInput_Tracks.X);
    Utils::ConnectBranch(fEventsTree.get(), "Track_Y", &fInput_Tracks.Y);
    Utils::ConnectBranch(fEventsTree.get(), "Track_Z", &fInput_Tracks.Z);
    Utils::ConnectBranch(fEventsTree.get(), "Track_Px", &fInput_Tracks.Px);
    Utils::ConnectBranch(fEventsTree.get(), "Track_Py", &fInput_Tracks.Py);
    Utils::ConnectBranch(fEventsTree.get(), "Track_Pz", &fInput_Tracks.Pz);
    Utils::ConnectBranch(fEventsTree.get(), "Track_Charge", &fInput_Tracks.Charge);
    Utils::ConnectBranch(fEventsTree.get(), "Track_NSigmaPion", &fInput_Tracks.NSigmaPion);
    Utils::ConnectBranch(fEventsTree.get(), "Track_NSigmaKaon", &fInput_Tracks.NSigmaKaon);
    Utils::ConnectBranch(fEventsTree.get(), "Track_NSigmaProton", &fInput_Tracks.NSigmaProton);

    Utils::ConnectBranch(fEventsTree.get(), "Track_Alpha", &fInput_Tracks.Alpha);
    Utils::ConnectBranch(fEventsTree.get(), "Track_Snp", &fInput_Tracks.Snp);
    Utils::ConnectBranch(fEventsTree.get(), "Track_Tgl", &fInput_Tracks.Tgl);
    Utils::ConnectBranch(fEventsTree.get(), "Track_Signed1Pt", &fInput_Tracks.Signed1Pt);
    Utils::ConnectBranch(fEventsTree.get(), "Track_SigmaY2", &fInput_Tracks.SigmaY2);
    Utils::ConnectBranch(fEventsTree.get(), "Track_SigmaZY", &fInput_Tracks.SigmaZY);
    Utils::ConnectBranch(fEventsTree.get(), "Track_SigmaZ2", &fInput_Tracks.SigmaZ2);
    Utils::ConnectBranch(fEventsTree.get(), "Track_SigmaSnpY", &fInput_Tracks.SigmaSnpY);
    Utils::ConnectBranch(fEventsTree.get(), "Track_SigmaSnpZ", &fInput_Tracks.SigmaSnpZ);
    Utils::ConnectBranch(fEventsTree.get(), "Track_SigmaSnp2", &fInput_Tracks.SigmaSnp2);
    Utils::ConnectBranch(fEventsTree.get(), "Track_SigmaTglY", &fInput_Tracks.SigmaTglY);
    Utils::ConnectBranch(fEventsTree.get(), "Track_SigmaTglZ", &fInput_Tracks.SigmaTglZ);
    Utils::ConnectBranch(fEventsTree.get(), "Track_SigmaTglSnp", &fInput_Tracks.SigmaTglSnp);
    Utils::ConnectBranch(fEventsTree.get(), "Track_SigmaTgl2", &fInput_Tracks.SigmaTgl2);
    Utils::ConnectBranch(fEventsTree.get(), "Track_Sigma1PtY", &fInput_Tracks.Sigma1PtY);
    Utils::ConnectBranch(fEventsTree.get(), "Track_Sigma1PtZ", &fInput_Tracks.Sigma1PtZ);
    Utils::ConnectBranch(fEventsTree.get(), "Track_Sigma1PtSnp", &fInput_Tracks.Sigma1PtSnp);
    Utils::ConnectBranch(fEventsTree.get(), "Track_Sigma1PtTgl", &fInput_Tracks.Sigma1PtTgl);
    Utils::ConnectBranch(fEventsTree.get(), "Track_Sigma1Pt2", &fInput_Tracks.Sigma1Pt2);

    if (IsMC()) Utils::ConnectBranch(fEventsTree.get(), "Track_McEntry", &fInput_Tracks.McEntry);

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
}

// ## OUTPUT ZONE ## //

bool Packager::PrepareOutputFile() {
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

bool Packager::PrepareOutputTree() {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    fOutputTree = std::make_unique<TTree>("PackedEvents", "Packed Events");
    if (!fOutputTree) {
        Logger::Error(__FUNCTION__, "Couldn't create TTree \"PackedEvents\"");
        return false;
    }

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
    return true;
}

void Packager::CreateOutputBranches() {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    CreateOutputBranches_Events();
    if (IsMC()) CreateOutputBranches_Injected();

    switch (GetReactionChannel()) {
        // standard channels //
        case EReactionChannel::A:
            CreateOutputBranches_V0s(EParticle::AntiLambda, fOutput_AntiLambdas);
            CreateOutputBranches_V0s(EParticle::KaonZeroShort, fOutput_KaonsZeroShort);
            if (IsMC()) {
                CreateOutputBranches_MC_V0s(EParticle::AntiLambda, fOutput_MC_AntiLambdas);
                CreateOutputBranches_MC_V0s(EParticle::KaonZeroShort, fOutput_MC_KaonsZeroShort);
            }
            break;
        case EReactionChannel::D:
            CreateOutputBranches_V0s(EParticle::AntiLambda, fOutput_AntiLambdas);
            CreateOutputBranches_Tracks(EParticle::PosKaon, fOutput_PosKaons);
            if (IsMC()) {
                CreateOutputBranches_MC_V0s(EParticle::AntiLambda, fOutput_MC_AntiLambdas);
                CreateOutputBranches_MC_Tracks(EParticle::PosKaon, fOutput_MC_PosKaons);
            }
            break;
        case EReactionChannel::E:
            CreateOutputBranches_V0s(EParticle::AntiLambda, fOutput_AntiLambdas);
            CreateOutputBranches_Tracks(EParticle::PosKaon, fOutput_PosKaons);
            CreateOutputBranches_Tracks(EParticle::PiMinus, fOutput_PiMinus);
            CreateOutputBranches_Tracks(EParticle::PiPlus, fOutput_PiPlus);
            if (IsMC()) {
                CreateOutputBranches_MC_V0s(EParticle::AntiLambda, fOutput_MC_AntiLambdas);
                CreateOutputBranches_MC_Tracks(EParticle::PosKaon, fOutput_MC_PosKaons);
                CreateOutputBranches_MC_Tracks(EParticle::PiMinus, fOutput_MC_PiMinus);
                CreateOutputBranches_MC_Tracks(EParticle::PiPlus, fOutput_MC_PiPlus);
            }
            break;
        case EReactionChannel::H:
            CreateOutputBranches_Tracks(EParticle::PosKaon, fOutput_PosKaons);
            if (IsMC()) CreateOutputBranches_MC_Tracks(EParticle::PosKaon, fOutput_MC_PosKaons);
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
                CreateOutputBranches_MC_V0s(EParticle::Lambda, fOutput_MC_Lambdas);
                CreateOutputBranches_MC_Tracks(EParticle::NegKaon, fOutput_MC_NegKaons);
            }
            break;
        case EReactionChannel::AntiE:
            CreateOutputBranches_V0s(EParticle::Lambda, fOutput_Lambdas);
            CreateOutputBranches_Tracks(EParticle::NegKaon, fOutput_NegKaons);
            CreateOutputBranches_Tracks(EParticle::PiMinus, fOutput_PiMinus);
            CreateOutputBranches_Tracks(EParticle::PiPlus, fOutput_PiPlus);
            if (IsMC()) {
                CreateOutputBranches_MC_V0s(EParticle::Lambda, fOutput_MC_Lambdas);
                CreateOutputBranches_MC_Tracks(EParticle::NegKaon, fOutput_MC_NegKaons);
                CreateOutputBranches_MC_Tracks(EParticle::PiMinus, fOutput_MC_PiMinus);
                CreateOutputBranches_MC_Tracks(EParticle::PiPlus, fOutput_MC_PiPlus);
            }
            break;
        case EReactionChannel::AntiH:
            CreateOutputBranches_Tracks(EParticle::NegKaon, fOutput_NegKaons);
            if (IsMC()) CreateOutputBranches_MC_Tracks(EParticle::NegKaon, fOutput_MC_NegKaons);
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
                CreateOutputBranches_MC_V0s(EParticle::AntiLambda, fOutput_MC_AntiLambdas);
                CreateOutputBranches_MC_V0s(EParticle::Lambda, fOutput_MC_Lambdas);
                CreateOutputBranches_MC_V0s(EParticle::KaonZeroShort, fOutput_MC_KaonsZeroShort);
                CreateOutputBranches_MC_Tracks(EParticle::NegKaon, fOutput_MC_NegKaons);
                CreateOutputBranches_MC_Tracks(EParticle::PosKaon, fOutput_MC_PosKaons);
                CreateOutputBranches_MC_Tracks(EParticle::PiMinus, fOutput_MC_PiMinus);
                CreateOutputBranches_MC_Tracks(EParticle::PiPlus, fOutput_MC_PiPlus);
            }
            break;
    }  // end of switch statement

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
}

void Packager::CreateOutputBranches_Events() {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    Utils::CreateBranch(fOutputTree.get(), "RunNumber", &fOutput_Event.RunNumber);
    Utils::CreateBranch(fOutputTree.get(), "DirNumber", &fOutput_Event.DirNumber);
    if (!IsMC()) Utils::CreateBranch(fOutputTree.get(), "DirNumberB", &fOutput_Event.DirNumberB);
    Utils::CreateBranch(fOutputTree.get(), "EventNumber", &fOutput_Event.EventNumber);
    Utils::CreateBranch(fOutputTree.get(), "Centrality", &fOutput_Event.Centrality);
    Utils::CreateBranch(fOutputTree.get(), "MagneticField", &fOutput_Event.MagneticField);
    Utils::CreateBranch(fOutputTree.get(), "PV_Xv", &fOutput_Event.PV_Xv);
    Utils::CreateBranch(fOutputTree.get(), "PV_Yv", &fOutput_Event.PV_Yv);
    Utils::CreateBranch(fOutputTree.get(), "PV_Zv", &fOutput_Event.PV_Zv);

    if (IsMC()) {
        Utils::CreateBranch(fOutputTree.get(), "MC_PV_Xv", &fOutput_Event.MC_PV_Xv);
        Utils::CreateBranch(fOutputTree.get(), "MC_PV_Yv", &fOutput_Event.MC_PV_Yv);
        Utils::CreateBranch(fOutputTree.get(), "MC_PV_Zv", &fOutput_Event.MC_PV_Zv);
    }

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
}

void Packager::CreateOutputBranches_Injected() {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

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

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
}

void Packager::CreateOutputBranches_V0s(EParticle pid, PackedEvents::V0s& sov) {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Entry", Particle::Acronym[pid]), &sov.Entry);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_X", Particle::Acronym[pid]), &sov.X);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Y", Particle::Acronym[pid]), &sov.Y);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Z", Particle::Acronym[pid]), &sov.Z);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Px", Particle::Acronym[pid]), &sov.Px);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Py", Particle::Acronym[pid]), &sov.Py);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Pz", Particle::Acronym[pid]), &sov.Pz);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_E", Particle::Acronym[pid]), &sov.E);

    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaX2", Particle::Acronym[pid]), &sov.Sigma.X2);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaXY", Particle::Acronym[pid]), &sov.Sigma.XY);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaY2", Particle::Acronym[pid]), &sov.Sigma.Y2);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaXZ", Particle::Acronym[pid]), &sov.Sigma.XZ);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaYZ", Particle::Acronym[pid]), &sov.Sigma.YZ);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaZ2", Particle::Acronym[pid]), &sov.Sigma.Z2);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaXPx", Particle::Acronym[pid]), &sov.Sigma.XPx);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaYPx", Particle::Acronym[pid]), &sov.Sigma.YPx);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaZPx", Particle::Acronym[pid]), &sov.Sigma.ZPx);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaPx2", Particle::Acronym[pid]), &sov.Sigma.Px2);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaXPy", Particle::Acronym[pid]), &sov.Sigma.XPy);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaYPy", Particle::Acronym[pid]), &sov.Sigma.YPy);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaZPy", Particle::Acronym[pid]), &sov.Sigma.ZPy);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaPxPy", Particle::Acronym[pid]), &sov.Sigma.PxPy);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaPy2", Particle::Acronym[pid]), &sov.Sigma.Py2);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaXPz", Particle::Acronym[pid]), &sov.Sigma.XPz);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaYPz", Particle::Acronym[pid]), &sov.Sigma.YPz);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaZPz", Particle::Acronym[pid]), &sov.Sigma.ZPz);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaPxPz", Particle::Acronym[pid]), &sov.Sigma.PxPz);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaPyPz", Particle::Acronym[pid]), &sov.Sigma.PyPz);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaPz2", Particle::Acronym[pid]), &sov.Sigma.Pz2);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaXE", Particle::Acronym[pid]), &sov.Sigma.XE);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaYE", Particle::Acronym[pid]), &sov.Sigma.YE);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaZE", Particle::Acronym[pid]), &sov.Sigma.ZE);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaPxE", Particle::Acronym[pid]), &sov.Sigma.PxE);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaPyE", Particle::Acronym[pid]), &sov.Sigma.PyE);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaPzE", Particle::Acronym[pid]), &sov.Sigma.PzE);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaE2", Particle::Acronym[pid]), &sov.Sigma.E2);

    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Neg_Entry", Particle::Acronym[pid]), &sov.Neg.Entry);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Neg_X", Particle::Acronym[pid]), &sov.Neg.X);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Neg_Y", Particle::Acronym[pid]), &sov.Neg.Y);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Neg_Z", Particle::Acronym[pid]), &sov.Neg.Z);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Neg_Px", Particle::Acronym[pid]), &sov.Neg.Px);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Neg_Py", Particle::Acronym[pid]), &sov.Neg.Py);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Neg_Pz", Particle::Acronym[pid]), &sov.Neg.Pz);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Neg_E", Particle::Acronym[pid]), &sov.Neg.E);

    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Neg_X_AtPCA", Particle::Acronym[pid]), &sov.Neg_X_AtPCA);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Neg_Y_AtPCA", Particle::Acronym[pid]), &sov.Neg_Y_AtPCA);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Neg_Z_AtPCA", Particle::Acronym[pid]), &sov.Neg_Z_AtPCA);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Neg_Px_AtPCA", Particle::Acronym[pid]), &sov.Neg_Px_AtPCA);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Neg_Py_AtPCA", Particle::Acronym[pid]), &sov.Neg_Py_AtPCA);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Neg_Pz_AtPCA", Particle::Acronym[pid]), &sov.Neg_Pz_AtPCA);

    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Pos_Entry", Particle::Acronym[pid]), &sov.Pos.Entry);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Pos_X", Particle::Acronym[pid]), &sov.Pos.X);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Pos_Y", Particle::Acronym[pid]), &sov.Pos.Y);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Pos_Z", Particle::Acronym[pid]), &sov.Pos.Z);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Pos_Px", Particle::Acronym[pid]), &sov.Pos.Px);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Pos_Py", Particle::Acronym[pid]), &sov.Pos.Py);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Pos_Pz", Particle::Acronym[pid]), &sov.Pos.Pz);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Pos_E", Particle::Acronym[pid]), &sov.Pos.E);

    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Pos_X_AtPCA", Particle::Acronym[pid]), &sov.Pos_X_AtPCA);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Pos_Y_AtPCA", Particle::Acronym[pid]), &sov.Pos_Y_AtPCA);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Pos_Z_AtPCA", Particle::Acronym[pid]), &sov.Pos_Z_AtPCA);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Pos_Px_AtPCA", Particle::Acronym[pid]), &sov.Pos_Px_AtPCA);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Pos_Py_AtPCA", Particle::Acronym[pid]), &sov.Pos_Py_AtPCA);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Pos_Pz_AtPCA", Particle::Acronym[pid]), &sov.Pos_Pz_AtPCA);

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
}

void Packager::CreateOutputBranches_Tracks(EParticle pid, PackedEvents::Tracks& sov) {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Entry", Particle::Acronym[pid]), &sov.Entry);

    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_X", Particle::Acronym[pid]), &sov.X);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Y", Particle::Acronym[pid]), &sov.Y);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Z", Particle::Acronym[pid]), &sov.Z);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Px", Particle::Acronym[pid]), &sov.Px);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Py", Particle::Acronym[pid]), &sov.Py);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_Pz", Particle::Acronym[pid]), &sov.Pz);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_E", Particle::Acronym[pid]), &sov.E);

    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaX2", Particle::Acronym[pid]), &sov.Sigma.X2);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaXY", Particle::Acronym[pid]), &sov.Sigma.XY);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaY2", Particle::Acronym[pid]), &sov.Sigma.Y2);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaXZ", Particle::Acronym[pid]), &sov.Sigma.XZ);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaYZ", Particle::Acronym[pid]), &sov.Sigma.YZ);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaZ2", Particle::Acronym[pid]), &sov.Sigma.Z2);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaXPx", Particle::Acronym[pid]), &sov.Sigma.XPx);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaYPx", Particle::Acronym[pid]), &sov.Sigma.YPx);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaZPx", Particle::Acronym[pid]), &sov.Sigma.ZPx);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaPx2", Particle::Acronym[pid]), &sov.Sigma.Px2);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaXPy", Particle::Acronym[pid]), &sov.Sigma.XPy);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaYPy", Particle::Acronym[pid]), &sov.Sigma.YPy);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaZPy", Particle::Acronym[pid]), &sov.Sigma.ZPy);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaPxPy", Particle::Acronym[pid]), &sov.Sigma.PxPy);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaPy2", Particle::Acronym[pid]), &sov.Sigma.Py2);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaXPz", Particle::Acronym[pid]), &sov.Sigma.XPz);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaYPz", Particle::Acronym[pid]), &sov.Sigma.YPz);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaZPz", Particle::Acronym[pid]), &sov.Sigma.ZPz);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaPxPz", Particle::Acronym[pid]), &sov.Sigma.PxPz);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaPyPz", Particle::Acronym[pid]), &sov.Sigma.PyPz);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaPz2", Particle::Acronym[pid]), &sov.Sigma.Pz2);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaXE", Particle::Acronym[pid]), &sov.Sigma.XE);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaYE", Particle::Acronym[pid]), &sov.Sigma.YE);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaZE", Particle::Acronym[pid]), &sov.Sigma.ZE);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaPxE", Particle::Acronym[pid]), &sov.Sigma.PxE);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaPyE", Particle::Acronym[pid]), &sov.Sigma.PyE);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaPzE", Particle::Acronym[pid]), &sov.Sigma.PzE);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_SigmaE2", Particle::Acronym[pid]), &sov.Sigma.E2);

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
}

void Packager::CreateOutputBranches_MC_V0s(EParticle pid, PackedEvents::MC_V0s& sov) {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Entry", Particle::Acronym[pid]), &sov.Entry);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_X", Particle::Acronym[pid]), &sov.X);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Y", Particle::Acronym[pid]), &sov.Y);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Z", Particle::Acronym[pid]), &sov.Z);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Px", Particle::Acronym[pid]), &sov.Px);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Py", Particle::Acronym[pid]), &sov.Py);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Pz", Particle::Acronym[pid]), &sov.Pz);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_E", Particle::Acronym[pid]), &sov.E);

    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_DecayX", Particle::Acronym[pid]), &sov.DecayX);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_DecayY", Particle::Acronym[pid]), &sov.DecayY);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_DecayZ", Particle::Acronym[pid]), &sov.DecayZ);

    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_PdgCode", Particle::Acronym[pid]), &sov.PdgCode);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Mother_Entry", Particle::Acronym[pid]), &sov.Mother_Entry);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Mother_PdgCode", Particle::Acronym[pid]), &sov.Mother_PdgCode);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_IsTrue", Particle::Acronym[pid]), &sov.IsTrue);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_IsSignal", Particle::Acronym[pid]), &sov.IsSignal);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_IsSecondary", Particle::Acronym[pid]), &sov.IsSecondary);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_ReactionID", Particle::Acronym[pid]), &sov.ReactionID);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_IsHybrid", Particle::Acronym[pid]), &sov.IsHybrid);

    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Neg_Entry", Particle::Acronym[pid]), &sov.Neg_Entry);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Neg_Px", Particle::Acronym[pid]), &sov.Neg_Px);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Neg_Py", Particle::Acronym[pid]), &sov.Neg_Py);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Neg_Pz", Particle::Acronym[pid]), &sov.Neg_Pz);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Neg_PdgCode", Particle::Acronym[pid]), &sov.Neg_PdgCode);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Neg_IsTrue", Particle::Acronym[pid]), &sov.Neg_IsTrue);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Neg_IsSignal", Particle::Acronym[pid]), &sov.Neg_IsSignal);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Neg_IsSecondary", Particle::Acronym[pid]), &sov.Neg_IsSecondary);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Neg_ReactionID", Particle::Acronym[pid]), &sov.Neg_ReactionID);

    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Pos_Entry", Particle::Acronym[pid]), &sov.Pos_Entry);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Pos_Px", Particle::Acronym[pid]), &sov.Pos_Px);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Pos_Py", Particle::Acronym[pid]), &sov.Pos_Py);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Pos_Pz", Particle::Acronym[pid]), &sov.Pos_Pz);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Pos_PdgCode", Particle::Acronym[pid]), &sov.Pos_PdgCode);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Pos_IsTrue", Particle::Acronym[pid]), &sov.Pos_IsTrue);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Pos_IsSignal", Particle::Acronym[pid]), &sov.Pos_IsSignal);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Pos_IsSecondary", Particle::Acronym[pid]), &sov.Pos_IsSecondary);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Pos_ReactionID", Particle::Acronym[pid]), &sov.Pos_ReactionID);

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
}

void Packager::CreateOutputBranches_MC_Tracks(EParticle pid, PackedEvents::MC_Tracks& sov) {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Entry", Particle::Acronym[pid]), &sov.Entry);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_X", Particle::Acronym[pid]), &sov.X);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Y", Particle::Acronym[pid]), &sov.Y);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Z", Particle::Acronym[pid]), &sov.Z);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Px", Particle::Acronym[pid]), &sov.Px);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Py", Particle::Acronym[pid]), &sov.Py);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Pz", Particle::Acronym[pid]), &sov.Pz);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_E", Particle::Acronym[pid]), &sov.E);

    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_PdgCode", Particle::Acronym[pid]), &sov.PdgCode);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Mother_Entry", Particle::Acronym[pid]), &sov.Mother_Entry);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_Mother_PdgCode", Particle::Acronym[pid]), &sov.Mother_PdgCode);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_GrandMother_Entry", Particle::Acronym[pid]), &sov.GrandMother_Entry);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_GrandMother_PdgCode", Particle::Acronym[pid]), &sov.GrandMother_PdgCode);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_IsTrue", Particle::Acronym[pid]), &sov.IsTrue);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_IsSignal", Particle::Acronym[pid]), &sov.IsSignal);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_IsSecondary", Particle::Acronym[pid]), &sov.IsSecondary);
    Utils::CreateBranch(fOutputTree.get(), fmt::format("{}_MC_ReactionID", Particle::Acronym[pid]), &sov.ReactionID);

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
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

void Packager::ProcessEvent() {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif
    fOutput_Event.RunNumber = fInput_Event.RunNumber;
    fOutput_Event.DirNumber = fInput_Event.DirNumber;
    if (!IsMC()) fOutput_Event.DirNumberB = fInput_Event.DirNumberB;
    fOutput_Event.EventNumber = fInput_Event.EventNumber;
    fOutput_Event.Centrality = fInput_Event.Centrality;
    fOutput_Event.MagneticField = fInput_Event.MagneticField;
    fOutput_Event.PV_Xv = fInput_Event.PV_Xv;
    fOutput_Event.PV_Yv = fInput_Event.PV_Yv;
    fOutput_Event.PV_Zv = fInput_Event.PV_Zv;

    if (IsMC()) {
        fOutput_Event.MC_PV_Xv = fInput_Event.MC_PV_Xv;
        fOutput_Event.MC_PV_Yv = fInput_Event.MC_PV_Yv;
        fOutput_Event.MC_PV_Zv = fInput_Event.MC_PV_Zv;
    }
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
}

// ## MC/Injected ZONE ## //

void Packager::Injected_GetSecondaryVertex() {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

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

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif
}

void Packager::Injected_Store() {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

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

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
}

// ## Tracks ZONE ## //

// Store tracks' ESD indices into vectors, according to their respective track PID and charge.
void Packager::ProcessTracks() {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    for (int esd_track{0}; esd_track < NumberTracks(); ++esd_track) {
        // get charge //
        int charge{fInput_Tracks.Charge->at(esd_track)};
        // PID //
        if (std::abs(fInput_Tracks.NSigmaProton->at(esd_track)) < Cuts::Track::AbsMax_PID_NSigma) {
            if (charge < 0) fVec_AntiProtons.push_back(esd_track);
            if (charge > 0) fVec_Protons.push_back(esd_track);
        }
        if (std::abs(fInput_Tracks.NSigmaKaon->at(esd_track)) < Cuts::Track::AbsMax_PID_NSigma) {
            if (charge < 0) fVec_NegKaons.push_back(esd_track);
            if (charge > 0) fVec_PosKaons.push_back(esd_track);
        }
        if (std::abs(fInput_Tracks.NSigmaPion->at(esd_track)) < Cuts::Track::AbsMax_PID_NSigma) {
            if (charge < 0) fVec_PiMinus.push_back(esd_track);
            if (charge > 0) fVec_PiPlus.push_back(esd_track);
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

// Note: intended for light particles only, i.e., kaons and pions.
void Packager::PackTracks(EParticle pid) {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    // determine rules based on particle pdg code //
    const std::vector<int>* vec;
    PackedEvents::Tracks* out;
    PackedEvents::MC_Tracks* mc_out{nullptr};
    double mass{Particle::Mass[pid]};
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

    // loop over selected tracks //
    for (auto esd_idx : *vec) {

        // prepare kf object //
        KF::Vector<6> neg_kf_params = KF::PackParams(fInput_Tracks, esd_idx);
        std::array<float, 5> neg_alice_params = ALICE::PackParams<float>(fInput_Tracks, esd_idx);
        std::array<float, 15> neg_alice_cov = ALICE::PackCovMatrix<float>(fInput_Tracks, esd_idx);
        KF::Track kf_track{KF::CreateParticle(neg_kf_params, neg_alice_params, neg_alice_cov, fInput_Tracks.Alpha->at(esd_idx),
                                              fInput_Tracks.Charge->at(esd_idx), mass),
                           esd_idx};

        // store //
        Store(kf_track, *out);
        if (IsMC()) {
            MC::Track mc_track{fInput_MC, fInput_Tracks.McEntry->at(esd_idx), pid};
            StoreMC(mc_track, *mc_out);
        }
    }  // end of loop over selected tracks

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
}

void Packager::Store(const KF::Track& track, PackedEvents::Tracks& sov) {

    sov.Entry->push_back(track.idx);
    sov.X->push_back(static_cast<float>(track.GetParameter(0)));
    sov.Y->push_back(static_cast<float>(track.GetParameter(1)));
    sov.Z->push_back(static_cast<float>(track.GetParameter(2)));
    sov.Px->push_back(static_cast<float>(track.GetParameter(3)));
    sov.Py->push_back(static_cast<float>(track.GetParameter(4)));
    sov.Pz->push_back(static_cast<float>(track.GetParameter(5)));
    sov.E->push_back(static_cast<float>(track.GetParameter(6)));

    sov.Sigma.X2->push_back(static_cast<float>(track.GetCovariance(0)));
    sov.Sigma.XY->push_back(static_cast<float>(track.GetCovariance(1)));
    sov.Sigma.Y2->push_back(static_cast<float>(track.GetCovariance(2)));
    sov.Sigma.XZ->push_back(static_cast<float>(track.GetCovariance(3)));
    sov.Sigma.YZ->push_back(static_cast<float>(track.GetCovariance(4)));
    sov.Sigma.Z2->push_back(static_cast<float>(track.GetCovariance(5)));
    sov.Sigma.XPx->push_back(static_cast<float>(track.GetCovariance(6)));
    sov.Sigma.YPx->push_back(static_cast<float>(track.GetCovariance(7)));
    sov.Sigma.ZPx->push_back(static_cast<float>(track.GetCovariance(8)));
    sov.Sigma.Px2->push_back(static_cast<float>(track.GetCovariance(9)));
    sov.Sigma.XPy->push_back(static_cast<float>(track.GetCovariance(10)));
    sov.Sigma.YPy->push_back(static_cast<float>(track.GetCovariance(11)));
    sov.Sigma.ZPy->push_back(static_cast<float>(track.GetCovariance(12)));
    sov.Sigma.PxPy->push_back(static_cast<float>(track.GetCovariance(13)));
    sov.Sigma.Py2->push_back(static_cast<float>(track.GetCovariance(14)));
    sov.Sigma.XPz->push_back(static_cast<float>(track.GetCovariance(15)));
    sov.Sigma.YPz->push_back(static_cast<float>(track.GetCovariance(16)));
    sov.Sigma.ZPz->push_back(static_cast<float>(track.GetCovariance(17)));
    sov.Sigma.PxPz->push_back(static_cast<float>(track.GetCovariance(18)));
    sov.Sigma.PyPz->push_back(static_cast<float>(track.GetCovariance(19)));
    sov.Sigma.Pz2->push_back(static_cast<float>(track.GetCovariance(20)));
    sov.Sigma.XE->push_back(static_cast<float>(track.GetCovariance(21)));
    sov.Sigma.YE->push_back(static_cast<float>(track.GetCovariance(22)));
    sov.Sigma.ZE->push_back(static_cast<float>(track.GetCovariance(23)));
    sov.Sigma.PxE->push_back(static_cast<float>(track.GetCovariance(24)));
    sov.Sigma.PyE->push_back(static_cast<float>(track.GetCovariance(25)));
    sov.Sigma.PzE->push_back(static_cast<float>(track.GetCovariance(26)));
    sov.Sigma.E2->push_back(static_cast<float>(track.GetCovariance(27)));
}

void Packager::StoreMC(const MC::Track& mc_track, PackedEvents::MC_Tracks& sov) {

    sov.Entry->push_back(mc_track.Entry);
    sov.X->push_back(mc_track.X);
    sov.Y->push_back(mc_track.Y);
    sov.Z->push_back(mc_track.Z);
    sov.Px->push_back(mc_track.Px);
    sov.Py->push_back(mc_track.Py);
    sov.Pz->push_back(mc_track.Pz);
    sov.E->push_back(mc_track.Energy);

    sov.PdgCode->push_back(mc_track.PdgCode);
    sov.Mother_Entry->push_back(mc_track.Mother_Entry);
    sov.Mother_PdgCode->push_back(mc_track.Mother_PdgCode);
    sov.GrandMother_Entry->push_back(mc_track.GrandMother_Entry);
    sov.GrandMother_PdgCode->push_back(mc_track.GrandMother_PdgCode);
    sov.IsTrue->push_back(static_cast<char>(mc_track.IsTrue));
    sov.IsSignal->push_back(static_cast<char>(mc_track.IsSignal));
    sov.IsSecondary->push_back(static_cast<char>(mc_track.IsSecondary));
    sov.ReactionID->push_back(mc_track.ReactionID);
}

// ## V0s ZONE ## //

void Packager::FindV0s(EParticle pid) {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    // determine rules based on V0 pdg code //
    const std::vector<int>* vec_neg;
    const std::vector<int>* vec_pos;
    PackedEvents::V0s* out;
    PackedEvents::MC_V0s* mc_out{nullptr};
    EParticle pid_neg;
    EParticle pid_pos;
    switch (pid) {
        case EParticle::AntiLambda:
            vec_neg = &fVec_AntiProtons;
            vec_pos = &fVec_PiPlus;
            out = &fOutput_AntiLambdas;
            if (IsMC()) mc_out = &fOutput_MC_AntiLambdas;
            pid_neg = EParticle::AntiProton;
            pid_pos = EParticle::PiPlus;
            break;
        case EParticle::Lambda:
            vec_neg = &fVec_PiMinus;
            vec_pos = &fVec_Protons;
            out = &fOutput_Lambdas;
            if (IsMC()) mc_out = &fOutput_MC_Lambdas;
            pid_neg = EParticle::PiMinus;
            pid_pos = EParticle::Proton;
            break;
        case EParticle::KaonZeroShort:
            vec_neg = &fVec_PiMinus;
            vec_pos = &fVec_PiPlus;
            out = &fOutput_KaonsZeroShort;
            if (IsMC()) mc_out = &fOutput_MC_KaonsZeroShort;
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

#ifdef T2S_USE_ALICE
            // ----------------------------- //
            // start of **Custom V0 Finder** //

            std::array<double, 5> neg_esd_params = ALICE::PackParams<double>(fInput_Tracks, esd_neg);
            std::array<double, 15> neg_esd_cov = ALICE::PackCovMatrix<double>(fInput_Tracks, esd_neg);
            ALICE::Track neg{fInput_Tracks.X->at(esd_neg), neg_esd_params, neg_esd_cov, fInput_Tracks.Alpha->at(esd_neg), -1};

            std::array<double, 5> pos_esd_params = ALICE::PackParams<double>(fInput_Tracks, esd_pos);
            std::array<double, 15> pos_esd_cov = ALICE::PackCovMatrix<double>(fInput_Tracks, esd_pos);
            ALICE::Track pos{fInput_Tracks.X->at(esd_pos), pos_esd_params, pos_esd_cov, fInput_Tracks.Alpha->at(esd_pos), +1};

            // double neg_d_wrt_pv{neg.GetDCAxy(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.MagneticField)};
            // if (TMath::Abs(neg_d_wrt_pv) > ALICE::Const::CustomV0Finder_Rmax) continue;
            // double pos_d_wrt_pv{pos.GetDCAxy(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.MagneticField)};
            // if (TMath::Abs(pos_d_wrt_pv) > ALICE::Const::CustomV0Finder_Rmax) continue;

            double xn;
            double xp;
            if (ALICE::Vertexer::Preoptimize(neg, pos, xn, xp, fInput_Event.MagneticField)) {
                if (!neg.PropagateTo(xn, fInput_Event.MagneticField)) continue;
                if (!pos.PropagateTo(xp, fInput_Event.MagneticField)) continue;
            } else {
                double dca{ALICE::Vertexer::Preoptimize_Numerically(neg, pos, xn, xp, fInput_Event.MagneticField)};
                if (!neg.PropagateTo(xn, fInput_Event.MagneticField)) continue;
                if (!pos.PropagateTo(xp, fInput_Event.MagneticField)) continue;
            }
#ifdef T2S_DEBUG
            Logger::Info(__FUNCTION__, "DEBUG NEG: Snp=" << neg.GetSnp() << " SigmaY2=" << neg.GetSigmaY2());
            Logger::Info(__FUNCTION__, "DEBUG POS: Snp=" << pos.GetSnp() << " SigmaY2=" << pos.GetSigmaY2());
#endif

            // if (dca_after_preopt > ALICE::Const::CustomV0Finder_DCAmax) continue;
            // if ((xn + xp) > 2. * ALICE::Const::CustomV0Finder_Rmax) continue;
            // if ((xn + xp) < 2. * ALICE::Const::CustomV0Finder_Rmin) continue;

            ALICE::V0 v0{neg, pos};
            v0.Refit();
            // if (v0.fStatus == 1) continue;

            // end of **Custom V0 Finder** //
            // --------------------------- //
#ifdef T2S_DEBUG
            Logger::Info(__FUNCTION__,
                         "   "
                         " :: idx,neg,pos="
                             << v0_entry << "," << esd_neg << "," << esd_pos);
            Logger::Info(__FUNCTION__, ";x,y,z=" << v0.X() << "," << v0.Y() << "," << v0.Z());
            // Logger::Info(__FUNCTION__,  ";x,y,z(neg)=" << v0.Neg_PCA_XYZ()[0] << "," << v0.Neg_PCA_XYZ()[1] << "," << v0.Neg_PCA_XYZ()[2]);
            // Logger::Info(__FUNCTION__,  ";x,y,z(pos)=" << v0.Pos_PCA_XYZ()[0] << "," << v0.Pos_PCA_XYZ()[1] << "," << v0.Pos_PCA_XYZ()[2]);
            // Logger::Info(__FUNCTION__,  ";mass=" << v0.Mass());
            // Logger::Info(__FUNCTION__,  ";dca_dau=" << v0.DCA_Daughters());
            Logger::Info(__FUNCTION__, ";radius=" << v0.Radius2D());
            // Logger::Info(__FUNCTION__,  ";dca_neg=" << v0.DCA_Neg_V0());
            // Logger::Info(__FUNCTION__,  ";dca_pos=" << v0.DCA_Pos_V0());
            // Logger::Info(__FUNCTION__,  ";pt=" << v0.Pt());
            // Logger::Info(__FUNCTION__,  ";eta=" << v0.Eta());
            // Logger::Info(__FUNCTION__,  ";qt=" << v0.ArmenterosQt());
            // Logger::Info(__FUNCTION__,  ";alpha=" << v0.ArmenterosAlpha());
            // Logger::Info(__FUNCTION__,  ";cpa_pv=" << v0.CPA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv));
            // Logger::Info(__FUNCTION__,  ";dca_pv=" << v0.DCA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv));
            Logger::Info(__FUNCTION__, '\n');
#endif
            // store //
            Store(v0, *out);
#else
            // prepare neg //
            std::array<double, 6> neg_kf_params = KF::PackParams(fInput_Tracks, esd_neg);
            std::array<float, 5> neg_alice_params = ALICE::PackParams<float>(fInput_Tracks, esd_neg);
            std::array<float, 15> neg_alice_cov = ALICE::PackCovMatrix<float>(fInput_Tracks, esd_neg);
            KF::Track neg{KF::CreateParticle(neg_kf_params, neg_alice_params, neg_alice_cov, fInput_Tracks.Alpha->at(esd_neg),
                                             fInput_Tracks.Charge->at(esd_neg), mass_neg),
                          esd_neg};

            // prepare pos //
            std::array<double, 6> pos_kf_params = KF::PackParams(fInput_Tracks, esd_pos);
            std::array<float, 5> pos_alice_params = ALICE::PackParams<float>(fInput_Tracks, esd_pos);
            std::array<float, 15> pos_alice_cov = ALICE::PackCovMatrix<float>(fInput_Tracks, esd_pos);
            KF::Track pos{KF::CreateParticle(pos_kf_params, pos_alice_params, pos_alice_cov, fInput_Tracks.Alpha->at(esd_pos),
                                             fInput_Tracks.Charge->at(esd_pos), mass_pos),
                          esd_pos};

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
            Logger::Debug(__FUNCTION__, ";cpa_pv={}", v0.CPA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv));
            Logger::Debug(__FUNCTION__, ";dca_pv={}", v0.DCA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv));
            Logger::Debug(__FUNCTION__, "");
#endif
            // store //
            Store(v0, *out);
#endif

            if (IsMC()) {
                MC::V0 mc_v0{fInput_MC, fInput_Tracks.McEntry->at(esd_neg), fInput_Tracks.McEntry->at(esd_pos), pid, pid_neg, pid_pos};
                StoreMC(mc_v0, *mc_out);
            }
            ++v0_entry;
        }  // end of loop over pos
    }  // end of loop over neg

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
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
    // if (v0.AbsArmQtOverAlpha() > Cuts::Lambda::AbsMax_ArmQtOverAlpha) return false;
    CutFlowHist->Fill(9.);
    if (v0.CPA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv) < Cuts::Lambda::Min_CPAwrtPV ||
        v0.CPA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv) > Cuts::Lambda::Max_CPAwrtPV) {
        return false;
    }
    CutFlowHist->Fill(10.);
    if (v0.DCA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv) < Cuts::Lambda::Min_DCAwrtPV) return false;
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
    if (v0.CPA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv) < Cuts::KaonZeroShort::Min_CPAwrtPV ||
        v0.CPA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv) > Cuts::KaonZeroShort::Max_CPAwrtPV) {
        return false;
    }
    fCutFlowHist_KaonsZeroShort->Fill(9.);
    if (v0.DCA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv) < Cuts::KaonZeroShort::Min_DCAwrtPV) return false;
    fCutFlowHist_KaonsZeroShort->Fill(10.);

    return true;
}

void Packager::Store(const ALICE::V0& v0, PackedEvents::V0s& sov) {

    // sov.Entry->push_back(v0.idx);
    sov.X->push_back(static_cast<float>(v0.X()));
    sov.Y->push_back(static_cast<float>(v0.Y()));
    sov.Z->push_back(static_cast<float>(v0.Z()));
    // sov.Px->push_back(static_cast<float>(v0.Px()));
    // sov.Py->push_back(static_cast<float>(v0.Py()));
    // sov.Pz->push_back(static_cast<float>(v0.Pz()));
    // sov.E->push_back(static_cast<float>(v0.E()));

    // sov.Sigma.X2->push_back(static_cast<float>(v0.GetCovariance(0)));
    // sov.Sigma.XY->push_back(static_cast<float>(v0.GetCovariance(1)));
    // sov.Sigma.Y2->push_back(static_cast<float>(v0.GetCovariance(2)));
    // sov.Sigma.XZ->push_back(static_cast<float>(v0.GetCovariance(3)));
    // sov.Sigma.YZ->push_back(static_cast<float>(v0.GetCovariance(4)));
    // sov.Sigma.Z2->push_back(static_cast<float>(v0.GetCovariance(5)));
    // sov.Sigma.XPx->push_back(static_cast<float>(v0.GetCovariance(6)));
    // sov.Sigma.YPx->push_back(static_cast<float>(v0.GetCovariance(7)));
    // sov.Sigma.ZPx->push_back(static_cast<float>(v0.GetCovariance(8)));
    // sov.Sigma.Px2->push_back(static_cast<float>(v0.GetCovariance(9)));
    // sov.Sigma.XPy->push_back(static_cast<float>(v0.GetCovariance(10)));
    // sov.Sigma.YPy->push_back(static_cast<float>(v0.GetCovariance(11)));
    // sov.Sigma.ZPy->push_back(static_cast<float>(v0.GetCovariance(12)));
    // sov.Sigma.PxPy->push_back(static_cast<float>(v0.GetCovariance(13)));
    // sov.Sigma.Py2->push_back(static_cast<float>(v0.GetCovariance(14)));
    // sov.Sigma.XPz->push_back(static_cast<float>(v0.GetCovariance(15)));
    // sov.Sigma.YPz->push_back(static_cast<float>(v0.GetCovariance(16)));
    // sov.Sigma.ZPz->push_back(static_cast<float>(v0.GetCovariance(17)));
    // sov.Sigma.PxPz->push_back(static_cast<float>(v0.GetCovariance(18)));
    // sov.Sigma.PyPz->push_back(static_cast<float>(v0.GetCovariance(19)));
    // sov.Sigma.Pz2->push_back(static_cast<float>(v0.GetCovariance(20)));
    // sov.Sigma.XE->push_back(static_cast<float>(v0.GetCovariance(21)));
    // sov.Sigma.YE->push_back(static_cast<float>(v0.GetCovariance(22)));
    // sov.Sigma.ZE->push_back(static_cast<float>(v0.GetCovariance(23)));
    // sov.Sigma.PxE->push_back(static_cast<float>(v0.GetCovariance(24)));
    // sov.Sigma.PyE->push_back(static_cast<float>(v0.GetCovariance(25)));
    // sov.Sigma.PzE->push_back(static_cast<float>(v0.GetCovariance(26)));
    // sov.Sigma.E2->push_back(static_cast<float>(v0.GetCovariance(27)));

    // sov.Neg.Entry->push_back(v0.Neg.idx);
    // sov.Neg.X->push_back(static_cast<float>(v0.Neg.X()));
    // sov.Neg.Y->push_back(static_cast<float>(v0.Neg.Y()));
    // sov.Neg.Z->push_back(static_cast<float>(v0.Neg.Z()));
    // sov.Neg.Px->push_back(static_cast<float>(v0.Neg.Px()));
    // sov.Neg.Py->push_back(static_cast<float>(v0.Neg.Py()));
    // sov.Neg.Pz->push_back(static_cast<float>(v0.Neg.Pz()));
    // sov.Neg.E->push_back(static_cast<float>(v0.Neg.E()));

    sov.Neg_X_AtPCA->push_back(static_cast<float>(v0.fParamN.GetX()));
    sov.Neg_Y_AtPCA->push_back(static_cast<float>(v0.fParamN.GetY()));
    sov.Neg_Z_AtPCA->push_back(static_cast<float>(v0.fParamN.GetZ()));
    // sov.Neg_Px_AtPCA->push_back(static_cast<float>(v0.Neg_PCA_PxPyPz()[0]));
    // sov.Neg_Py_AtPCA->push_back(static_cast<float>(v0.Neg_PCA_PxPyPz()[1]));
    // sov.Neg_Pz_AtPCA->push_back(static_cast<float>(v0.Neg_PCA_PxPyPz()[2]));

    // sov.Pos.Entry->push_back(v0.Pos.idx);
    // sov.Pos.X->push_back(static_cast<float>(v0.Pos.X()));
    // sov.Pos.Y->push_back(static_cast<float>(v0.Pos.Y()));
    // sov.Pos.Z->push_back(static_cast<float>(v0.Pos.Z()));
    // sov.Pos.Px->push_back(static_cast<float>(v0.Pos.Px()));
    // sov.Pos.Py->push_back(static_cast<float>(v0.Pos.Py()));
    // sov.Pos.Pz->push_back(static_cast<float>(v0.Pos.Pz()));
    // sov.Pos.E->push_back(static_cast<float>(v0.Pos.E()));

    sov.Pos_X_AtPCA->push_back(static_cast<float>(v0.fParamP.GetX()));
    sov.Pos_Y_AtPCA->push_back(static_cast<float>(v0.fParamP.GetY()));
    sov.Pos_Z_AtPCA->push_back(static_cast<float>(v0.fParamP.GetZ()));
    // sov.Pos_Px_AtPCA->push_back(static_cast<float>(v0.Pos_PCA_PxPyPz()[0]));
    // sov.Pos_Py_AtPCA->push_back(static_cast<float>(v0.Pos_PCA_PxPyPz()[1]));
    // sov.Pos_Pz_AtPCA->push_back(static_cast<float>(v0.Pos_PCA_PxPyPz()[2]));
}

void Packager::Store(const KF::V0& v0, PackedEvents::V0s& sov) {

    sov.Entry->push_back(v0.idx);
    sov.X->push_back(static_cast<float>(v0.X()));
    sov.Y->push_back(static_cast<float>(v0.Y()));
    sov.Z->push_back(static_cast<float>(v0.Z()));
    sov.Px->push_back(static_cast<float>(v0.Px()));
    sov.Py->push_back(static_cast<float>(v0.Py()));
    sov.Pz->push_back(static_cast<float>(v0.Pz()));
    sov.E->push_back(static_cast<float>(v0.E()));

    sov.Sigma.X2->push_back(static_cast<float>(v0.GetCovariance(0)));
    sov.Sigma.XY->push_back(static_cast<float>(v0.GetCovariance(1)));
    sov.Sigma.Y2->push_back(static_cast<float>(v0.GetCovariance(2)));
    sov.Sigma.XZ->push_back(static_cast<float>(v0.GetCovariance(3)));
    sov.Sigma.YZ->push_back(static_cast<float>(v0.GetCovariance(4)));
    sov.Sigma.Z2->push_back(static_cast<float>(v0.GetCovariance(5)));
    sov.Sigma.XPx->push_back(static_cast<float>(v0.GetCovariance(6)));
    sov.Sigma.YPx->push_back(static_cast<float>(v0.GetCovariance(7)));
    sov.Sigma.ZPx->push_back(static_cast<float>(v0.GetCovariance(8)));
    sov.Sigma.Px2->push_back(static_cast<float>(v0.GetCovariance(9)));
    sov.Sigma.XPy->push_back(static_cast<float>(v0.GetCovariance(10)));
    sov.Sigma.YPy->push_back(static_cast<float>(v0.GetCovariance(11)));
    sov.Sigma.ZPy->push_back(static_cast<float>(v0.GetCovariance(12)));
    sov.Sigma.PxPy->push_back(static_cast<float>(v0.GetCovariance(13)));
    sov.Sigma.Py2->push_back(static_cast<float>(v0.GetCovariance(14)));
    sov.Sigma.XPz->push_back(static_cast<float>(v0.GetCovariance(15)));
    sov.Sigma.YPz->push_back(static_cast<float>(v0.GetCovariance(16)));
    sov.Sigma.ZPz->push_back(static_cast<float>(v0.GetCovariance(17)));
    sov.Sigma.PxPz->push_back(static_cast<float>(v0.GetCovariance(18)));
    sov.Sigma.PyPz->push_back(static_cast<float>(v0.GetCovariance(19)));
    sov.Sigma.Pz2->push_back(static_cast<float>(v0.GetCovariance(20)));
    sov.Sigma.XE->push_back(static_cast<float>(v0.GetCovariance(21)));
    sov.Sigma.YE->push_back(static_cast<float>(v0.GetCovariance(22)));
    sov.Sigma.ZE->push_back(static_cast<float>(v0.GetCovariance(23)));
    sov.Sigma.PxE->push_back(static_cast<float>(v0.GetCovariance(24)));
    sov.Sigma.PyE->push_back(static_cast<float>(v0.GetCovariance(25)));
    sov.Sigma.PzE->push_back(static_cast<float>(v0.GetCovariance(26)));
    sov.Sigma.E2->push_back(static_cast<float>(v0.GetCovariance(27)));

    sov.Neg.Entry->push_back(v0.Neg.idx);
    sov.Neg.X->push_back(static_cast<float>(v0.Neg.X()));
    sov.Neg.Y->push_back(static_cast<float>(v0.Neg.Y()));
    sov.Neg.Z->push_back(static_cast<float>(v0.Neg.Z()));
    sov.Neg.Px->push_back(static_cast<float>(v0.Neg.Px()));
    sov.Neg.Py->push_back(static_cast<float>(v0.Neg.Py()));
    sov.Neg.Pz->push_back(static_cast<float>(v0.Neg.Pz()));
    sov.Neg.E->push_back(static_cast<float>(v0.Neg.E()));

    sov.Neg_X_AtPCA->push_back(static_cast<float>(v0.Neg_PCA_XYZ()[0]));
    sov.Neg_Y_AtPCA->push_back(static_cast<float>(v0.Neg_PCA_XYZ()[1]));
    sov.Neg_Z_AtPCA->push_back(static_cast<float>(v0.Neg_PCA_XYZ()[2]));
    sov.Neg_Px_AtPCA->push_back(static_cast<float>(v0.Neg_PCA_PxPyPz()[0]));
    sov.Neg_Py_AtPCA->push_back(static_cast<float>(v0.Neg_PCA_PxPyPz()[1]));
    sov.Neg_Pz_AtPCA->push_back(static_cast<float>(v0.Neg_PCA_PxPyPz()[2]));

    sov.Pos.Entry->push_back(v0.Pos.idx);
    sov.Pos.X->push_back(static_cast<float>(v0.Pos.X()));
    sov.Pos.Y->push_back(static_cast<float>(v0.Pos.Y()));
    sov.Pos.Z->push_back(static_cast<float>(v0.Pos.Z()));
    sov.Pos.Px->push_back(static_cast<float>(v0.Pos.Px()));
    sov.Pos.Py->push_back(static_cast<float>(v0.Pos.Py()));
    sov.Pos.Pz->push_back(static_cast<float>(v0.Pos.Pz()));
    sov.Pos.E->push_back(static_cast<float>(v0.Pos.E()));

    sov.Pos_X_AtPCA->push_back(static_cast<float>(v0.Pos_PCA_XYZ()[0]));
    sov.Pos_Y_AtPCA->push_back(static_cast<float>(v0.Pos_PCA_XYZ()[1]));
    sov.Pos_Z_AtPCA->push_back(static_cast<float>(v0.Pos_PCA_XYZ()[2]));
    sov.Pos_Px_AtPCA->push_back(static_cast<float>(v0.Pos_PCA_PxPyPz()[0]));
    sov.Pos_Py_AtPCA->push_back(static_cast<float>(v0.Pos_PCA_PxPyPz()[1]));
    sov.Pos_Pz_AtPCA->push_back(static_cast<float>(v0.Pos_PCA_PxPyPz()[2]));
}

void Packager::StoreMC(const MC::V0& mc_v0, PackedEvents::MC_V0s& sov) {

    sov.Entry->push_back(mc_v0.Entry);
    sov.X->push_back(mc_v0.X);
    sov.Y->push_back(mc_v0.Y);
    sov.Z->push_back(mc_v0.Z);
    sov.Px->push_back(mc_v0.Px);
    sov.Py->push_back(mc_v0.Py);
    sov.Pz->push_back(mc_v0.Pz);
    sov.E->push_back(mc_v0.Energy);

    sov.DecayX->push_back(mc_v0.DecayX());
    sov.DecayY->push_back(mc_v0.DecayY());
    sov.DecayZ->push_back(mc_v0.DecayZ());

    sov.PdgCode->push_back(mc_v0.PdgCode);
    sov.Mother_Entry->push_back(mc_v0.Mother_Entry);
    sov.Mother_PdgCode->push_back(mc_v0.Mother_PdgCode);
    sov.IsTrue->push_back(static_cast<char>(mc_v0.IsTrue));
    sov.IsSignal->push_back(static_cast<char>(mc_v0.IsSignal));
    sov.IsSecondary->push_back(static_cast<char>(mc_v0.IsSecondary));
    sov.ReactionID->push_back(mc_v0.ReactionID);
    sov.IsHybrid->push_back(static_cast<char>(mc_v0.IsHybrid));

    // neg //
    sov.Neg_Entry->push_back(mc_v0.neg.Entry);
    sov.Neg_Px->push_back(mc_v0.neg.Px);
    sov.Neg_Py->push_back(mc_v0.neg.Py);
    sov.Neg_Pz->push_back(mc_v0.neg.Pz);
    sov.Neg_PdgCode->push_back(mc_v0.neg.PdgCode);
    sov.Neg_IsTrue->push_back(static_cast<char>(mc_v0.neg.IsTrue));
    sov.Neg_IsSignal->push_back(static_cast<char>(mc_v0.neg.IsSignal));
    sov.Neg_IsSecondary->push_back(static_cast<char>(mc_v0.neg.IsSecondary));
    sov.Neg_ReactionID->push_back(mc_v0.neg.ReactionID);

    // pos //
    sov.Pos_Entry->push_back(mc_v0.pos.Entry);
    sov.Pos_Px->push_back(mc_v0.pos.Px);
    sov.Pos_Py->push_back(mc_v0.pos.Py);
    sov.Pos_Pz->push_back(mc_v0.pos.Pz);
    sov.Pos_PdgCode->push_back(mc_v0.pos.PdgCode);
    sov.Pos_IsTrue->push_back(static_cast<char>(mc_v0.pos.IsTrue));
    sov.Pos_IsSignal->push_back(static_cast<char>(mc_v0.pos.IsSignal));
    sov.Pos_IsSecondary->push_back(static_cast<char>(mc_v0.pos.IsSecondary));
    sov.Pos_ReactionID->push_back(mc_v0.pos.ReactionID);
}

// ## END OF CYCLES ## //

void Packager::EndOfEvent() {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

    // fill tree //
    fOutputTree->Fill();

    // clear temporary containers //
    fVec_SV_X.clear();
    fVec_SV_Y.clear();
    fVec_SV_Z.clear();

    fVec_AntiProtons.clear();
    fVec_Protons.clear();
    fVec_NegKaons.clear();
    fVec_PosKaons.clear();
    fVec_PiMinus.clear();
    fVec_PiPlus.clear();

    // clear output branches //
    if (IsMC()) fOutput_Injected.Clear();
    switch (GetReactionChannel()) {
        // standard channels //
        case EReactionChannel::A:
            fOutput_AntiLambdas.Clear();
            fOutput_KaonsZeroShort.Clear();
            if (IsMC()) {
                fOutput_MC_AntiLambdas.Clear();
                fOutput_MC_KaonsZeroShort.Clear();
            }
            break;
        case EReactionChannel::D:
            fOutput_AntiLambdas.Clear();
            fOutput_PosKaons.Clear();
            if (IsMC()) {
                fOutput_MC_AntiLambdas.Clear();
                fOutput_MC_PosKaons.Clear();
            }
            break;
        case EReactionChannel::E:
            fOutput_AntiLambdas.Clear();
            fOutput_PosKaons.Clear();
            fOutput_PiMinus.Clear();
            fOutput_PiPlus.Clear();
            if (IsMC()) {
                fOutput_MC_AntiLambdas.Clear();
                fOutput_MC_PosKaons.Clear();
                fOutput_MC_PiMinus.Clear();
                fOutput_MC_PiPlus.Clear();
            }
            break;
        case EReactionChannel::H:
            fOutput_PosKaons.Clear();
            if (IsMC()) fOutput_MC_PosKaons.Clear();
            break;
        // anti-channels //
        case EReactionChannel::AntiA:
            fOutput_Lambdas.Clear();
            fOutput_KaonsZeroShort.Clear();
            if (IsMC()) {
                fOutput_MC_Lambdas.Clear();
                fOutput_MC_KaonsZeroShort.Clear();
            }
            break;
        case EReactionChannel::AntiD:
            fOutput_Lambdas.Clear();
            fOutput_NegKaons.Clear();
            if (IsMC()) {
                fOutput_MC_Lambdas.Clear();
                fOutput_MC_NegKaons.Clear();
            }
            break;
        case EReactionChannel::AntiE:
            fOutput_Lambdas.Clear();
            fOutput_NegKaons.Clear();
            fOutput_PiMinus.Clear();
            fOutput_PiPlus.Clear();
            if (IsMC()) {
                fOutput_MC_Lambdas.Clear();
                fOutput_MC_NegKaons.Clear();
                fOutput_MC_PiMinus.Clear();
                fOutput_MC_PiPlus.Clear();
            }
            break;
        case EReactionChannel::AntiH:
            fOutput_NegKaons.Clear();
            if (IsMC()) fOutput_MC_NegKaons.Clear();
            break;
        // for data //
        case EReactionChannel::All:
            fOutput_AntiLambdas.Clear();
            fOutput_Lambdas.Clear();
            fOutput_KaonsZeroShort.Clear();
            fOutput_NegKaons.Clear();
            fOutput_PosKaons.Clear();
            fOutput_PiMinus.Clear();
            fOutput_PiPlus.Clear();
            if (IsMC()) {
                fOutput_MC_AntiLambdas.Clear();
                fOutput_MC_Lambdas.Clear();
                fOutput_MC_KaonsZeroShort.Clear();
                fOutput_MC_NegKaons.Clear();
                fOutput_MC_PosKaons.Clear();
                fOutput_MC_PiMinus.Clear();
                fOutput_MC_PiPlus.Clear();
            }
            break;
    }

#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Finished.");
#endif
}

void Packager::EndOfAnalysis() {
#ifdef T2S_DEBUG
    Logger::Debug(__FUNCTION__, "Starting.");
#endif

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

    fEventsTree->ResetBranchAddresses();
    fOutputTree->ResetBranchAddresses();

    Logger::Info(__FUNCTION__, "All done.");
}

}  // namespace Tree2Secondaries
