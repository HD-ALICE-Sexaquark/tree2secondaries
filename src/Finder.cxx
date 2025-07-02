#include <filesystem>
#include <iostream>
#include <memory>
#include <set>

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
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    fTree_PackedEvents = std::make_unique<TChain>("PackedEvents");
    for (const auto& path : fSettings.PathInputFiles) {
        if (fTree_PackedEvents->Add(path.c_str()) == 0) {
            std::cerr << "   " << __FUNCTION__ << " :: Couldn't add TFile \"" << path << "\"" << '\n';
        }
    }
    if (!fTree_PackedEvents->GetEntries()) {
        std::cerr << "   " << __FUNCTION__ << " :: Couldn't manage to read any entry." << '\n';
        return false;
    }
    std::cout << "   " << __FUNCTION__ << " :: TChain \"" << fTree_PackedEvents->GetName() << "\" loaded successfully with "
              << fTree_PackedEvents->GetNtrees() << " trees and " << fTree_PackedEvents->GetEntries() << " total entries." << '\n';

    ConnectInputBranches();

    if (!PrepareOutputFile()) return false;

    if (!PrepareOutputTree()) return false;
    CreateOutputBranches();

    if (IsSignalMC()) {
        if (!Injected_PrepareOutputTree()) return false;
        Injected_CreateOutputBranches();
    }

    std::cout << "   " << __FUNCTION__ << " :: Finder initialized successfully" << '\n';

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
    return true;
}

// ## INPUT ZONE ## //

void Finder::ConnectInputBranches() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    fTree_PackedEvents->SetBranchStatus("*", false);

    ConnectBranches_Events();
    if (IsMC() && IsSignalMC()) ConnectBranches_Injected();

    switch (GetReactionChannel()) {
        // standard channels //
        case ReactionChannel::A:
            ConnectBranches_V0s(static_cast<std::string>(Acronym::AntiLambda), fPacked_AntiLambdas);
            ConnectBranches_V0s(static_cast<std::string>(Acronym::KaonZeroShort), fPacked_KaonsZeroShort);
            if (IsMC()) {
                ConnectBranches_MC_V0s(static_cast<std::string>(Acronym::AntiLambda), fPacked_MC_AntiLambdas);
                ConnectBranches_MC_V0s(static_cast<std::string>(Acronym::KaonZeroShort), fPacked_MC_KaonsZeroShort);
            }
            break;
        case ReactionChannel::D:
            ConnectBranches_V0s(static_cast<std::string>(Acronym::AntiLambda), fPacked_AntiLambdas);
            ConnectBranches_Tracks(static_cast<std::string>(Acronym::PosKaon), fPacked_PosKaons);
            if (IsMC()) {
                ConnectBranches_MC_V0s(static_cast<std::string>(Acronym::AntiLambda), fPacked_MC_AntiLambdas);
                ConnectBranches_MC_Tracks(static_cast<std::string>(Acronym::PosKaon), fPacked_MC_PosKaons);
            }
            break;
        case ReactionChannel::E:
            ConnectBranches_V0s(static_cast<std::string>(Acronym::AntiLambda), fPacked_AntiLambdas);
            ConnectBranches_Tracks(static_cast<std::string>(Acronym::PosKaon), fPacked_PosKaons);
            ConnectBranches_Tracks(static_cast<std::string>(Acronym::PiMinus), fPacked_PiMinus);
            ConnectBranches_Tracks(static_cast<std::string>(Acronym::PiPlus), fPacked_PiPlus);
            if (IsMC()) {
                ConnectBranches_MC_V0s(static_cast<std::string>(Acronym::AntiLambda), fPacked_MC_AntiLambdas);
                ConnectBranches_MC_Tracks(static_cast<std::string>(Acronym::PosKaon), fPacked_MC_PosKaons);
                ConnectBranches_MC_Tracks(static_cast<std::string>(Acronym::PiMinus), fPacked_MC_PiMinus);
                ConnectBranches_MC_Tracks(static_cast<std::string>(Acronym::PiPlus), fPacked_MC_PiPlus);
            }
            break;
        case ReactionChannel::H:
            ConnectBranches_Tracks(static_cast<std::string>(Acronym::PosKaon), fPacked_PosKaons);
            if (IsMC()) ConnectBranches_MC_Tracks(static_cast<std::string>(Acronym::PosKaon), fPacked_MC_PosKaons);
            break;
        // anti-channels //
        case ReactionChannel::AntiA:
            ConnectBranches_V0s(static_cast<std::string>(Acronym::Lambda), fPacked_Lambdas);
            ConnectBranches_V0s(static_cast<std::string>(Acronym::KaonZeroShort), fPacked_KaonsZeroShort);
            if (IsMC()) {
                ConnectBranches_MC_V0s(static_cast<std::string>(Acronym::Lambda), fPacked_MC_Lambdas);
                ConnectBranches_MC_V0s(static_cast<std::string>(Acronym::KaonZeroShort), fPacked_MC_KaonsZeroShort);
            }
            break;
        case ReactionChannel::AntiD:
            ConnectBranches_V0s(static_cast<std::string>(Acronym::Lambda), fPacked_Lambdas);
            ConnectBranches_Tracks(static_cast<std::string>(Acronym::NegKaon), fPacked_NegKaons);
            if (IsMC()) {
                ConnectBranches_MC_V0s(static_cast<std::string>(Acronym::Lambda), fPacked_MC_Lambdas);
                ConnectBranches_MC_Tracks(static_cast<std::string>(Acronym::NegKaon), fPacked_MC_NegKaons);
            }
            break;
        case ReactionChannel::AntiE:
            ConnectBranches_V0s(static_cast<std::string>(Acronym::Lambda), fPacked_Lambdas);
            ConnectBranches_Tracks(static_cast<std::string>(Acronym::NegKaon), fPacked_NegKaons);
            ConnectBranches_Tracks(static_cast<std::string>(Acronym::PiMinus), fPacked_PiMinus);
            ConnectBranches_Tracks(static_cast<std::string>(Acronym::PiPlus), fPacked_PiPlus);
            if (IsMC()) {
                ConnectBranches_MC_V0s(static_cast<std::string>(Acronym::Lambda), fPacked_MC_Lambdas);
                ConnectBranches_MC_Tracks(static_cast<std::string>(Acronym::NegKaon), fPacked_MC_NegKaons);
                ConnectBranches_MC_Tracks(static_cast<std::string>(Acronym::PiMinus), fPacked_MC_PiMinus);
                ConnectBranches_MC_Tracks(static_cast<std::string>(Acronym::PiPlus), fPacked_MC_PiPlus);
            }
            break;
        case ReactionChannel::AntiH:
            ConnectBranches_Tracks(static_cast<std::string>(Acronym::NegKaon), fPacked_NegKaons);
            if (IsMC()) ConnectBranches_MC_Tracks(static_cast<std::string>(Acronym::NegKaon), fPacked_MC_NegKaons);
            break;
        default:
            break;
    }  // end of switch statement

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Finder::ConnectBranches_Events() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
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
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Finder::ConnectBranches_Injected() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
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
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Finder::ConnectBranches_V0s(const std::string& name_v0, PackedEvents::V0s& vec_v0s) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Entry", &vec_v0s.Entry, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_X", &vec_v0s.X, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Y", &vec_v0s.Y, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Z", &vec_v0s.Z, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Px", &vec_v0s.Px, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Py", &vec_v0s.Py, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Pz", &vec_v0s.Pz, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_E", &vec_v0s.E, name_v0);

    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaX2", &vec_v0s.Sigma.X2, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaXY", &vec_v0s.Sigma.XY, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaY2", &vec_v0s.Sigma.Y2, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaXZ", &vec_v0s.Sigma.XZ, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaYZ", &vec_v0s.Sigma.YZ, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaZ2", &vec_v0s.Sigma.Z2, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaXPx", &vec_v0s.Sigma.XPx, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaYPx", &vec_v0s.Sigma.YPx, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaZPx", &vec_v0s.Sigma.ZPx, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaPx2", &vec_v0s.Sigma.Px2, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaXPy", &vec_v0s.Sigma.XPy, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaYPy", &vec_v0s.Sigma.YPy, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaZPy", &vec_v0s.Sigma.ZPy, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaPxPy", &vec_v0s.Sigma.PxPy, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaPy2", &vec_v0s.Sigma.Py2, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaXPz", &vec_v0s.Sigma.XPz, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaYPz", &vec_v0s.Sigma.YPz, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaZPz", &vec_v0s.Sigma.ZPz, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaPxPz", &vec_v0s.Sigma.PxPz, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaPyPz", &vec_v0s.Sigma.PyPz, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaPz2", &vec_v0s.Sigma.Pz2, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaXE", &vec_v0s.Sigma.XE, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaYE", &vec_v0s.Sigma.YE, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaZE", &vec_v0s.Sigma.ZE, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaPxE", &vec_v0s.Sigma.PxE, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaPyE", &vec_v0s.Sigma.PyE, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaPzE", &vec_v0s.Sigma.PzE, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaE2", &vec_v0s.Sigma.E2, name_v0);

    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Neg_Entry", &vec_v0s.Neg.Entry, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Neg_X", &vec_v0s.Neg.X, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Neg_Y", &vec_v0s.Neg.Y, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Neg_Z", &vec_v0s.Neg.Z, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Neg_Px", &vec_v0s.Neg.Px, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Neg_Py", &vec_v0s.Neg.Py, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Neg_Pz", &vec_v0s.Neg.Pz, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Neg_E", &vec_v0s.Neg.E, name_v0);

    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Neg_X_AtPCA", &vec_v0s.Neg_X_AtPCA, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Neg_Y_AtPCA", &vec_v0s.Neg_Y_AtPCA, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Neg_Z_AtPCA", &vec_v0s.Neg_Z_AtPCA, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Neg_Px_AtPCA", &vec_v0s.Neg_Px_AtPCA, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Neg_Py_AtPCA", &vec_v0s.Neg_Py_AtPCA, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Neg_Pz_AtPCA", &vec_v0s.Neg_Pz_AtPCA, name_v0);

    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Pos_Entry", &vec_v0s.Pos.Entry, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Pos_X", &vec_v0s.Pos.X, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Pos_Y", &vec_v0s.Pos.Y, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Pos_Z", &vec_v0s.Pos.Z, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Pos_Px", &vec_v0s.Pos.Px, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Pos_Py", &vec_v0s.Pos.Py, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Pos_Pz", &vec_v0s.Pos.Pz, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Pos_E", &vec_v0s.Pos.E, name_v0);

    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Pos_X_AtPCA", &vec_v0s.Pos_X_AtPCA, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Pos_Y_AtPCA", &vec_v0s.Pos_Y_AtPCA, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Pos_Z_AtPCA", &vec_v0s.Pos_Z_AtPCA, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Pos_Px_AtPCA", &vec_v0s.Pos_Px_AtPCA, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Pos_Py_AtPCA", &vec_v0s.Pos_Py_AtPCA, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Pos_Pz_AtPCA", &vec_v0s.Pos_Pz_AtPCA, name_v0);

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Finder::ConnectBranches_Tracks(const std::string& name_part, PackedEvents::Tracks& vec_tracks) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Entry", &vec_tracks.Entry, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_X", &vec_tracks.X, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Y", &vec_tracks.Y, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Z", &vec_tracks.Z, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Px", &vec_tracks.Px, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Py", &vec_tracks.Py, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_Pz", &vec_tracks.Pz, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_E", &vec_tracks.E, name_part);

    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaX2", &vec_tracks.Sigma.X2, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaXY", &vec_tracks.Sigma.XY, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaY2", &vec_tracks.Sigma.Y2, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaXZ", &vec_tracks.Sigma.XZ, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaYZ", &vec_tracks.Sigma.YZ, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaZ2", &vec_tracks.Sigma.Z2, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaXPx", &vec_tracks.Sigma.XPx, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaYPx", &vec_tracks.Sigma.YPx, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaZPx", &vec_tracks.Sigma.ZPx, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaPx2", &vec_tracks.Sigma.Px2, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaXPy", &vec_tracks.Sigma.XPy, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaYPy", &vec_tracks.Sigma.YPy, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaZPy", &vec_tracks.Sigma.ZPy, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaPxPy", &vec_tracks.Sigma.PxPy, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaPy2", &vec_tracks.Sigma.Py2, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaXPz", &vec_tracks.Sigma.XPz, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaYPz", &vec_tracks.Sigma.YPz, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaZPz", &vec_tracks.Sigma.ZPz, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaPxPz", &vec_tracks.Sigma.PxPz, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaPyPz", &vec_tracks.Sigma.PyPz, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaPz2", &vec_tracks.Sigma.Pz2, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaXE", &vec_tracks.Sigma.XE, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaYE", &vec_tracks.Sigma.YE, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaZE", &vec_tracks.Sigma.ZE, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaPxE", &vec_tracks.Sigma.PxE, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaPyE", &vec_tracks.Sigma.PyE, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaPzE", &vec_tracks.Sigma.PzE, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_SigmaE2", &vec_tracks.Sigma.E2, name_part);

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Finder::ConnectBranches_MC_V0s(const std::string& name_v0, PackedEvents::MC_V0s& sov) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Entry", &sov.Entry, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_X", &sov.X, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Y", &sov.Y, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Z", &sov.Z, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Px", &sov.Px, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Py", &sov.Py, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Pz", &sov.Pz, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_E", &sov.E, name_v0);

    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_DecayX", &sov.DecayX, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_DecayY", &sov.DecayY, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_DecayZ", &sov.DecayZ, name_v0);

    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_PdgCode", &sov.PdgCode, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Mother_Entry", &sov.Mother_Entry, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Mother_PdgCode", &sov.Mother_PdgCode, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_IsTrue", &sov.IsTrue, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_IsSignal", &sov.IsSignal, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_IsSecondary", &sov.IsSecondary, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_ReactionID", &sov.ReactionID, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_IsHybrid", &sov.IsHybrid, name_v0);

    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Neg_Entry", &sov.Neg_Entry, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Neg_Px", &sov.Neg_Px, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Neg_Py", &sov.Neg_Py, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Neg_Pz", &sov.Neg_Pz, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Neg_PdgCode", &sov.Neg_PdgCode, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Neg_IsTrue", &sov.Neg_IsTrue, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Neg_IsSignal", &sov.Neg_IsSignal, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Neg_IsSecondary", &sov.Neg_IsSecondary, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Neg_ReactionID", &sov.Neg_ReactionID, name_v0);

    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Pos_Entry", &sov.Pos_Entry, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Pos_Px", &sov.Pos_Px, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Pos_Py", &sov.Pos_Py, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Pos_Pz", &sov.Pos_Pz, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Pos_PdgCode", &sov.Pos_PdgCode, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Pos_IsTrue", &sov.Pos_IsTrue, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Pos_IsSignal", &sov.Pos_IsSignal, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Pos_IsSecondary", &sov.Pos_IsSecondary, name_v0);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Pos_ReactionID", &sov.Pos_ReactionID, name_v0);

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Finder::ConnectBranches_MC_Tracks(const std::string& name_part, PackedEvents::MC_Tracks& sov) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Entry", &sov.Entry, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_X", &sov.X, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Y", &sov.Y, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Z", &sov.Z, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Px", &sov.Px, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Py", &sov.Py, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Pz", &sov.Pz, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_E", &sov.E, name_part);

    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_PdgCode", &sov.PdgCode, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Mother_Entry", &sov.Mother_Entry, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_Mother_PdgCode", &sov.Mother_PdgCode, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_GrandMother_Entry", &sov.GrandMother_Entry, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_GrandMother_PdgCode", &sov.GrandMother_PdgCode, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_IsTrue", &sov.IsTrue, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_IsSignal", &sov.IsSignal, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_IsSecondary", &sov.IsSecondary, name_part);
    Utils::ConnectBranch(fTree_PackedEvents.get(), "_MC_ReactionID", &sov.ReactionID, name_part);

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

// ## OUTPUT ZONE ## //

bool Finder::PrepareOutputFile() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    const std::filesystem::path output_path(fSettings.PathOutputFile);
    if (output_path.has_parent_path()) std::filesystem::create_directories(output_path.parent_path());

    fOutputFile = std::unique_ptr<TFile>(TFile::Open(fSettings.PathOutputFile.c_str(), "RECREATE"));
    if (!fOutputFile) {
        std::cerr << "   " << __FUNCTION__ << " :: TFile \"" << fSettings.PathOutputFile << "\" couldn't be created" << '\n';
        return false;
    }

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
    return true;
}

bool Finder::PrepareOutputTree() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    std::string tree_name{"Candidates"};
    switch (GetReactionChannel()) {
        case ReactionChannel::A:
            tree_name += static_cast<std::string>(Name::ChannelA);
            break;
        case ReactionChannel::D:
            tree_name += static_cast<std::string>(Name::ChannelD);
            break;
        case ReactionChannel::E:
            tree_name += static_cast<std::string>(Name::ChannelE);
            break;
        case ReactionChannel::H:
            tree_name += static_cast<std::string>(Name::ChannelH);
            break;
        case ReactionChannel::AntiA:
            tree_name += static_cast<std::string>(Name::AntiA);
            break;
        case ReactionChannel::AntiD:
            tree_name += static_cast<std::string>(Name::AntiD);
            break;
        case ReactionChannel::AntiE:
            tree_name += static_cast<std::string>(Name::AntiE);
            break;
        case ReactionChannel::AntiH:
            tree_name += static_cast<std::string>(Name::AntiH);
            break;
        default:
            return false;
    }

    fOutputTree = std::make_unique<TTree>(tree_name.c_str(), "Final results");
    if (!fOutputTree) {
        std::cerr << "   " << __FUNCTION__ << " :: TTree \"" << tree_name << "\" couldn't be created" << '\n';
        return false;
    }

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
    return true;
}

void Finder::CreateOutputBranches(Found::ChannelA& out_branches) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    fOutputTree->Branch("RunNumber", &out_branches.RunNumber);
    fOutputTree->Branch("DirNumber", &out_branches.DirNumber);
    if (!IsMC()) fOutputTree->Branch("DirNumberB", &out_branches.DirNumberB);
    fOutputTree->Branch("EventNumber", &out_branches.EventNumber);
    fOutputTree->Branch("MagneticField", &out_branches.MagneticField);
    fOutputTree->Branch("PV_Xv", &out_branches.PV_Xv);
    fOutputTree->Branch("PV_Yv", &out_branches.PV_Yv);
    fOutputTree->Branch("PV_Zv", &out_branches.PV_Zv);
    if (IsMC()) {
        fOutputTree->Branch("MC_PV_Xv", &out_branches.MC_PV_Xv);
        fOutputTree->Branch("MC_PV_Yv", &out_branches.MC_PV_Yv);
        fOutputTree->Branch("MC_PV_Zv", &out_branches.MC_PV_Zv);
    }

    fOutputTree->Branch("X", &out_branches.X);
    fOutputTree->Branch("Y", &out_branches.Y);
    fOutputTree->Branch("Z", &out_branches.Z);
    fOutputTree->Branch("Px", &out_branches.Px);
    fOutputTree->Branch("Py", &out_branches.Py);
    fOutputTree->Branch("Pz", &out_branches.Pz);
    fOutputTree->Branch("E", &out_branches.E);
    fOutputTree->Branch("E_MinusNucleon", &out_branches.E_MinusNucleon);

    fOutputTree->Branch("V0A_X_AtPCA", &out_branches.V0A_X_AtPCA);
    fOutputTree->Branch("V0A_Y_AtPCA", &out_branches.V0A_Y_AtPCA);
    fOutputTree->Branch("V0A_Z_AtPCA", &out_branches.V0A_Z_AtPCA);

    fOutputTree->Branch("V0B_X_AtPCA", &out_branches.V0B_X_AtPCA);
    fOutputTree->Branch("V0B_Y_AtPCA", &out_branches.V0B_Y_AtPCA);
    fOutputTree->Branch("V0B_Z_AtPCA", &out_branches.V0B_Z_AtPCA);

    fOutputTree->Branch("V0A_Entry", &out_branches.V0A_Entry);
    fOutputTree->Branch("V0A_X_AtDecay", &out_branches.V0A.X);
    fOutputTree->Branch("V0A_Y_AtDecay", &out_branches.V0A.Y);
    fOutputTree->Branch("V0A_Z_AtDecay", &out_branches.V0A.Z);
    fOutputTree->Branch("V0A_Px", &out_branches.V0A.Px);
    fOutputTree->Branch("V0A_Py", &out_branches.V0A.Py);
    fOutputTree->Branch("V0A_Pz", &out_branches.V0A.Pz);
    fOutputTree->Branch("V0A_E", &out_branches.V0A.E);

    fOutputTree->Branch("V0A_Neg_Entry", &out_branches.V0A_Neg_Entry);
    fOutputTree->Branch("V0A_Pos_Entry", &out_branches.V0A_Pos_Entry);

    fOutputTree->Branch("V0B_Entry", &out_branches.V0B_Entry);
    fOutputTree->Branch("V0B_X_AtDecay", &out_branches.V0B.X);
    fOutputTree->Branch("V0B_Y_AtDecay", &out_branches.V0B.Y);
    fOutputTree->Branch("V0B_Z_AtDecay", &out_branches.V0B.Z);
    fOutputTree->Branch("V0B_Px", &out_branches.V0B.Px);
    fOutputTree->Branch("V0B_Py", &out_branches.V0B.Py);
    fOutputTree->Branch("V0B_Pz", &out_branches.V0B.Pz);
    fOutputTree->Branch("V0B_E", &out_branches.V0B.E);

    fOutputTree->Branch("V0B_Neg_Entry", &out_branches.V0B_Neg_Entry);
    fOutputTree->Branch("V0B_Pos_Entry", &out_branches.V0B_Pos_Entry);

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Finder::CreateOutputBranches(Found::ChannelD& out_branches) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    fOutputTree->Branch("RunNumber", &out_branches.RunNumber);
    fOutputTree->Branch("DirNumber", &out_branches.DirNumber);
    if (!IsMC()) fOutputTree->Branch("DirNumberB", &out_branches.DirNumberB);
    fOutputTree->Branch("EventNumber", &out_branches.EventNumber);
    fOutputTree->Branch("MagneticField", &out_branches.MagneticField);
    fOutputTree->Branch("PV_Xv", &out_branches.PV_Xv);
    fOutputTree->Branch("PV_Yv", &out_branches.PV_Yv);
    fOutputTree->Branch("PV_Zv", &out_branches.PV_Zv);
    if (IsMC()) {
        fOutputTree->Branch("MC_PV_Xv", &out_branches.MC_PV_Xv);
        fOutputTree->Branch("MC_PV_Yv", &out_branches.MC_PV_Yv);
        fOutputTree->Branch("MC_PV_Zv", &out_branches.MC_PV_Zv);
    }

    fOutputTree->Branch("X", &out_branches.X);
    fOutputTree->Branch("Y", &out_branches.Y);
    fOutputTree->Branch("Z", &out_branches.Z);
    fOutputTree->Branch("Px", &out_branches.Px);
    fOutputTree->Branch("Py", &out_branches.Py);
    fOutputTree->Branch("Pz", &out_branches.Pz);
    fOutputTree->Branch("E", &out_branches.E);
    fOutputTree->Branch("E_MinusNucleon", &out_branches.E_MinusNucleon);

    fOutputTree->Branch("V0_X_AtPCA", &out_branches.V0_X_AtPCA);
    fOutputTree->Branch("V0_Y_AtPCA", &out_branches.V0_Y_AtPCA);
    fOutputTree->Branch("V0_Z_AtPCA", &out_branches.V0_Z_AtPCA);

    fOutputTree->Branch("Kaon_X_AtPCA", &out_branches.Kaon_X_AtPCA);
    fOutputTree->Branch("Kaon_Y_AtPCA", &out_branches.Kaon_Y_AtPCA);
    fOutputTree->Branch("Kaon_Z_AtPCA", &out_branches.Kaon_Z_AtPCA);
    fOutputTree->Branch("Kaon_Px_AtPCA", &out_branches.Kaon_Px_AtPCA);
    fOutputTree->Branch("Kaon_Py_AtPCA", &out_branches.Kaon_Py_AtPCA);
    fOutputTree->Branch("Kaon_Pz_AtPCA", &out_branches.Kaon_Pz_AtPCA);

    fOutputTree->Branch("V0_Entry", &out_branches.V0_Entry);
    fOutputTree->Branch("V0_X_AtDecay", &out_branches.V0.X);
    fOutputTree->Branch("V0_Y_AtDecay", &out_branches.V0.Y);
    fOutputTree->Branch("V0_Z_AtDecay", &out_branches.V0.Z);
    fOutputTree->Branch("V0_Px", &out_branches.V0.Px);
    fOutputTree->Branch("V0_Py", &out_branches.V0.Py);
    fOutputTree->Branch("V0_Pz", &out_branches.V0.Pz);
    fOutputTree->Branch("V0_E", &out_branches.V0.E);

    fOutputTree->Branch("V0_Neg_Entry", &out_branches.V0_Neg_Entry);
    fOutputTree->Branch("V0_Pos_Entry", &out_branches.V0_Pos_Entry);

    fOutputTree->Branch("Kaon_Entry", &out_branches.Kaon_Entry);
    fOutputTree->Branch("Kaon_X", &out_branches.Kaon.X);
    fOutputTree->Branch("Kaon_Y", &out_branches.Kaon.Y);
    fOutputTree->Branch("Kaon_Z", &out_branches.Kaon.Z);
    fOutputTree->Branch("Kaon_Px", &out_branches.Kaon.Px);
    fOutputTree->Branch("Kaon_Py", &out_branches.Kaon.Py);
    fOutputTree->Branch("Kaon_Pz", &out_branches.Kaon.Pz);
    fOutputTree->Branch("Kaon_E", &out_branches.Kaon.E);

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
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
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    fOutputTree->Branch("MC_X", &sov.X);
    fOutputTree->Branch("MC_Y", &sov.Y);
    fOutputTree->Branch("MC_Z", &sov.Z);
    fOutputTree->Branch("MC_Px_Before", &sov.BeforePx);
    fOutputTree->Branch("MC_Py_Before", &sov.BeforePy);
    fOutputTree->Branch("MC_Pz_Before", &sov.BeforePz);
    fOutputTree->Branch("MC_E_Before", &sov.BeforeE);

    fOutputTree->Branch("MC_Px_After", &sov.AfterPx);
    fOutputTree->Branch("MC_Py_After", &sov.AfterPy);
    fOutputTree->Branch("MC_Pz_After", &sov.AfterPz);
    fOutputTree->Branch("MC_E_After", &sov.AfterE);

    fOutputTree->Branch("MC_IsSignal", &sov.IsSignal);
    fOutputTree->Branch("MC_ReactionID", &sov.ReactionID);
    fOutputTree->Branch("MC_IsHybrid", &sov.IsHybrid);

    fOutputTree->Branch("MC_Px_Nucleon", &sov.NucleonPx);
    fOutputTree->Branch("MC_Py_Nucleon", &sov.NucleonPy);
    fOutputTree->Branch("MC_Pz_Nucleon", &sov.NucleonPz);
    fOutputTree->Branch("MC_E_Nucleon", &sov.NucleonPz);

    fOutputTree->Branch("MC_V0A_Entry", &sov.V0A_Entry);
    fOutputTree->Branch("MC_V0A_PdgCode", &sov.V0A_PdgCode);
    fOutputTree->Branch("MC_V0A_Mother_Entry", &sov.V0A_Mother_Entry);
    fOutputTree->Branch("MC_V0A_Mother_PdgCode", &sov.V0A_Mother_PdgCode);
    fOutputTree->Branch("MC_V0A_Neg_Entry", &sov.V0A_Neg_Entry);
    fOutputTree->Branch("MC_V0A_Pos_Entry", &sov.V0A_Pos_Entry);
    fOutputTree->Branch("MC_V0A_Px", &sov.V0A_Px);
    fOutputTree->Branch("MC_V0A_Py", &sov.V0A_Py);
    fOutputTree->Branch("MC_V0A_Pz", &sov.V0A_Pz);
    fOutputTree->Branch("MC_V0A_E", &sov.V0A_E);
    fOutputTree->Branch("MC_V0A_IsTrue", &sov.V0A_IsTrue);
    fOutputTree->Branch("MC_V0A_IsSignal", &sov.V0A_IsSignal);
    fOutputTree->Branch("MC_V0A_IsSecondary", &sov.V0A_IsSecondary);
    fOutputTree->Branch("MC_V0A_ReactionID", &sov.V0A_ReactionID);
    fOutputTree->Branch("MC_V0A_IsHybrid", &sov.V0A_IsHybrid);

    fOutputTree->Branch("MC_V0B_Entry", &sov.V0B_Entry);
    fOutputTree->Branch("MC_V0B_PdgCode", &sov.V0B_PdgCode);
    fOutputTree->Branch("MC_V0B_Mother_Entry", &sov.V0B_Mother_Entry);
    fOutputTree->Branch("MC_V0B_Mother_PdgCode", &sov.V0B_Mother_PdgCode);
    fOutputTree->Branch("MC_V0B_Neg_Entry", &sov.V0B_Neg_Entry);
    fOutputTree->Branch("MC_V0B_Pos_Entry", &sov.V0B_Pos_Entry);
    fOutputTree->Branch("MC_V0B_Px", &sov.V0B_Px);
    fOutputTree->Branch("MC_V0B_Py", &sov.V0B_Py);
    fOutputTree->Branch("MC_V0B_Pz", &sov.V0B_Pz);
    fOutputTree->Branch("MC_V0B_E", &sov.V0B_E);
    fOutputTree->Branch("MC_V0B_IsTrue", &sov.V0B_IsTrue);
    fOutputTree->Branch("MC_V0B_IsSignal", &sov.V0B_IsSignal);
    fOutputTree->Branch("MC_V0B_IsSecondary", &sov.V0B_IsSecondary);
    fOutputTree->Branch("MC_V0B_ReactionID", &sov.V0B_ReactionID);
    fOutputTree->Branch("MC_V0B_IsHybrid", &sov.V0B_IsHybrid);

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Finder::CreateOutputBranches(Found::MC_ChannelD& sov) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    fOutputTree->Branch("MC_X", &sov.X);
    fOutputTree->Branch("MC_Y", &sov.Y);
    fOutputTree->Branch("MC_Z", &sov.Z);
    fOutputTree->Branch("MC_Px_Before", &sov.BeforePx);
    fOutputTree->Branch("MC_Py_Before", &sov.BeforePy);
    fOutputTree->Branch("MC_Pz_Before", &sov.BeforePz);
    fOutputTree->Branch("MC_E_Before", &sov.BeforeE);

    fOutputTree->Branch("MC_Px_After", &sov.AfterPx);
    fOutputTree->Branch("MC_Py_After", &sov.AfterPy);
    fOutputTree->Branch("MC_Pz_After", &sov.AfterPz);
    fOutputTree->Branch("MC_E_After", &sov.AfterE);

    fOutputTree->Branch("MC_IsSignal", &sov.IsSignal);
    fOutputTree->Branch("MC_ReactionID", &sov.ReactionID);
    fOutputTree->Branch("MC_IsHybrid", &sov.IsHybrid);

    fOutputTree->Branch("MC_Nucleon_Px", &sov.NucleonPx);
    fOutputTree->Branch("MC_Nucleon_Py", &sov.NucleonPy);
    fOutputTree->Branch("MC_Nucleon_Pz", &sov.NucleonPz);
    fOutputTree->Branch("MC_Nucleon_E", &sov.NucleonPz);

    fOutputTree->Branch("MC_V0_Entry", &sov.V0_Entry);
    fOutputTree->Branch("MC_V0_PdgCode", &sov.V0_PdgCode);
    fOutputTree->Branch("MC_V0_Mother_Entry", &sov.V0_Mother_Entry);
    fOutputTree->Branch("MC_V0_Mother_PdgCode", &sov.V0_Mother_PdgCode);
    fOutputTree->Branch("MC_V0_Neg_Entry", &sov.V0_Neg_Entry);
    fOutputTree->Branch("MC_V0_Pos_Entry", &sov.V0_Pos_Entry);
    fOutputTree->Branch("MC_V0_Px", &sov.V0_Px);
    fOutputTree->Branch("MC_V0_Py", &sov.V0_Py);
    fOutputTree->Branch("MC_V0_Pz", &sov.V0_Pz);
    fOutputTree->Branch("MC_V0_E", &sov.V0_E);
    fOutputTree->Branch("MC_V0_IsTrue", &sov.V0_IsTrue);
    fOutputTree->Branch("MC_V0_IsSignal", &sov.V0_IsSignal);
    fOutputTree->Branch("MC_V0_IsSecondary", &sov.V0_IsSecondary);
    fOutputTree->Branch("MC_V0_ReactionID", &sov.V0_ReactionID);
    fOutputTree->Branch("MC_V0_IsHybrid", &sov.V0_IsHybrid);

    fOutputTree->Branch("MC_Kaon_Entry", &sov.Kaon_Entry);
    fOutputTree->Branch("MC_Kaon_PdgCode", &sov.Kaon_PdgCode);
    fOutputTree->Branch("MC_Kaon_Mother_Entry", &sov.Kaon_Mother_Entry);
    fOutputTree->Branch("MC_Kaon_Mother_PdgCode", &sov.Kaon_Mother_PdgCode);
    fOutputTree->Branch("MC_Kaon_GrandMother_Entry", &sov.Kaon_GrandMother_Entry);
    fOutputTree->Branch("MC_Kaon_GrandMother_PdgCode", &sov.Kaon_GrandMother_PdgCode);
    fOutputTree->Branch("MC_Kaon_Px", &sov.Kaon_Px);
    fOutputTree->Branch("MC_Kaon_Py", &sov.Kaon_Py);
    fOutputTree->Branch("MC_Kaon_Pz", &sov.Kaon_Pz);
    fOutputTree->Branch("MC_Kaon_E", &sov.Kaon_E);
    fOutputTree->Branch("MC_Kaon_IsTrue", &sov.Kaon_IsTrue);
    fOutputTree->Branch("MC_Kaon_IsSignal", &sov.Kaon_IsSignal);
    fOutputTree->Branch("MC_Kaon_IsSecondary", &sov.Kaon_IsSecondary);
    fOutputTree->Branch("MC_Kaon_ReactionID", &sov.Kaon_ReactionID);

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

// ## Injected ZONE ## //

bool Finder::Injected_PrepareOutputTree() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    std::string tree_name{"Injected"};

    fOutputTree_Injected = std::make_unique<TTree>(tree_name.c_str(), "");
    if (!fOutputTree_Injected) {
        std::cerr << "   " << __FUNCTION__ << " :: TTree \"" << tree_name << "\" couldn't be created" << '\n';
        return false;
    }

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
    return true;
}

void Finder::Injected_CreateOutputBranches() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
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
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
};

void Finder::Injected_FlattenAndStore() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
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
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

// ## Channel A ZONE ## //

void Finder::FindSexaquarks_ChannelA(bool anti_channel) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    const PackedEvents::V0s* Packed_Lambdas{anti_channel ? &fPacked_Lambdas : &fPacked_AntiLambdas};
    const PackedEvents::MC_V0s* MC_Lambdas{anti_channel ? &fPacked_MC_Lambdas : &fPacked_MC_AntiLambdas};
    PdgCode PdgCode_Lambda{anti_channel ? PdgCode::Lambda : PdgCode::AntiLambda};

    // loop over all possible pairs of (anti)lambda + K0S //
    auto n_lambdas = static_cast<int>(Packed_Lambdas->Entry->size());
    auto n_k0s = static_cast<int>(fPacked_KaonsZeroShort.Entry->size());
    for (int idx_v0a{0}; idx_v0a < n_lambdas; ++idx_v0a) {

        // unpack (anti)lambda //
        auto v0a = KF::UnpackV0(*Packed_Lambdas, idx_v0a, PdgCode_Lambda);

        for (int idx_v0b{0}; idx_v0b < n_k0s; ++idx_v0b) {

            // unpack K0S //
            auto v0b = KF::UnpackV0(fPacked_KaonsZeroShort, idx_v0b, PdgCode::KaonZeroShort);

            // sanity check //
            std::set<int> unique_track_entries{v0a.Neg.idx, v0a.Pos.idx, v0b.Neg.idx, v0b.Pos.idx};
            if (unique_track_entries.size() < 4) continue;

            // fit //
            KF::ChannelA sexa{PdgMass::Neutron, v0a, v0b, fInput_Event.MagneticField};

#ifdef T2S_DEBUG
            std::cout << "   " << __FUNCTION__ << " :: idx(v0a,neg,pos)=" << v0a.idx << "," << v0a.Neg.idx << "," << v0a.Pos.idx;
            std::cout << ";x,y,z=" << v0a.X() << "," << v0a.Y() << "," << v0a.Z();
            std::cout << ";px,py,pz=" << v0a.Px() << "," << v0a.Py() << "," << v0a.Pz();
            std::cout << ";mass=" << v0a.Mass() << '\n';

            std::cout << "   " << __FUNCTION__ << " :: idx(v0b,neg,pos)=" << v0b.idx << "," << v0b.Neg.idx << "," << v0b.Pos.idx;
            std::cout << ";x,y,z=" << v0b.X() << "," << v0b.Y() << "," << v0b.Z();
            std::cout << ";px,py,pz=" << v0b.Px() << "," << v0b.Py() << "," << v0b.Pz();
            std::cout << ";mass=" << v0b.Mass() << '\n';

            std::cout << "   " << __FUNCTION__ << " :: x,y,z=" << sexa.X() << "," << sexa.Y() << "," << sexa.Z();
            std::cout << ";x,y,z(v0a)=" << sexa.V0A_PCA_XYZ()[0] << "," << sexa.V0A_PCA_XYZ()[1] << "," << sexa.V0A_PCA_XYZ()[2];
            std::cout << ";x,y,z(v0b)=" << sexa.V0B_PCA_XYZ()[0] << "," << sexa.V0B_PCA_XYZ()[1] << "," << sexa.V0B_PCA_XYZ()[2];
            std::cout << ";mass=" << sexa.Mass();
            std::cout << ";mass_minus_n=" << sexa.Mass_MinusNucleon();
            std::cout << ";dca_btw_v0s=" << sexa.DCA_btw_V0s();
            std::cout << ";radius=" << sexa.Radius();
            std::cout << ";dca_v0a=" << sexa.DCA_V0A_wrt_SV();
            std::cout << ";dca_v0b=" << sexa.DCA_V0B_wrt_SV();
            std::cout << ";dca_v0a_neg=" << sexa.DCA_V0ANeg_wrt_SV(fInput_Event.MagneticField);
            std::cout << ";dca_v0a_pos=" << sexa.DCA_V0APos_wrt_SV(fInput_Event.MagneticField);
            std::cout << ";dca_v0b_neg=" << sexa.DCA_V0BNeg_wrt_SV(fInput_Event.MagneticField);
            std::cout << ";dca_v0b_pos=" << sexa.DCA_V0BPos_wrt_SV(fInput_Event.MagneticField);
            std::cout << ";pt=" << sexa.Pt();
            std::cout << ";eta=" << sexa.Eta();
            std::cout << ";decay_length(v0a,v0b)=" << sexa.DecayLength_V0A() << "," << sexa.DecayLength_V0B();
            std::cout << ";cpa_pv=" << sexa.CPA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv) << '\n';

            MC::V0 mc_v0a_debug{*MC_Lambdas, v0a.idx};
            MC::V0 mc_v0b_debug{fPacked_MC_KaonsZeroShort, v0b.idx};
            MC::ChannelA mc_sexa_debug{fInput_Injected, fSettings.SexaquarkMass, mc_v0a_debug, mc_v0b_debug};
            if (IsMC()) {
                std::cout << "   " << __FUNCTION__ << " :: ";
                std::cout << "sexa(is_signal,reaction_id,is_hybrid)=" << mc_sexa_debug.IsSignal << "," << mc_sexa_debug.ReactionID << ","
                          << mc_sexa_debug.IsHybrid;
                std::cout << ";v0a(entry,pdg_code)=" << mc_sexa_debug.V0A.Entry << ", " << mc_sexa_debug.V0A.PdgCode;
                std::cout << ";v0a(is_signal,reaction_id,is_hybrid)=" << mc_sexa_debug.V0A.IsSignal << "," << mc_sexa_debug.V0A.ReactionID << ","
                          << mc_sexa_debug.V0A.IsHybrid;
                std::cout << ";v0b(entry,pdg_code)=" << mc_sexa_debug.V0B.Entry << ", " << mc_sexa_debug.V0B.PdgCode;
                std::cout << ";v0b(is_signal,reaction_id,is_hybrid)=" << mc_sexa_debug.V0B.IsSignal << "," << mc_sexa_debug.V0B.ReactionID << ","
                          << mc_sexa_debug.V0B.IsHybrid << '\n';
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
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

bool Finder::PassesCuts(const KF::ChannelA& sexa) const {

    if (sexa.Radius2D() < Cuts::ChannelA::Min_Radius2D || sexa.Radius2D() > Cuts::ChannelA::Max_Radius2D) return false;
    if (sexa.DecayLength_V0A() > Cuts::ChannelA::Max_DecayLengthLa) return false;
    if (sexa.DecayLength_V0B() > Cuts::ChannelA::Max_DecayLengthK0) return false;
    if (sexa.AbsRapidity_MinusNucleon() > Cuts::ChannelA::AbsMax_Rapidity) return false;  // PENDING: kinematic cut, affected by Fermi motion
    if (sexa.Mass_MinusNucleon() < Cuts::ChannelA::Min_MassMinusNucleon || sexa.Mass_MinusNucleon() > Cuts::ChannelA::Max_MassMinusNucleon)
        return false;  // PENDING: kinematic cut, affected by Fermi motion
    if (sexa.CPA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv) < Cuts::ChannelA::Min_CPAwrtPV ||
        sexa.CPA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv) > Cuts::ChannelA::Max_CPAwrtPV)
        return false;  // PENDING: kinematic cut, affected by Fermi motion
    if (sexa.DCA_V0ANeg_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelA::Max_DCALaNegSV) return false;
    if (sexa.DCA_V0APos_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelA::Max_DCALaPosSV) return false;
    if (sexa.DCA_V0BNeg_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelA::Max_DCAK0NegSV) return false;
    if (sexa.DCA_V0BPos_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelA::Max_DCAK0PosSV) return false;
    if (sexa.DCA_V0A_wrt_SV() > Cuts::ChannelA::Max_DCALaSV) return false;
    if (sexa.DCA_V0B_wrt_SV() > Cuts::ChannelA::Max_DCAK0SV) return false;
    if (sexa.DCA_btw_V0s() > Cuts::ChannelA::Max_DCAbtwV0s) return false;

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
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    const PackedEvents::V0s* Packed_Lambdas{anti_channel ? &fPacked_Lambdas : &fPacked_AntiLambdas};
    const PackedEvents::Tracks* Packed_Kaons{anti_channel ? &fPacked_NegKaons : &fPacked_PosKaons};
    const PackedEvents::MC_V0s* MC_Lambdas{anti_channel ? &fPacked_MC_Lambdas : &fPacked_MC_AntiLambdas};
    const PackedEvents::MC_Tracks* MC_Kaons{anti_channel ? &fPacked_MC_NegKaons : &fPacked_MC_PosKaons};
    PdgCode PdgCode_Lambda{anti_channel ? PdgCode::Lambda : PdgCode::AntiLambda};
    int Charge_Kaon{anti_channel ? -1 : +1};

    // loop over all possible pairs of (anti)lambda + (pos/neg)kaon //
    auto n_lambdas = static_cast<int>(Packed_Lambdas->Entry->size());
    auto n_kaons = static_cast<int>(Packed_Kaons->Entry->size());
    for (int idx_v0{0}; idx_v0 < n_lambdas; ++idx_v0) {

        // unpack (anti)lambda //
        auto v0 = KF::UnpackV0(*Packed_Lambdas, idx_v0, PdgCode_Lambda);

        for (int idx_kaon{0}; idx_kaon < n_kaons; ++idx_kaon) {

            // unpack kaon //
            auto kaon = KF::UnpackTrack(*Packed_Kaons, idx_kaon, Charge_Kaon);

            // sanity check //
            std::set<int> unique_track_entries{v0.Neg.idx, v0.Pos.idx, kaon.idx};
            if (unique_track_entries.size() < 3) continue;

            // fit //
            KF::ChannelD sexa{PdgMass::Proton, v0, kaon, fInput_Event.MagneticField};

            // apply cuts //
            if (!PassesCuts(sexa)) continue;
#ifdef T2S_DEBUG
            std::cout << "   " << __FUNCTION__ << " :: idx(sexa,v0,neg,pos,kaon)="  //
                      << sexa.V0.idx << "," << sexa.V0.Neg.idx << "," << sexa.V0.Pos.idx << "," << sexa.Kaon.idx;
            std::cout << ";x,y,z=" << sexa.X() << "," << sexa.Y() << "," << sexa.Z();
            std::cout << ";x,y,z(v0)=" << sexa.V0_PCA_XYZ()[0] << "," << sexa.V0_PCA_XYZ()[1] << "," << sexa.V0_PCA_XYZ()[2];
            std::cout << ";x,y,z(kaon)=" << sexa.Kaon_PCA_XYZ()[0] << "," << sexa.Kaon_PCA_XYZ()[1] << "," << sexa.Kaon_PCA_XYZ()[2];
            std::cout << ";mass=" << sexa.Mass();
            std::cout << ";dca_v0_kaon=" << sexa.DCA_btw_V0_Kaon();
            std::cout << ";radius=" << sexa.Radius();
            std::cout << ";dca_v0=" << sexa.DCA_V0_wrt_SV();
            std::cout << ";dca_kaon=" << sexa.DCA_Kaon_wrt_SV();
            std::cout << ";pt=" << sexa.Pt();
            std::cout << ";eta=" << sexa.Eta();
            std::cout << ";cpa_pv=" << sexa.CPA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv) << '\n';
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
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

bool Finder::PassesCuts(const KF::ChannelD& sexa) const {

    if (sexa.Radius2D() < Cuts::ChannelD::Min_Radius2D || sexa.Radius2D() > Cuts::ChannelD::Max_Radius2D) return false;
    if (sexa.AbsRapidity_MinusNucleon() > Cuts::ChannelD::AbsMax_Rapidity) return false;  // PENDING: kinematics, affected by Fermi motion
    if (sexa.CPA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv) < Cuts::ChannelD::Min_CPAwrtPV ||
        sexa.CPA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv) > Cuts::ChannelD::Max_CPAwrtPV)
        return false;  // PENDING: kinematics, affected by Fermi motion
    if (sexa.DCA_V0_wrt_SV() > Cuts::ChannelD::Max_DCALaSV) return false;
    if (sexa.DCA_Kaon_wrt_SV() > Cuts::ChannelD::Max_DCAKaSV) return false;
    if (sexa.DCA_V0Neg_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelD::Max_DCALaNegSV) return false;
    if (sexa.DCA_V0Pos_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelD::Max_DCALaPosSV) return false;
    if (sexa.DCA_btw_V0_Kaon() > Cuts::ChannelD::Max_DCAKaLa) return false;
    if (sexa.DCA_btw_V0Neg_Kaon(fInput_Event.MagneticField) > Cuts::ChannelD::Max_DCALaNegKa) return false;
    if (sexa.DCA_btw_V0Pos_Kaon(fInput_Event.MagneticField) > Cuts::ChannelD::Max_DCALaPosKa) return false;

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
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    if (IsSignalMC()) {
        fOutputTree_Injected->Write();
        std::cout << "TTree \"" << fOutputTree_Injected->GetName() << "\" has been written onto TFile \"" << fSettings.PathOutputFile << "\"" << '\n';
    }

    fOutputTree->Write();
    std::cout << "TTree \"" << fOutputTree->GetName() << "\" has been written onto TFile \"" << fSettings.PathOutputFile << "\"" << '\n';

    fTree_PackedEvents->ResetBranchAddresses();
    fOutputTree->ResetBranchAddresses();
    std::cout << "Done." << '\n';

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

}  // namespace Tree2Secondaries
