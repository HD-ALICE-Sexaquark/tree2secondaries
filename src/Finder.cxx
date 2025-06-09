#include <iostream>
#include <memory>
#include <set>

#include "Finder/Cuts.hxx"
#include "Finder/Finder.hxx"
#include "Math/Constants.hxx"
#include "Math/KFWrapper.hxx"
#include "Structures/PackedEvents.hxx"

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

    switch (GetReactionChannel()) {
        // standard channels //
        case ReactionChannel::A:
            ConnectBranches_V0s(static_cast<std::string>(Acronym::AntiLambda), fPacked_AntiLambdas);
            ConnectBranches_V0s(static_cast<std::string>(Acronym::KaonZeroShort), fPacked_KaonsZeroShort);
            break;
        case ReactionChannel::D:
            ConnectBranches_V0s(static_cast<std::string>(Acronym::AntiLambda), fPacked_AntiLambdas);
            ConnectBranches_Tracks(static_cast<std::string>(Acronym::PosKaon), fPacked_PosKaons);
            break;
        case ReactionChannel::E:
            ConnectBranches_V0s(static_cast<std::string>(Acronym::AntiLambda), fPacked_AntiLambdas);
            ConnectBranches_Tracks(static_cast<std::string>(Acronym::PosKaon), fPacked_PosKaons);
            ConnectBranches_Tracks(static_cast<std::string>(Acronym::PiMinus), fPacked_PiMinus);
            ConnectBranches_Tracks(static_cast<std::string>(Acronym::PiPlus), fPacked_PiPlus);
            break;
        case ReactionChannel::H:
            ConnectBranches_Tracks(static_cast<std::string>(Acronym::PosKaon), fPacked_PosKaons);
            break;
        // anti-channels //
        case ReactionChannel::AntiA:
            ConnectBranches_V0s(static_cast<std::string>(Acronym::Lambda), fPacked_Lambdas);
            ConnectBranches_V0s(static_cast<std::string>(Acronym::KaonZeroShort), fPacked_KaonsZeroShort);
            break;
        case ReactionChannel::AntiD:
            ConnectBranches_V0s(static_cast<std::string>(Acronym::Lambda), fPacked_Lambdas);
            ConnectBranches_Tracks(static_cast<std::string>(Acronym::NegKaon), fPacked_NegKaons);
            break;
        case ReactionChannel::AntiE:
            ConnectBranches_V0s(static_cast<std::string>(Acronym::Lambda), fPacked_Lambdas);
            ConnectBranches_Tracks(static_cast<std::string>(Acronym::NegKaon), fPacked_NegKaons);
            ConnectBranches_Tracks(static_cast<std::string>(Acronym::PiMinus), fPacked_PiMinus);
            ConnectBranches_Tracks(static_cast<std::string>(Acronym::PiPlus), fPacked_PiPlus);
            break;
        case ReactionChannel::AntiH:
            ConnectBranches_Tracks(static_cast<std::string>(Acronym::NegKaon), fPacked_NegKaons);
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

    fTree_PackedEvents->SetBranchStatus("RunNumber", true);
    fTree_PackedEvents->SetBranchStatus("DirNumber", true);
    fTree_PackedEvents->SetBranchStatus("EventNumber", true);
    fTree_PackedEvents->SetBranchStatus("Centrality", true);
    fTree_PackedEvents->SetBranchStatus("MagneticField", true);
    fTree_PackedEvents->SetBranchStatus("PV_Xv", true);
    fTree_PackedEvents->SetBranchStatus("PV_Yv", true);
    fTree_PackedEvents->SetBranchStatus("PV_Zv", true);

    fTree_PackedEvents->SetBranchAddress("RunNumber", &fInput_Event.RunNumber);
    fTree_PackedEvents->SetBranchAddress("DirNumber", &fInput_Event.DirNumber);
    fTree_PackedEvents->SetBranchAddress("EventNumber", &fInput_Event.EventNumber);
    fTree_PackedEvents->SetBranchAddress("Centrality", &fInput_Event.Centrality);
    fTree_PackedEvents->SetBranchAddress("MagneticField", &fInput_Event.MagneticField);
    fTree_PackedEvents->SetBranchAddress("PV_Xv", &fInput_Event.PV_Xv);
    fTree_PackedEvents->SetBranchAddress("PV_Yv", &fInput_Event.PV_Yv);
    fTree_PackedEvents->SetBranchAddress("PV_Zv", &fInput_Event.PV_Zv);

    if (IsMC()) {
        fTree_PackedEvents->SetBranchStatus("MC_PV_Xv", true);
        fTree_PackedEvents->SetBranchStatus("MC_PV_Yv", true);
        fTree_PackedEvents->SetBranchStatus("MC_PV_Zv", true);

        fTree_PackedEvents->SetBranchAddress("MC_PV_Xv", &fInput_Event.MC_PV_Xv);
        fTree_PackedEvents->SetBranchAddress("MC_PV_Yv", &fInput_Event.MC_PV_Yv);
        fTree_PackedEvents->SetBranchAddress("MC_PV_Zv", &fInput_Event.MC_PV_Zv);
    }
#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Finder::ConnectBranches_V0s(const std::string& name_v0, PackedEvents::V0s& vec_v0s) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    fTree_PackedEvents->SetBranchStatus((name_v0 + "_Entry").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_X").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_Y").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_Z").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_Px").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_Py").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_Pz").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_E").c_str(), true);

    fTree_PackedEvents->SetBranchStatus((name_v0 + "_SigmaX2").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_SigmaXY").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_SigmaY2").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_SigmaXZ").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_SigmaYZ").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_SigmaZ2").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_SigmaXPx").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_SigmaYPx").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_SigmaZPx").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_SigmaPx2").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_SigmaXPy").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_SigmaYPy").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_SigmaZPy").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_SigmaPxPy").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_SigmaPy2").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_SigmaXPz").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_SigmaYPz").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_SigmaZPz").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_SigmaPxPz").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_SigmaPyPz").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_SigmaPz2").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_SigmaXE").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_SigmaYE").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_SigmaZE").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_SigmaPxE").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_SigmaPyE").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_SigmaPzE").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_SigmaE2").c_str(), true);

    fTree_PackedEvents->SetBranchStatus((name_v0 + "_Neg_Entry").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_Neg_X").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_Neg_Y").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_Neg_Z").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_Neg_Px").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_Neg_Py").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_Neg_Pz").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_Neg_E").c_str(), true);

    fTree_PackedEvents->SetBranchStatus((name_v0 + "_Pos_Entry").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_Pos_X").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_Pos_Y").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_Pos_Z").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_Pos_Px").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_Pos_Py").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_Pos_Pz").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_v0 + "_Pos_E").c_str(), true);

    fTree_PackedEvents->SetBranchAddress((name_v0 + "_Entry").c_str(), &vec_v0s.Entry);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_X").c_str(), &vec_v0s.X);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_Y").c_str(), &vec_v0s.Y);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_Z").c_str(), &vec_v0s.Z);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_Px").c_str(), &vec_v0s.Px);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_Py").c_str(), &vec_v0s.Py);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_Pz").c_str(), &vec_v0s.Pz);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_E").c_str(), &vec_v0s.E);

    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaX2").c_str(), &vec_v0s.Sigma.X2);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaXY").c_str(), &vec_v0s.Sigma.XY);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaY2").c_str(), &vec_v0s.Sigma.Y2);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaXZ").c_str(), &vec_v0s.Sigma.XZ);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaYZ").c_str(), &vec_v0s.Sigma.YZ);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaZ2").c_str(), &vec_v0s.Sigma.Z2);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaXPx").c_str(), &vec_v0s.Sigma.XPx);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaYPx").c_str(), &vec_v0s.Sigma.YPx);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaZPx").c_str(), &vec_v0s.Sigma.ZPx);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaPx2").c_str(), &vec_v0s.Sigma.Px2);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaXPy").c_str(), &vec_v0s.Sigma.XPy);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaYPy").c_str(), &vec_v0s.Sigma.YPy);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaZPy").c_str(), &vec_v0s.Sigma.ZPy);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaPxPy").c_str(), &vec_v0s.Sigma.PxPy);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaPy2").c_str(), &vec_v0s.Sigma.Py2);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaXPz").c_str(), &vec_v0s.Sigma.XPz);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaYPz").c_str(), &vec_v0s.Sigma.YPz);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaZPz").c_str(), &vec_v0s.Sigma.ZPz);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaPxPz").c_str(), &vec_v0s.Sigma.PxPz);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaPyPz").c_str(), &vec_v0s.Sigma.PyPz);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaPz2").c_str(), &vec_v0s.Sigma.Pz2);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaXE").c_str(), &vec_v0s.Sigma.XE);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaYE").c_str(), &vec_v0s.Sigma.YE);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaZE").c_str(), &vec_v0s.Sigma.ZE);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaPxE").c_str(), &vec_v0s.Sigma.PxE);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaPyE").c_str(), &vec_v0s.Sigma.PyE);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaPzE").c_str(), &vec_v0s.Sigma.PzE);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaE2").c_str(), &vec_v0s.Sigma.E2);

    fTree_PackedEvents->SetBranchAddress((name_v0 + "_Neg_Entry").c_str(), &vec_v0s.Neg.Entry);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_Neg_X").c_str(), &vec_v0s.Neg.X);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_Neg_Y").c_str(), &vec_v0s.Neg.Y);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_Neg_Z").c_str(), &vec_v0s.Neg.Z);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_Neg_Px").c_str(), &vec_v0s.Neg.Px);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_Neg_Py").c_str(), &vec_v0s.Neg.Py);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_Neg_Pz").c_str(), &vec_v0s.Neg.Pz);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_Neg_E").c_str(), &vec_v0s.Neg.E);

    fTree_PackedEvents->SetBranchAddress((name_v0 + "_Pos_Entry").c_str(), &vec_v0s.Pos.Entry);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_Pos_X").c_str(), &vec_v0s.Pos.X);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_Pos_Y").c_str(), &vec_v0s.Pos.Y);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_Pos_Z").c_str(), &vec_v0s.Pos.Z);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_Pos_Px").c_str(), &vec_v0s.Pos.Px);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_Pos_Py").c_str(), &vec_v0s.Pos.Py);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_Pos_Pz").c_str(), &vec_v0s.Pos.Pz);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_Pos_E").c_str(), &vec_v0s.Pos.E);

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Finder::ConnectBranches_Tracks(const std::string& name_part, PackedEvents::Tracks& vec_tracks) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif
    fTree_PackedEvents->SetBranchStatus((name_part + "_Entry").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_X").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_Y").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_Z").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_Px").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_Py").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_Pz").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_E").c_str(), true);

    fTree_PackedEvents->SetBranchStatus((name_part + "_SigmaX2").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_SigmaXY").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_SigmaY2").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_SigmaXZ").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_SigmaYZ").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_SigmaZ2").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_SigmaXPx").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_SigmaYPx").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_SigmaZPx").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_SigmaPx2").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_SigmaXPy").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_SigmaYPy").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_SigmaZPy").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_SigmaPxPy").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_SigmaPy2").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_SigmaXPz").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_SigmaYPz").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_SigmaZPz").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_SigmaPxPz").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_SigmaPyPz").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_SigmaPz2").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_SigmaXE").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_SigmaYE").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_SigmaZE").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_SigmaPxE").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_SigmaPyE").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_SigmaPzE").c_str(), true);
    fTree_PackedEvents->SetBranchStatus((name_part + "_SigmaE2").c_str(), true);

    fTree_PackedEvents->SetBranchAddress((name_part + "_Entry").c_str(), &vec_tracks.Entry);
    fTree_PackedEvents->SetBranchAddress((name_part + "_X").c_str(), &vec_tracks.X);
    fTree_PackedEvents->SetBranchAddress((name_part + "_Y").c_str(), &vec_tracks.Y);
    fTree_PackedEvents->SetBranchAddress((name_part + "_Z").c_str(), &vec_tracks.Z);
    fTree_PackedEvents->SetBranchAddress((name_part + "_Px").c_str(), &vec_tracks.Px);
    fTree_PackedEvents->SetBranchAddress((name_part + "_Py").c_str(), &vec_tracks.Py);
    fTree_PackedEvents->SetBranchAddress((name_part + "_Pz").c_str(), &vec_tracks.Pz);
    fTree_PackedEvents->SetBranchAddress((name_part + "_E").c_str(), &vec_tracks.E);

    fTree_PackedEvents->SetBranchAddress((name_part + "_SigmaX2").c_str(), &vec_tracks.Sigma.X2);
    fTree_PackedEvents->SetBranchAddress((name_part + "_SigmaXY").c_str(), &vec_tracks.Sigma.XY);
    fTree_PackedEvents->SetBranchAddress((name_part + "_SigmaY2").c_str(), &vec_tracks.Sigma.Y2);
    fTree_PackedEvents->SetBranchAddress((name_part + "_SigmaXZ").c_str(), &vec_tracks.Sigma.XZ);
    fTree_PackedEvents->SetBranchAddress((name_part + "_SigmaYZ").c_str(), &vec_tracks.Sigma.YZ);
    fTree_PackedEvents->SetBranchAddress((name_part + "_SigmaZ2").c_str(), &vec_tracks.Sigma.Z2);
    fTree_PackedEvents->SetBranchAddress((name_part + "_SigmaXPx").c_str(), &vec_tracks.Sigma.XPx);
    fTree_PackedEvents->SetBranchAddress((name_part + "_SigmaYPx").c_str(), &vec_tracks.Sigma.YPx);
    fTree_PackedEvents->SetBranchAddress((name_part + "_SigmaZPx").c_str(), &vec_tracks.Sigma.ZPx);
    fTree_PackedEvents->SetBranchAddress((name_part + "_SigmaPx2").c_str(), &vec_tracks.Sigma.Px2);
    fTree_PackedEvents->SetBranchAddress((name_part + "_SigmaXPy").c_str(), &vec_tracks.Sigma.XPy);
    fTree_PackedEvents->SetBranchAddress((name_part + "_SigmaYPy").c_str(), &vec_tracks.Sigma.YPy);
    fTree_PackedEvents->SetBranchAddress((name_part + "_SigmaZPy").c_str(), &vec_tracks.Sigma.ZPy);
    fTree_PackedEvents->SetBranchAddress((name_part + "_SigmaPxPy").c_str(), &vec_tracks.Sigma.PxPy);
    fTree_PackedEvents->SetBranchAddress((name_part + "_SigmaPy2").c_str(), &vec_tracks.Sigma.Py2);
    fTree_PackedEvents->SetBranchAddress((name_part + "_SigmaXPz").c_str(), &vec_tracks.Sigma.XPz);
    fTree_PackedEvents->SetBranchAddress((name_part + "_SigmaYPz").c_str(), &vec_tracks.Sigma.YPz);
    fTree_PackedEvents->SetBranchAddress((name_part + "_SigmaZPz").c_str(), &vec_tracks.Sigma.ZPz);
    fTree_PackedEvents->SetBranchAddress((name_part + "_SigmaPxPz").c_str(), &vec_tracks.Sigma.PxPz);
    fTree_PackedEvents->SetBranchAddress((name_part + "_SigmaPyPz").c_str(), &vec_tracks.Sigma.PyPz);
    fTree_PackedEvents->SetBranchAddress((name_part + "_SigmaPz2").c_str(), &vec_tracks.Sigma.Pz2);
    fTree_PackedEvents->SetBranchAddress((name_part + "_SigmaXE").c_str(), &vec_tracks.Sigma.XE);
    fTree_PackedEvents->SetBranchAddress((name_part + "_SigmaYE").c_str(), &vec_tracks.Sigma.YE);
    fTree_PackedEvents->SetBranchAddress((name_part + "_SigmaZE").c_str(), &vec_tracks.Sigma.ZE);
    fTree_PackedEvents->SetBranchAddress((name_part + "_SigmaPxE").c_str(), &vec_tracks.Sigma.PxE);
    fTree_PackedEvents->SetBranchAddress((name_part + "_SigmaPyE").c_str(), &vec_tracks.Sigma.PyE);
    fTree_PackedEvents->SetBranchAddress((name_part + "_SigmaPzE").c_str(), &vec_tracks.Sigma.PzE);
    fTree_PackedEvents->SetBranchAddress((name_part + "_SigmaE2").c_str(), &vec_tracks.Sigma.E2);

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

// ## OUTPUT ZONE ## //

bool Finder::PrepareOutputFile() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif
    fOutputFile = std::unique_ptr<TFile>(TFile::Open(fSettings.PathOutputFile.c_str(), "RECREATE"));
    if (!fOutputFile) {
        std::cerr << "TFile \"" << fSettings.PathOutputFile << "\" couldn't be created" << '\n';
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

    std::string tree_name;
    switch (GetReactionChannel()) {
        case ReactionChannel::A:
            tree_name = static_cast<std::string>(Name::ChannelA);
            break;
        case ReactionChannel::D:
            tree_name = static_cast<std::string>(Name::ChannelD);
            break;
        case ReactionChannel::E:
            tree_name = static_cast<std::string>(Name::ChannelE);
            break;
        case ReactionChannel::H:
            tree_name = static_cast<std::string>(Name::ChannelH);
            break;
        case ReactionChannel::AntiA:
            tree_name = static_cast<std::string>(Name::AntiA);
            break;
        case ReactionChannel::AntiD:
            tree_name = static_cast<std::string>(Name::AntiD);
            break;
        case ReactionChannel::AntiE:
            tree_name = static_cast<std::string>(Name::AntiE);
            break;
        case ReactionChannel::AntiH:
            tree_name = static_cast<std::string>(Name::AntiH);
            break;
        default:
            tree_name = "";
            break;
    }

    fOutputTree = std::make_unique<TTree>(tree_name.c_str(), "Final results");
    if (!fOutputTree) {
        std::cerr << "TTree \"" << tree_name << "\" couldn't be created" << '\n';
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
    fOutputTree->Branch("Entry", &out_branches.Entry);

    fOutputTree->Branch("X", &out_branches.X);
    fOutputTree->Branch("Y", &out_branches.Y);
    fOutputTree->Branch("Z", &out_branches.Z);
    fOutputTree->Branch("Px", &out_branches.Px);
    fOutputTree->Branch("Py", &out_branches.Py);
    fOutputTree->Branch("Pz", &out_branches.Pz);
    fOutputTree->Branch("E", &out_branches.E);
    fOutputTree->Branch("E_MinusNucleon", &out_branches.E_MinusNucleon);

    fOutputTree->Branch("V0A_Entry", &out_branches.V0A.Entry);
    fOutputTree->Branch("V0A_Neg_Entry", &out_branches.V0A_Neg_Entry);
    fOutputTree->Branch("V0A_Pos_Entry", &out_branches.V0A_Pos_Entry);
    fOutputTree->Branch("V0A_X", &out_branches.V0A.X);
    fOutputTree->Branch("V0A_Y", &out_branches.V0A.Y);
    fOutputTree->Branch("V0A_Z", &out_branches.V0A.Z);
    fOutputTree->Branch("V0A_Px", &out_branches.V0A.Px);
    fOutputTree->Branch("V0A_Py", &out_branches.V0A.Py);
    fOutputTree->Branch("V0A_Pz", &out_branches.V0A.Pz);
    fOutputTree->Branch("V0A_E", &out_branches.V0A.E);
    fOutputTree->Branch("V0A_DecayX", &out_branches.V0A_DecayX);
    fOutputTree->Branch("V0A_DecayY", &out_branches.V0A_DecayY);
    fOutputTree->Branch("V0A_DecayZ", &out_branches.V0A_DecayZ);

    fOutputTree->Branch("V0B_Entry", &out_branches.V0B.Entry);
    fOutputTree->Branch("V0B_Neg_Entry", &out_branches.V0B_Neg_Entry);
    fOutputTree->Branch("V0B_Pos_Entry", &out_branches.V0B_Pos_Entry);
    fOutputTree->Branch("V0B_X", &out_branches.V0B.X);
    fOutputTree->Branch("V0B_Y", &out_branches.V0B.Y);
    fOutputTree->Branch("V0B_Z", &out_branches.V0B.Z);
    fOutputTree->Branch("V0B_Px", &out_branches.V0B.Px);
    fOutputTree->Branch("V0B_Py", &out_branches.V0B.Py);
    fOutputTree->Branch("V0B_Pz", &out_branches.V0B.Pz);
    fOutputTree->Branch("V0B_E", &out_branches.V0B.E);
    fOutputTree->Branch("V0B_DecayX", &out_branches.V0B_DecayX);
    fOutputTree->Branch("V0B_DecayY", &out_branches.V0B_DecayY);
    fOutputTree->Branch("V0B_DecayZ", &out_branches.V0B_DecayZ);

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
    fOutputTree->Branch("Entry", &out_branches.Entry);

    fOutputTree->Branch("X", &out_branches.X);
    fOutputTree->Branch("Y", &out_branches.Y);
    fOutputTree->Branch("Z", &out_branches.Z);
    fOutputTree->Branch("Px", &out_branches.Px);
    fOutputTree->Branch("Py", &out_branches.Py);
    fOutputTree->Branch("Pz", &out_branches.Pz);
    fOutputTree->Branch("E", &out_branches.E);
    fOutputTree->Branch("E_MinusNucleon", &out_branches.E_MinusNucleon);

    fOutputTree->Branch("V0_Entry", &out_branches.V0.Entry);
    fOutputTree->Branch("V0_X", &out_branches.V0.X);
    fOutputTree->Branch("V0_Y", &out_branches.V0.Y);
    fOutputTree->Branch("V0_Z", &out_branches.V0.Z);
    fOutputTree->Branch("V0_Px", &out_branches.V0.Px);
    fOutputTree->Branch("V0_Py", &out_branches.V0.Py);
    fOutputTree->Branch("V0_Pz", &out_branches.V0.Pz);
    fOutputTree->Branch("V0_E", &out_branches.V0.E);

    fOutputTree->Branch("Kaon_Entry", &out_branches.Kaon.Entry);
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

// ## Intermediate ZONE ## //

void Finder::UnpackV0s(PdgCode pdg_code) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    // determine rules based on particle pdg code //
    const PackedEvents::V0s* in;
    std::vector<KF::V0>* vec_kf;
    switch (pdg_code) {
        case PdgCode::AntiLambda:
            in = &fPacked_AntiLambdas;
            vec_kf = &fKF_AntiLambdas;
            break;
        case PdgCode::Lambda:
            in = &fPacked_Lambdas;
            vec_kf = &fKF_Lambdas;
            break;
        case PdgCode::KaonZeroShort:
            in = &fPacked_KaonsZeroShort;
            vec_kf = &fKF_KaonsZeroShort;
            break;
        default:
            return;
    }

    // loop over v0s //
    size_t n_v0s{in->Entry->size()};
    vec_kf->reserve(n_v0s);
    for (size_t idx_v0{0}; idx_v0 < n_v0s; ++idx_v0) {
        auto params = KF::UnpackParams(*in, idx_v0);
        auto cov_matrix = KF::UnpackCovMatrix(*in, idx_v0);
        vec_kf->emplace_back(params, cov_matrix, 0, idx_v0, in->Neg.Entry->at(idx_v0), in->Pos.Entry->at(idx_v0), static_cast<int>(pdg_code));
    }  // end of loop over v0s

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Finder::UnpackTracks(PdgCode pdg_code) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    // determine rules based on particle pdg code //
    const PackedEvents::Tracks* in;
    std::vector<KF::Track>* vec_kf;
    int charge;
    switch (pdg_code) {
        case PdgCode::NegKaon:
            in = &fPacked_NegKaons;
            vec_kf = &fKF_NegKaons;
            charge = -1;
            break;
        case PdgCode::PosKaon:
            in = &fPacked_PosKaons;
            vec_kf = &fKF_PosKaons;
            charge = +1;
            break;
        case PdgCode::PiMinus:
            in = &fPacked_PiMinus;
            vec_kf = &fKF_PiMinus;
            charge = -1;
            break;
        case PdgCode::PiPlus:
            in = &fPacked_PiPlus;
            vec_kf = &fKF_PiPlus;
            charge = +1;
            break;
        default:
            return;
    }

    // loop over tracks //
    size_t n_tracks{in->Entry->size()};
    vec_kf->reserve(n_tracks);
    for (size_t idx{0}; idx < n_tracks; ++idx) {
        auto params = KF::UnpackParams(*in, idx);
        auto cov_matrix = KF::UnpackCovMatrix(*in, idx);
        vec_kf->emplace_back(params, cov_matrix, charge, idx);
    }  // end of loop over tracks
#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

// ## Channel A ZONE ## //

void Finder::FindSexaquarks_ChannelA(bool anti_channel) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    const std::vector<KF::V0>* KF_Lambdas{anti_channel ? &fKF_Lambdas : &fKF_AntiLambdas};

    // loop over all possible pairs of (anti)lambda + K0S //
    size_t n_sexa{0};
    for (const auto& v0a : *KF_Lambdas) {
        for (const auto& v0b : fKF_KaonsZeroShort) {

            // sanity check //
            std::set<size_t> unique_track_entries{v0a.idx_neg, v0a.idx_pos, v0b.idx_neg, v0b.idx_pos};
            if (unique_track_entries.size() < 4) continue;

            // fit //
            KF::ChannelA sexa;
            sexa.AddDaughter(v0a, fInput_Event.MagneticField);
            sexa.AddDaughter(v0b, fInput_Event.MagneticField);

            // apply cuts //
            if (!PassesCuts(sexa)) continue;

            // add info //
            sexa.SetIndices({n_sexa, v0a.idx, v0a.idx_neg, v0a.idx_pos, v0b.idx, v0b.idx_neg, v0b.idx_pos});
            sexa.SetAdditionalInfo_V0A({v0a.X(), v0a.Y(), v0a.Z()}, v0a.E());
            sexa.SetAdditionalInfo_V0B({v0b.X(), v0b.Y(), v0b.Z()}, v0b.E());
            sexa.Nucleon_Mass = PdgMass::Neutron;
#ifdef T2S_DEBUG
            std::cout << "   " << __FUNCTION__ << " :: idx,v0a,v0a_neg,v0a_pos,v0b,v0b_neg,v0b_pos="                            //
                      << sexa.idx << "," << sexa.idx_lambda << "," << sexa.idx_lambda_neg << "," << sexa.idx_lambda_pos << ","  //
                      << sexa.idx_k0s << "," << sexa.idx_k0s_neg << "," << sexa.idx_k0s_pos;
            std::cout << ";x,y,z=" << sexa.X() << "," << sexa.Y() << "," << sexa.Z();
            std::cout << ";x,y,z(v0a)=" << sexa.PCA_V0A()[0] << "," << sexa.PCA_V0A()[1] << "," << sexa.PCA_V0A()[2];
            std::cout << ";x,y,z(v0b)=" << sexa.PCA_V0B()[0] << "," << sexa.PCA_V0B()[1] << "," << sexa.PCA_V0B()[2];
            std::cout << ";mass=" << sexa.Mass();
            std::cout << ";dca_dau=" << sexa.DCA_V0s();
            std::cout << ";radius=" << sexa.Radius();
            std::cout << ";dca_neg=" << sexa.DCA_V0A();
            std::cout << ";dca_pos=" << sexa.DCA_V0B();
            std::cout << ";pt=" << sexa.Pt();
            std::cout << ";eta=" << sexa.Eta();
            std::cout << ";cpa_pv=" << sexa.CPA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv) << '\n';
#endif

            // store //
            Store(sexa);
            ++n_sexa;
        }
    }

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

bool Finder::PassesCuts(const KF::ChannelA& sexa) const {

    if (sexa.Radius() < Cuts::ChannelA::Min_Radius || sexa.Radius() > Cuts::ChannelA::Max_Radius) return false;
    if (sexa.DecayLength_V0A() > Cuts::ChannelA::Max_DecayLengthLa) return false;
    if (sexa.DecayLength_V0B() > Cuts::ChannelA::Max_DecayLengthK0) return false;
    if (sexa.AbsRapidity_MinusNucleon() > Cuts::ChannelA::AbsMax_Rapidity) return false;  // PENDING: kinematic cut, affected by Fermi motion
    if (sexa.Mass_MinusNucleon() < Cuts::ChannelA::Min_MassMinusNucleon || sexa.Mass_MinusNucleon() > Cuts::ChannelA::Max_MassMinusNucleon)
        return false;  // PENDING: kinematic cut, affected by Fermi motion
    if (sexa.CPA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv) < Cuts::ChannelA::Min_CPAwrtPV ||
        sexa.CPA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv) > Cuts::ChannelA::Max_CPAwrtPV)
        return false;  // PENDING: kinematic cut, affected by Fermi motion
    // if (sexa.DCA_V0A_Neg() > Cuts::ChannelA::Max_DCALaNegSV) return false;  // PENDING: not trivial
    // if (sexa.DCA_V0A_Pos() > Cuts::ChannelA::Max_DCALaPosSV) return false;  // PENDING: not trivial
    // if (sexa.DCA_V0B_Neg() > Cuts::ChannelA::Max_DCAK0NegSV) return false;  // PENDING: not trivial
    // if (sexa.DCA_V0B_Pos() > Cuts::ChannelA::Max_DCAK0PosSV) return false;  // PENDING: not trivial
    if (sexa.DCA_V0A() > Cuts::ChannelA::Max_DCALaSV) return false;
    if (sexa.DCA_V0B() > Cuts::ChannelA::Max_DCAK0SV) return false;
    if (sexa.DCA_V0s() > Cuts::ChannelA::Max_DCAbtwV0s) return false;

    return true;
}

void Finder::Store(const KF::ChannelA& sexa) {

    fOutput_ChannelA.RunNumber = fInput_Event.RunNumber;
    fOutput_ChannelA.DirNumber = fInput_Event.DirNumber;
    if (!IsMC()) fOutput_ChannelA.DirNumberB = fInput_Event.DirNumberB;
    fOutput_ChannelA.EventNumber = fInput_Event.EventNumber;
    fOutput_ChannelA.Entry = sexa.idx;

    fOutput_ChannelA.X = static_cast<float>(sexa.X());
    fOutput_ChannelA.Y = static_cast<float>(sexa.Y());
    fOutput_ChannelA.Z = static_cast<float>(sexa.Z());
    fOutput_ChannelA.Px = static_cast<float>(sexa.Px());
    fOutput_ChannelA.Py = static_cast<float>(sexa.Py());
    fOutput_ChannelA.Pz = static_cast<float>(sexa.Pz());
    fOutput_ChannelA.E = static_cast<float>(sexa.E());
    fOutput_ChannelA.E_MinusNucleon = static_cast<float>(sexa.E_MinusNucleon());

    fOutput_ChannelA.V0A.Entry = sexa.idx_lambda;
    fOutput_ChannelA.V0A_Neg_Entry = sexa.idx_lambda_neg;
    fOutput_ChannelA.V0A_Pos_Entry = sexa.idx_lambda_pos;
    fOutput_ChannelA.V0A.X = static_cast<float>(sexa.PCA_V0A()[0]);
    fOutput_ChannelA.V0A.Y = static_cast<float>(sexa.PCA_V0A()[1]);
    fOutput_ChannelA.V0A.Z = static_cast<float>(sexa.PCA_V0A()[2]);
    fOutput_ChannelA.V0A.Px = static_cast<float>(sexa.Mom_V0A()[0]);
    fOutput_ChannelA.V0A.Py = static_cast<float>(sexa.Mom_V0A()[1]);
    fOutput_ChannelA.V0A.Pz = static_cast<float>(sexa.Mom_V0A()[2]);
    fOutput_ChannelA.V0A.E = static_cast<float>(sexa.V0A_Energy);
    fOutput_ChannelA.V0A_DecayX = static_cast<float>(sexa.V0A_DecayVtx[0]);
    fOutput_ChannelA.V0A_DecayY = static_cast<float>(sexa.V0A_DecayVtx[1]);
    fOutput_ChannelA.V0A_DecayZ = static_cast<float>(sexa.V0A_DecayVtx[2]);

    fOutput_ChannelA.V0B.Entry = sexa.idx_k0s;
    fOutput_ChannelA.V0B_Neg_Entry = sexa.idx_k0s_neg;
    fOutput_ChannelA.V0B_Pos_Entry = sexa.idx_k0s_pos;
    fOutput_ChannelA.V0B.X = static_cast<float>(sexa.PCA_V0B()[0]);
    fOutput_ChannelA.V0B.Y = static_cast<float>(sexa.PCA_V0B()[1]);
    fOutput_ChannelA.V0B.Z = static_cast<float>(sexa.PCA_V0B()[2]);
    fOutput_ChannelA.V0B.Px = static_cast<float>(sexa.Mom_V0B()[0]);
    fOutput_ChannelA.V0B.Py = static_cast<float>(sexa.Mom_V0B()[1]);
    fOutput_ChannelA.V0B.Pz = static_cast<float>(sexa.Mom_V0B()[2]);
    fOutput_ChannelA.V0B.E = static_cast<float>(sexa.V0B_Energy);
    fOutput_ChannelA.V0B_DecayX = static_cast<float>(sexa.V0B_DecayVtx[0]);
    fOutput_ChannelA.V0B_DecayY = static_cast<float>(sexa.V0B_DecayVtx[1]);
    fOutput_ChannelA.V0B_DecayZ = static_cast<float>(sexa.V0B_DecayVtx[2]);

    fOutputTree->Fill();
}

// ## Channel D ZONE ## //

void Finder::FindSexaquarks_ChannelD(bool anti_channel) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    const std::vector<KF::V0>* KF_Lambdas{anti_channel ? &fKF_Lambdas : &fKF_AntiLambdas};
    const std::vector<KF::Track>* KF_Kaons{anti_channel ? &fKF_NegKaons : &fKF_PosKaons};

    // loop over all possible pairs of (anti)lambda + (pos/neg)kaon //
    size_t n_sexa{0};
    for (const auto& v0 : *KF_Lambdas) {
        for (const auto& bach : *KF_Kaons) {
            // sanity check //
            std::set<size_t> unique_track_entries{v0.idx_neg, v0.idx_pos, bach.idx};
            if (unique_track_entries.size() < 3) continue;
            // fit //
            KF::ChannelD sexa;
            sexa.AddDaughter(v0, fInput_Event.MagneticField);
            sexa.AddDaughter(bach, fInput_Event.MagneticField);
            // apply cuts //
            if (!PassesCuts(sexa)) continue;
            // add info //
            sexa.SetIndices({n_sexa, v0.idx, v0.idx_neg, v0.idx_pos, bach.idx});
            sexa.SetAdditionalInfo_V0({v0.X(), v0.Y(), v0.Z()}, v0.E());
            sexa.Kaon_Energy = bach.E();
            sexa.Nucleon_Mass = PdgMass::Proton;
            // store //
            Store(sexa);
            ++n_sexa;
        }
    }

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

bool Finder::PassesCuts(const KF::ChannelD& sexa) const {

    if (sexa.Radius() < Cuts::ChannelD::Min_Radius || sexa.Radius() > Cuts::ChannelD::Max_Radius) return false;
    if (sexa.AbsRapidity_MinusNucleon() > Cuts::ChannelD::AbsMax_Rapidity) return false;  // PENDING: kinematics, affected by Fermi motion
    if (sexa.CPA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv) < Cuts::ChannelD::Min_CPAwrtPV ||
        sexa.CPA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv) > Cuts::ChannelD::Max_CPAwrtPV)
        return false;  // PENDING: kinematics, affected by Fermi motion
    if (sexa.DCA_V0() > Cuts::ChannelD::Max_DCALaSV) return false;
    if (sexa.DCA_Kaon() > Cuts::ChannelD::Max_DCAKaSV) return false;
    // if (sexa.DCA_V0_Neg() > Cuts::ChannelD::Max_DCALaNegSV) return false; // PENDING: not trivial
    // if (sexa.DCA_V0_Pos() > Cuts::ChannelD::Max_DCALaPosSV) return false; // PENDING: not trivial
    if (sexa.DCA_Kaon_V0() > Cuts::ChannelD::Max_DCAKaLa) return false;
    // if (sexa.DCA_V0_Neg_Kaon() > Cuts::ChannelD::Max_DCALaNegKa) return false; // PENDING: not trivial
    // if (sexa.DCA_V0_Pos_Kaon() > Cuts::ChannelD::Max_DCALaPosKa) return false; // PENDING: not trivial

    return true;
}

void Finder::Store(const KF::ChannelD& sexa) {

    fOutput_ChannelD.RunNumber = fInput_Event.RunNumber;
    fOutput_ChannelD.DirNumber = fInput_Event.DirNumber;
    if (!IsMC()) fOutput_ChannelD.DirNumberB = fInput_Event.DirNumberB;
    fOutput_ChannelD.EventNumber = fInput_Event.EventNumber;
    fOutput_ChannelD.Entry = sexa.idx;

    fOutput_ChannelD.X = static_cast<float>(sexa.X());
    fOutput_ChannelD.Y = static_cast<float>(sexa.Y());
    fOutput_ChannelD.Z = static_cast<float>(sexa.Z());
    fOutput_ChannelD.Px = static_cast<float>(sexa.Px());
    fOutput_ChannelD.Py = static_cast<float>(sexa.Py());
    fOutput_ChannelD.Pz = static_cast<float>(sexa.Pz());
    fOutput_ChannelD.E = static_cast<float>(sexa.E());
    fOutput_ChannelD.E_MinusNucleon = static_cast<float>(sexa.E_MinusNucleon());

    fOutput_ChannelD.V0.Entry = sexa.idx_lambda;
    fOutput_ChannelD.V0_Neg_Entry = sexa.idx_lambda_neg;
    fOutput_ChannelD.V0_Pos_Entry = sexa.idx_lambda_pos;
    fOutput_ChannelD.V0.X = static_cast<float>(sexa.PCA_First()[0]);
    fOutput_ChannelD.V0.Y = static_cast<float>(sexa.PCA_First()[1]);
    fOutput_ChannelD.V0.Z = static_cast<float>(sexa.PCA_First()[2]);
    fOutput_ChannelD.V0.Px = static_cast<float>(sexa.Mom_First()[0]);
    fOutput_ChannelD.V0.Py = static_cast<float>(sexa.Mom_First()[1]);
    fOutput_ChannelD.V0.Pz = static_cast<float>(sexa.Mom_First()[2]);
    fOutput_ChannelD.V0.E = static_cast<float>(sexa.V0_Energy);
    fOutput_ChannelD.V0_DecayX = static_cast<float>(sexa.V0_DecayVtx[0]);
    fOutput_ChannelD.V0_DecayY = static_cast<float>(sexa.V0_DecayVtx[1]);
    fOutput_ChannelD.V0_DecayZ = static_cast<float>(sexa.V0_DecayVtx[2]);

    fOutput_ChannelD.Kaon.Entry = sexa.idx_kaon;
    fOutput_ChannelD.Kaon.X = static_cast<float>(sexa.PCA_Second()[0]);
    fOutput_ChannelD.Kaon.Y = static_cast<float>(sexa.PCA_Second()[1]);
    fOutput_ChannelD.Kaon.Z = static_cast<float>(sexa.PCA_Second()[2]);
    fOutput_ChannelD.Kaon.Px = static_cast<float>(sexa.Mom_Second()[0]);
    fOutput_ChannelD.Kaon.Py = static_cast<float>(sexa.Mom_Second()[1]);
    fOutput_ChannelD.Kaon.Pz = static_cast<float>(sexa.Mom_Second()[2]);
    fOutput_ChannelD.Kaon.E = static_cast<float>(sexa.Kaon_Energy);

    if (IsMC()) {
        // PENDING
    }

    fOutputTree->Fill();
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

void Finder::EndOfEvent() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    // clear temporary containers //
    switch (GetReactionChannel()) {
        // standard channels //
        case ReactionChannel::A:
            fPacked_AntiLambdas.Clear();
            fPacked_KaonsZeroShort.Clear();
            fKF_AntiLambdas.clear();
            fKF_KaonsZeroShort.clear();
            break;
        case ReactionChannel::D:
            fPacked_AntiLambdas.Clear();
            fPacked_PosKaons.Clear();
            fKF_AntiLambdas.clear();
            fKF_PosKaons.clear();
            break;
        case ReactionChannel::E:
            fPacked_AntiLambdas.Clear();
            fPacked_PosKaons.Clear();
            fPacked_PiMinus.Clear();
            fPacked_PiPlus.Clear();
            fKF_AntiLambdas.clear();
            fKF_PosKaons.clear();
            fKF_PiMinus.clear();
            fKF_PiPlus.clear();
            break;
        case ReactionChannel::H:
            fPacked_PosKaons.Clear();
            fKF_PosKaons.clear();
            break;
        // anti-channels //
        case ReactionChannel::AntiA:
            fPacked_Lambdas.Clear();
            fPacked_KaonsZeroShort.Clear();
            fKF_Lambdas.clear();
            fKF_KaonsZeroShort.clear();
            break;
        case ReactionChannel::AntiD:
            fPacked_Lambdas.Clear();
            fPacked_NegKaons.Clear();
            fKF_Lambdas.clear();
            fKF_NegKaons.clear();
            break;
        case ReactionChannel::AntiE:
            fPacked_Lambdas.Clear();
            fPacked_NegKaons.Clear();
            fPacked_PiMinus.Clear();
            fPacked_PiPlus.Clear();
            fKF_Lambdas.clear();
            fKF_NegKaons.clear();
            fKF_PiMinus.clear();
            fKF_PiPlus.clear();
            break;
        case ReactionChannel::AntiH:
            fPacked_NegKaons.Clear();
            fKF_NegKaons.clear();
            break;
        default:
            break;
    }

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Finder::EndOfAnalysis() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

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
