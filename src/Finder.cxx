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
            std::cerr << "Couldn't add TFile \"" << path << "\"" << '\n';
        }
    }
    if (!fTree_PackedEvents->GetEntries()) {
        std::cerr << "Couldn't manage to read any entry." << '\n';
        return false;
    }
    std::cout << "TChain \"PackedEvents\" loaded successfully with " << fTree_PackedEvents->GetNtrees() << " trees and "
              << fTree_PackedEvents->GetEntries() << " total entries." << '\n';

    ConnectInputBranches();

    if (!PrepareOutputFile()) return false;
    if (!PrepareOutputTree()) return false;
    CreateOutputBranches();

    std::cout << "Finder initialized successfully" << '\n';

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
        // PENDING
    }
#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Finder::ConnectBranches_V0s(const std::string& name_v0, PackedEvents::V0s& vec_v0s) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

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

    fTree_PackedEvents->SetBranchAddress((name_v0 + "_Xv").c_str(), &vec_v0s.Xv);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_Yv").c_str(), &vec_v0s.Yv);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_Zv").c_str(), &vec_v0s.Zv);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_Px").c_str(), &vec_v0s.Px);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_Py").c_str(), &vec_v0s.Py);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_Pz").c_str(), &vec_v0s.Pz);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_E").c_str(), &vec_v0s.E);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaX2").c_str(), &vec_v0s.SigmaX2);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaXY").c_str(), &vec_v0s.SigmaXY);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaY2").c_str(), &vec_v0s.SigmaY2);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaXZ").c_str(), &vec_v0s.SigmaXZ);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaYZ").c_str(), &vec_v0s.SigmaYZ);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaZ2").c_str(), &vec_v0s.SigmaZ2);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaXPx").c_str(), &vec_v0s.SigmaXPx);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaYPx").c_str(), &vec_v0s.SigmaYPx);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaZPx").c_str(), &vec_v0s.SigmaZPx);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaPx2").c_str(), &vec_v0s.SigmaPx2);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaXPy").c_str(), &vec_v0s.SigmaXPy);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaYPy").c_str(), &vec_v0s.SigmaYPy);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaZPy").c_str(), &vec_v0s.SigmaZPy);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaPxPy").c_str(), &vec_v0s.SigmaPxPy);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaPy2").c_str(), &vec_v0s.SigmaPy2);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaXPz").c_str(), &vec_v0s.SigmaXPz);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaYPz").c_str(), &vec_v0s.SigmaYPz);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaZPz").c_str(), &vec_v0s.SigmaZPz);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaPxPz").c_str(), &vec_v0s.SigmaPxPz);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaPyPz").c_str(), &vec_v0s.SigmaPyPz);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaPz2").c_str(), &vec_v0s.SigmaPz2);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaXE").c_str(), &vec_v0s.SigmaXE);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaYE").c_str(), &vec_v0s.SigmaYE);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaZE").c_str(), &vec_v0s.SigmaZE);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaPxE").c_str(), &vec_v0s.SigmaPxE);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaPyE").c_str(), &vec_v0s.SigmaPyE);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaPzE").c_str(), &vec_v0s.SigmaPzE);
    fTree_PackedEvents->SetBranchAddress((name_v0 + "_SigmaE2").c_str(), &vec_v0s.SigmaE2);

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Finder::ConnectBranches_Tracks(const std::string& name_part, PackedEvents::Tracks& vec_tracks) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif
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

    fTree_PackedEvents->SetBranchAddress("Track_Xv", &vec_tracks.Xv);
    fTree_PackedEvents->SetBranchAddress("Track_Yv", &vec_tracks.Yv);
    fTree_PackedEvents->SetBranchAddress("Track_Zv", &vec_tracks.Zv);
    fTree_PackedEvents->SetBranchAddress("Track_Px", &vec_tracks.Px);
    fTree_PackedEvents->SetBranchAddress("Track_Py", &vec_tracks.Py);
    fTree_PackedEvents->SetBranchAddress("Track_Pz", &vec_tracks.Pz);
    fTree_PackedEvents->SetBranchAddress("Track_E", &vec_tracks.E);

    fTree_PackedEvents->SetBranchAddress("Track_SigmaX2", &vec_tracks.SigmaX2);
    fTree_PackedEvents->SetBranchAddress("Track_SigmaXY", &vec_tracks.SigmaXY);
    fTree_PackedEvents->SetBranchAddress("Track_SigmaY2", &vec_tracks.SigmaY2);
    fTree_PackedEvents->SetBranchAddress("Track_SigmaXZ", &vec_tracks.SigmaXZ);
    fTree_PackedEvents->SetBranchAddress("Track_SigmaYZ", &vec_tracks.SigmaYZ);
    fTree_PackedEvents->SetBranchAddress("Track_SigmaZ2", &vec_tracks.SigmaZ2);
    fTree_PackedEvents->SetBranchAddress("Track_SigmaXPx", &vec_tracks.SigmaXPx);
    fTree_PackedEvents->SetBranchAddress("Track_SigmaYPx", &vec_tracks.SigmaYPx);
    fTree_PackedEvents->SetBranchAddress("Track_SigmaZPx", &vec_tracks.SigmaZPx);
    fTree_PackedEvents->SetBranchAddress("Track_SigmaPx2", &vec_tracks.SigmaPx2);
    fTree_PackedEvents->SetBranchAddress("Track_SigmaXPy", &vec_tracks.SigmaXPy);
    fTree_PackedEvents->SetBranchAddress("Track_SigmaYPy", &vec_tracks.SigmaYPy);
    fTree_PackedEvents->SetBranchAddress("Track_SigmaZPy", &vec_tracks.SigmaZPy);
    fTree_PackedEvents->SetBranchAddress("Track_SigmaPxPy", &vec_tracks.SigmaPxPy);
    fTree_PackedEvents->SetBranchAddress("Track_SigmaPy2", &vec_tracks.SigmaPy2);
    fTree_PackedEvents->SetBranchAddress("Track_SigmaXPz", &vec_tracks.SigmaXPz);
    fTree_PackedEvents->SetBranchAddress("Track_SigmaYPz", &vec_tracks.SigmaYPz);
    fTree_PackedEvents->SetBranchAddress("Track_SigmaZPz", &vec_tracks.SigmaZPz);
    fTree_PackedEvents->SetBranchAddress("Track_SigmaPxPz", &vec_tracks.SigmaPxPz);
    fTree_PackedEvents->SetBranchAddress("Track_SigmaPyPz", &vec_tracks.SigmaPyPz);
    fTree_PackedEvents->SetBranchAddress("Track_SigmaPz2", &vec_tracks.SigmaPz2);
    fTree_PackedEvents->SetBranchAddress("Track_SigmaXE", &vec_tracks.SigmaXE);
    fTree_PackedEvents->SetBranchAddress("Track_SigmaYE", &vec_tracks.SigmaYE);
    fTree_PackedEvents->SetBranchAddress("Track_SigmaZE", &vec_tracks.SigmaZE);
    fTree_PackedEvents->SetBranchAddress("Track_SigmaPxE", &vec_tracks.SigmaPxE);
    fTree_PackedEvents->SetBranchAddress("Track_SigmaPyE", &vec_tracks.SigmaPyE);
    fTree_PackedEvents->SetBranchAddress("Track_SigmaPzE", &vec_tracks.SigmaPzE);
    fTree_PackedEvents->SetBranchAddress("Track_SigmaE2", &vec_tracks.SigmaE2);

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
    // PENDING
}
void Finder::CreateOutputBranches(Found::ChannelD& out_branches) {
    // PENDING
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
    for (size_t i{0}; i < in->Entry->size(); ++i) {
        KF::V0 kf_v0;
        // PENDING
        kf_v0.idx = in->Entry->at(i);
        kf_v0.idx_neg = in->Neg_Entry->at(i);
        kf_v0.idx_pos = in->Pos_Entry->at(i);
        kf_v0.pdg_code_hyp = static_cast<int>(pdg_code);
        vec_kf->push_back(kf_v0);
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
    switch (pdg_code) {
        case PdgCode::NegKaon:
            in = &fPacked_NegKaons;
            vec_kf = &fKF_NegKaons;
            break;
        case PdgCode::PosKaon:
            in = &fPacked_PosKaons;
            vec_kf = &fKF_PosKaons;
            break;
        case PdgCode::PiMinus:
            in = &fPacked_PiMinus;
            vec_kf = &fKF_PiMinus;
            break;
        case PdgCode::PiPlus:
            in = &fPacked_PiPlus;
            vec_kf = &fKF_PiPlus;
            break;
        default:
            return;
    }

    // loop over tracks //
    for (size_t i{0}; i < in->Entry->size(); ++i) {
        KF::Track kf_track;
        // PENDING
        vec_kf->push_back(kf_track);
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
    int n_sexa{0};
    for (const auto& v0a : *KF_Lambdas) {
        for (const auto& v0b : fKF_KaonsZeroShort) {
            // sanity check //
            std::set<size_t> unique_track_entries{v0a.idx_neg, v0a.idx_pos, v0b.idx_neg, v0b.idx_pos};
            if (unique_track_entries.size() < 4) continue;
            // fit //
            KF::ChannelA sexa;
            // PENDING
            sexa.AddDaughter(v0a, fInput_Event.MagneticField);
            sexa.AddDaughter(v0b, fInput_Event.MagneticField);
            // apply cuts //
            if (!PassesCuts(sexa)) continue;
            ++n_sexa;
            // store //
            Store(sexa);
        }
    }

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

bool Finder::PassesCuts(const KF::ChannelA& sexa) const {
    // PENDING
    return true;
}

void Finder::Store(const KF::ChannelA& sexa) {
    // PENDING
}

// ## Channel D ZONE ## //

void Finder::FindSexaquarks_ChannelD(bool anti_channel) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    const std::vector<KF::V0>* KF_Lambdas{anti_channel ? &fKF_Lambdas : &fKF_AntiLambdas};
    const std::vector<KF::Track>* KF_Kaons{anti_channel ? &fKF_NegKaons : &fKF_PosKaons};

    // loop over all possible pairs of (anti)lambda + (pos/neg)kaon //
    int n_sexa{0};
    for (const auto& v0 : *KF_Lambdas) {
        for (const auto& bach : *KF_Kaons) {
            // sanity check //
            std::set<size_t> unique_track_entries{v0.idx_neg, v0.idx_pos, bach.idx};
            if (unique_track_entries.size() < 3) continue;
            // fit //
            KF::ChannelD sexa;
            // PENDING
            sexa.AddDaughter(v0, fInput_Event.MagneticField);
            sexa.AddDaughter(bach, fInput_Event.MagneticField);
            // apply cuts //
            if (!PassesCuts(sexa)) continue;
            ++n_sexa;
            // store //
            Store(sexa);
        }
    }

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

bool Finder::PassesCuts(const KF::ChannelD& sexa) const {
    //
    return true;
}

void Finder::Store(const KF::ChannelD& sexa) {}

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

    // fill tree //
    fOutputTree->Fill();

    // clear temporary containers //
    fPacked_AntiLambdas.Clear(IsMC());
    fPacked_Lambdas.Clear(IsMC());
    fPacked_KaonsZeroShort.Clear(IsMC());

    fPacked_NegKaons.Clear(IsMC());
    fPacked_PosKaons.Clear(IsMC());
    fPacked_PiMinus.Clear(IsMC());
    fPacked_PiPlus.Clear(IsMC());

    // clear output branches //
    switch (GetReactionChannel()) {
        // standard channels //
        case ReactionChannel::A:
        case ReactionChannel::AntiA:
            fOutput_ChannelA.Clear(IsMC());
            break;
        case ReactionChannel::D:
        case ReactionChannel::AntiD:
            fOutput_ChannelD.Clear(IsMC());
            break;
        case ReactionChannel::E:
        case ReactionChannel::AntiE:
            // PENDING
            break;
        case ReactionChannel::H:
        case ReactionChannel::AntiH:
            // PENDING
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
