#include <iostream>
#include <memory>

#include "Math/Constants.hxx"
#include "Math/KFWrapper.hxx"
#include "Packager/Cuts.hxx"
#include "Packager/Packager.hxx"

namespace Tree2Secondaries {

bool Packager::Initialize() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    fEventsTree = std::make_unique<TChain>("Events");
    for (const auto& path : fSettings.PathInputFiles) {
        if (fEventsTree->Add(path.c_str()) == 0) {
            std::cerr << "Couldn't add TFile \"" << path << "\"" << '\n';
        }
    }
    if (!fEventsTree->GetEntries()) {
        std::cerr << "Couldn't manage to read any entry." << '\n';
        return false;
    }
    std::cout << "TChain \"" << fEventsTree->GetName() << "\" loaded successfully with " << fEventsTree->GetNtrees() << " trees and "
              << fEventsTree->GetEntries() << " total entries." << '\n';

    ConnectInputBranches();

    if (!PrepareOutputFile()) return false;
    if (!PrepareOutputTree()) return false;
    CreateOutputBranches();

    std::cout << "Packager initialized successfully" << '\n';

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
    return true;
}

// ## INPUT ZONE ## //

void Packager::ConnectInputBranches() {
    fEventsTree->SetBranchStatus("*", false);

    ConnectBranchesEvents();
    if (IsMC()) {
        ConnectBranchesMC();
        if (IsSignalMC()) ConnectBranchesInjected();
    }
    ConnectBranchesTracks();
}

void Packager::ConnectBranchesEvents() {
    fEventsTree->SetBranchStatus("RunNumber", true);
    fEventsTree->SetBranchStatus("DirNumber", true);
    if (!IsMC()) fEventsTree->SetBranchStatus("DirNumberB", true);
    fEventsTree->SetBranchStatus("EventNumber", true);
    fEventsTree->SetBranchStatus("Centrality", true);
    fEventsTree->SetBranchStatus("MagneticField", true);
    fEventsTree->SetBranchStatus("PV_Xv", true);
    fEventsTree->SetBranchStatus("PV_Yv", true);
    fEventsTree->SetBranchStatus("PV_Zv", true);

    fEventsTree->SetBranchAddress("RunNumber", &fInput_Event.RunNumber);
    fEventsTree->SetBranchAddress("DirNumber", &fInput_Event.DirNumber);
    if (!IsMC()) fEventsTree->SetBranchAddress("DirNumberB", &fInput_Event.DirNumberB);
    fEventsTree->SetBranchAddress("EventNumber", &fInput_Event.EventNumber);
    fEventsTree->SetBranchAddress("Centrality", &fInput_Event.Centrality);
    fEventsTree->SetBranchAddress("MagneticField", &fInput_Event.MagneticField);
    fEventsTree->SetBranchAddress("PV_Xv", &fInput_Event.PV_Xv);
    fEventsTree->SetBranchAddress("PV_Yv", &fInput_Event.PV_Yv);
    fEventsTree->SetBranchAddress("PV_Zv", &fInput_Event.PV_Zv);

    if (IsMC()) {
        fEventsTree->SetBranchStatus("MC_PV_Xv", true);
        fEventsTree->SetBranchStatus("MC_PV_Yv", true);
        fEventsTree->SetBranchStatus("MC_PV_Zv", true);

        fEventsTree->SetBranchAddress("MC_PV_Xv", &fInput_Event.MC_PV_Xv);
        fEventsTree->SetBranchAddress("MC_PV_Yv", &fInput_Event.MC_PV_Yv);
        fEventsTree->SetBranchAddress("MC_PV_Zv", &fInput_Event.MC_PV_Zv);
    }
}

void Packager::ConnectBranchesInjected() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif
    fEventsTree->SetBranchStatus("ReactionID", true);
    fEventsTree->SetBranchStatus("Sexaquark_Px", true);
    fEventsTree->SetBranchStatus("Sexaquark_Py", true);
    fEventsTree->SetBranchStatus("Sexaquark_Pz", true);
    fEventsTree->SetBranchStatus("Nucleon_Px", true);
    fEventsTree->SetBranchStatus("Nucleon_Py", true);
    fEventsTree->SetBranchStatus("Nucleon_Pz", true);

    fEventsTree->SetBranchAddress("ReactionID", &fInput_Injected.ReactionID);
    fEventsTree->SetBranchAddress("Sexaquark_Px", &fInput_Injected.Px);
    fEventsTree->SetBranchAddress("Sexaquark_Py", &fInput_Injected.Py);
    fEventsTree->SetBranchAddress("Sexaquark_Pz", &fInput_Injected.Pz);
    fEventsTree->SetBranchAddress("Nucleon_Px", &fInput_Injected.Nucleon_Px);
    fEventsTree->SetBranchAddress("Nucleon_Py", &fInput_Injected.Nucleon_Py);
    fEventsTree->SetBranchAddress("Nucleon_Pz", &fInput_Injected.Nucleon_Pz);
#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Packager::ConnectBranchesMC() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif
    fEventsTree->SetBranchStatus("MC_X", true);
    fEventsTree->SetBranchStatus("MC_Y", true);
    fEventsTree->SetBranchStatus("MC_Z", true);
    fEventsTree->SetBranchStatus("MC_Px", true);
    fEventsTree->SetBranchStatus("MC_Py", true);
    fEventsTree->SetBranchStatus("MC_Pz", true);
    fEventsTree->SetBranchStatus("MC_E", true);

    fEventsTree->SetBranchStatus("MC_PdgCode", true);
    fEventsTree->SetBranchStatus("MC_Mother_McEntry", true);
    fEventsTree->SetBranchStatus("MC_Status", true);
    fEventsTree->SetBranchStatus("MC_Generator", true);
    fEventsTree->SetBranchStatus("MC_IsPrimary", true);
    fEventsTree->SetBranchStatus("MC_IsSecFromMat", true);
    fEventsTree->SetBranchStatus("MC_IsSecFromWeak", true);

    fEventsTree->SetBranchAddress("MC_X", &fTruthHandler.fInput_MC.X);
    fEventsTree->SetBranchAddress("MC_Y", &fTruthHandler.fInput_MC.Y);
    fEventsTree->SetBranchAddress("MC_Z", &fTruthHandler.fInput_MC.Z);
    fEventsTree->SetBranchAddress("MC_Px", &fTruthHandler.fInput_MC.Px);
    fEventsTree->SetBranchAddress("MC_Py", &fTruthHandler.fInput_MC.Py);
    fEventsTree->SetBranchAddress("MC_Pz", &fTruthHandler.fInput_MC.Pz);
    fEventsTree->SetBranchAddress("MC_E", &fTruthHandler.fInput_MC.E);

    fEventsTree->SetBranchAddress("MC_PdgCode", &fTruthHandler.fInput_MC.PdgCode);
    fEventsTree->SetBranchAddress("MC_Mother_McEntry", &fTruthHandler.fInput_MC.MotherEntry);
    fEventsTree->SetBranchAddress("MC_Status", &fTruthHandler.fInput_MC.Status);
    fEventsTree->SetBranchAddress("MC_Generator", &fTruthHandler.fInput_MC.Generator);
    fEventsTree->SetBranchAddress("MC_IsPrimary", &fTruthHandler.fInput_MC.IsPrimary);
    fEventsTree->SetBranchAddress("MC_IsSecFromMat", &fTruthHandler.fInput_MC.IsSecFromMat);
    fEventsTree->SetBranchAddress("MC_IsSecFromWeak", &fTruthHandler.fInput_MC.IsSecFromWeak);
#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Packager::ConnectBranchesTracks() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif
    fEventsTree->SetBranchStatus("Track_X", true);
    fEventsTree->SetBranchStatus("Track_Y", true);
    fEventsTree->SetBranchStatus("Track_Z", true);
    fEventsTree->SetBranchStatus("Track_Px", true);
    fEventsTree->SetBranchStatus("Track_Py", true);
    fEventsTree->SetBranchStatus("Track_Pz", true);
    fEventsTree->SetBranchStatus("Track_Charge", true);
    fEventsTree->SetBranchStatus("Track_NSigmaPion", true);
    fEventsTree->SetBranchStatus("Track_NSigmaKaon", true);
    fEventsTree->SetBranchStatus("Track_NSigmaProton", true);

    fEventsTree->SetBranchStatus("Track_Alpha", true);
    fEventsTree->SetBranchStatus("Track_Snp", true);
    fEventsTree->SetBranchStatus("Track_Tgl", true);
    fEventsTree->SetBranchStatus("Track_Signed1Pt", true);
    fEventsTree->SetBranchStatus("Track_SigmaY2", true);
    fEventsTree->SetBranchStatus("Track_SigmaZY", true);
    fEventsTree->SetBranchStatus("Track_SigmaZ2", true);
    fEventsTree->SetBranchStatus("Track_SigmaSnpY", true);
    fEventsTree->SetBranchStatus("Track_SigmaSnpZ", true);
    fEventsTree->SetBranchStatus("Track_SigmaSnp2", true);
    fEventsTree->SetBranchStatus("Track_SigmaTglY", true);
    fEventsTree->SetBranchStatus("Track_SigmaTglZ", true);
    fEventsTree->SetBranchStatus("Track_SigmaTglSnp", true);
    fEventsTree->SetBranchStatus("Track_SigmaTgl2", true);
    fEventsTree->SetBranchStatus("Track_Sigma1PtY", true);
    fEventsTree->SetBranchStatus("Track_Sigma1PtZ", true);
    fEventsTree->SetBranchStatus("Track_Sigma1PtSnp", true);
    fEventsTree->SetBranchStatus("Track_Sigma1PtTgl", true);
    fEventsTree->SetBranchStatus("Track_Sigma1Pt2", true);

    fEventsTree->SetBranchAddress("Track_X", &fInput_Tracks.X);
    fEventsTree->SetBranchAddress("Track_Y", &fInput_Tracks.Y);
    fEventsTree->SetBranchAddress("Track_Z", &fInput_Tracks.Z);
    fEventsTree->SetBranchAddress("Track_Px", &fInput_Tracks.Px);
    fEventsTree->SetBranchAddress("Track_Py", &fInput_Tracks.Py);
    fEventsTree->SetBranchAddress("Track_Pz", &fInput_Tracks.Pz);
    fEventsTree->SetBranchAddress("Track_Charge", &fInput_Tracks.Charge);
    fEventsTree->SetBranchAddress("Track_NSigmaPion", &fInput_Tracks.NSigmaPion);
    fEventsTree->SetBranchAddress("Track_NSigmaKaon", &fInput_Tracks.NSigmaKaon);
    fEventsTree->SetBranchAddress("Track_NSigmaProton", &fInput_Tracks.NSigmaProton);

    fEventsTree->SetBranchAddress("Track_Alpha", &fInput_Tracks.Alpha);
    fEventsTree->SetBranchAddress("Track_Snp", &fInput_Tracks.Snp);
    fEventsTree->SetBranchAddress("Track_Tgl", &fInput_Tracks.Tgl);
    fEventsTree->SetBranchAddress("Track_Signed1Pt", &fInput_Tracks.Signed1Pt);
    fEventsTree->SetBranchAddress("Track_SigmaY2", &fInput_Tracks.SigmaY2);
    fEventsTree->SetBranchAddress("Track_SigmaZY", &fInput_Tracks.SigmaZY);
    fEventsTree->SetBranchAddress("Track_SigmaZ2", &fInput_Tracks.SigmaZ2);
    fEventsTree->SetBranchAddress("Track_SigmaSnpY", &fInput_Tracks.SigmaSnpY);
    fEventsTree->SetBranchAddress("Track_SigmaSnpZ", &fInput_Tracks.SigmaSnpZ);
    fEventsTree->SetBranchAddress("Track_SigmaSnp2", &fInput_Tracks.SigmaSnp2);
    fEventsTree->SetBranchAddress("Track_SigmaTglY", &fInput_Tracks.SigmaTglY);
    fEventsTree->SetBranchAddress("Track_SigmaTglZ", &fInput_Tracks.SigmaTglZ);
    fEventsTree->SetBranchAddress("Track_SigmaTglSnp", &fInput_Tracks.SigmaTglSnp);
    fEventsTree->SetBranchAddress("Track_SigmaTgl2", &fInput_Tracks.SigmaTgl2);
    fEventsTree->SetBranchAddress("Track_Sigma1PtY", &fInput_Tracks.Sigma1PtY);
    fEventsTree->SetBranchAddress("Track_Sigma1PtZ", &fInput_Tracks.Sigma1PtZ);
    fEventsTree->SetBranchAddress("Track_Sigma1PtSnp", &fInput_Tracks.Sigma1PtSnp);
    fEventsTree->SetBranchAddress("Track_Sigma1PtTgl", &fInput_Tracks.Sigma1PtTgl);
    fEventsTree->SetBranchAddress("Track_Sigma1Pt2", &fInput_Tracks.Sigma1Pt2);

    if (IsMC()) {
        fEventsTree->SetBranchStatus("Track_McEntry", true);

        fEventsTree->SetBranchAddress("Track_McEntry", &fInput_Tracks.McEntry);
    }
#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

// ## OUTPUT ZONE ## //

bool Packager::PrepareOutputFile() {
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

bool Packager::PrepareOutputTree() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    fOutputTree = std::make_unique<TTree>("PackedEvents", "Packed Events");
    if (!fOutputTree) {
        std::cerr << "TTree \"PackedEvents\" couldn't be created" << '\n';
        return false;
    }

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
    return true;
}

void Packager::CreateOutputBranches() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    CreateOutputBranchesEvents();
    if (IsSignalMC()) CreateOutputBranchesInjected();

    switch (GetReactionChannel()) {
        // standard channels //
        case ReactionChannel::A:
            CreateOutputBranchesV0s(static_cast<std::string>(Acronym::AntiLambda), fOutput_AntiLambdas);
            CreateOutputBranchesV0s(static_cast<std::string>(Acronym::KaonZeroShort), fOutput_KaonsZeroShort);
            break;
        case ReactionChannel::D:
            CreateOutputBranchesV0s(static_cast<std::string>(Acronym::AntiLambda), fOutput_AntiLambdas);
            CreateOutputBranchesTracks(static_cast<std::string>(Acronym::PosKaon), fOutput_PosKaons);
            break;
        case ReactionChannel::E:
            CreateOutputBranchesV0s(static_cast<std::string>(Acronym::AntiLambda), fOutput_AntiLambdas);
            CreateOutputBranchesTracks(static_cast<std::string>(Acronym::PosKaon), fOutput_PosKaons);
            CreateOutputBranchesTracks(static_cast<std::string>(Acronym::PiMinus), fOutput_PiMinus);
            CreateOutputBranchesTracks(static_cast<std::string>(Acronym::PiPlus), fOutput_PiPlus);
            break;
        case ReactionChannel::H:
            CreateOutputBranchesTracks(static_cast<std::string>(Acronym::PosKaon), fOutput_PosKaons);
            break;
        // anti-channels //
        case ReactionChannel::AntiA:
            CreateOutputBranchesV0s(static_cast<std::string>(Acronym::Lambda), fOutput_Lambdas);
            CreateOutputBranchesV0s(static_cast<std::string>(Acronym::KaonZeroShort), fOutput_KaonsZeroShort);
            break;
        case ReactionChannel::AntiD:
            CreateOutputBranchesV0s(static_cast<std::string>(Acronym::Lambda), fOutput_Lambdas);
            CreateOutputBranchesTracks(static_cast<std::string>(Acronym::NegKaon), fOutput_NegKaons);
            break;
        case ReactionChannel::AntiE:
            CreateOutputBranchesV0s(static_cast<std::string>(Acronym::Lambda), fOutput_Lambdas);
            CreateOutputBranchesTracks(static_cast<std::string>(Acronym::NegKaon), fOutput_NegKaons);
            CreateOutputBranchesTracks(static_cast<std::string>(Acronym::PiMinus), fOutput_PiMinus);
            CreateOutputBranchesTracks(static_cast<std::string>(Acronym::PiPlus), fOutput_PiPlus);
            break;
        case ReactionChannel::AntiH:
            CreateOutputBranchesTracks(static_cast<std::string>(Acronym::NegKaon), fOutput_NegKaons);
            break;
        // for data //
        case ReactionChannel::All:
            CreateOutputBranchesV0s(static_cast<std::string>(Acronym::AntiLambda), fOutput_AntiLambdas);
            CreateOutputBranchesV0s(static_cast<std::string>(Acronym::Lambda), fOutput_Lambdas);
            CreateOutputBranchesV0s(static_cast<std::string>(Acronym::KaonZeroShort), fOutput_KaonsZeroShort);
            CreateOutputBranchesTracks(static_cast<std::string>(Acronym::NegKaon), fOutput_NegKaons);
            CreateOutputBranchesTracks(static_cast<std::string>(Acronym::PosKaon), fOutput_PosKaons);
            CreateOutputBranchesTracks(static_cast<std::string>(Acronym::PiMinus), fOutput_PiMinus);
            CreateOutputBranchesTracks(static_cast<std::string>(Acronym::PiPlus), fOutput_PiPlus);
            break;
    }  // end of switch statement
#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Packager::CreateOutputBranchesEvents() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    fOutputTree->Branch("RunNumber", &fOutput_Event.RunNumber);
    fOutputTree->Branch("DirNumber", &fOutput_Event.DirNumber);
    if (!IsMC()) fOutputTree->Branch("DirNumberB", &fOutput_Event.DirNumberB);
    fOutputTree->Branch("EventNumber", &fOutput_Event.EventNumber);
    fOutputTree->Branch("Centrality", &fOutput_Event.Centrality);
    fOutputTree->Branch("MagneticField", &fOutput_Event.MagneticField);
    fOutputTree->Branch("PV_Xv", &fOutput_Event.PV_Xv);
    fOutputTree->Branch("PV_Yv", &fOutput_Event.PV_Yv);
    fOutputTree->Branch("PV_Zv", &fOutput_Event.PV_Zv);

    if (IsMC()) {
        fOutputTree->Branch("MC_PV_Xv", &fOutput_Event.MC_PV_Xv);
        fOutputTree->Branch("MC_PV_Yv", &fOutput_Event.MC_PV_Yv);
        fOutputTree->Branch("MC_PV_Zv", &fOutput_Event.MC_PV_Zv);
    }
#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Packager::CreateOutputBranchesInjected() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    fOutputTree->Branch("ReactionID", &fOutput_Injected.ReactionID);
    fOutputTree->Branch("Sexaquark_Px", &fOutput_Injected.Px);
    fOutputTree->Branch("Sexaquark_Py", &fOutput_Injected.Py);
    fOutputTree->Branch("Sexaquark_Pz", &fOutput_Injected.Pz);
    fOutputTree->Branch("Nucleon_Px", &fOutput_Injected.Nucleon_Px);
    fOutputTree->Branch("Nucleon_Py", &fOutput_Injected.Nucleon_Py);
    fOutputTree->Branch("Nucleon_Pz", &fOutput_Injected.Nucleon_Pz);
#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Packager::CreateOutputBranchesV0s(const std::string& name_v0, PackedEvents::V0s& out_branches) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    fOutputTree->Branch((name_v0 + "_Entry").c_str(), &out_branches.Entry);
    fOutputTree->Branch((name_v0 + "_X").c_str(), &out_branches.X);
    fOutputTree->Branch((name_v0 + "_Y").c_str(), &out_branches.Y);
    fOutputTree->Branch((name_v0 + "_Z").c_str(), &out_branches.Z);
    fOutputTree->Branch((name_v0 + "_Px").c_str(), &out_branches.Px);
    fOutputTree->Branch((name_v0 + "_Py").c_str(), &out_branches.Py);
    fOutputTree->Branch((name_v0 + "_Pz").c_str(), &out_branches.Pz);
    fOutputTree->Branch((name_v0 + "_E").c_str(), &out_branches.E);

    fOutputTree->Branch((name_v0 + "_SigmaX2").c_str(), &out_branches.Sigma.X2);
    fOutputTree->Branch((name_v0 + "_SigmaXY").c_str(), &out_branches.Sigma.XY);
    fOutputTree->Branch((name_v0 + "_SigmaY2").c_str(), &out_branches.Sigma.Y2);
    fOutputTree->Branch((name_v0 + "_SigmaXZ").c_str(), &out_branches.Sigma.XZ);
    fOutputTree->Branch((name_v0 + "_SigmaYZ").c_str(), &out_branches.Sigma.YZ);
    fOutputTree->Branch((name_v0 + "_SigmaZ2").c_str(), &out_branches.Sigma.Z2);
    fOutputTree->Branch((name_v0 + "_SigmaXPx").c_str(), &out_branches.Sigma.XPx);
    fOutputTree->Branch((name_v0 + "_SigmaYPx").c_str(), &out_branches.Sigma.YPx);
    fOutputTree->Branch((name_v0 + "_SigmaZPx").c_str(), &out_branches.Sigma.ZPx);
    fOutputTree->Branch((name_v0 + "_SigmaPx2").c_str(), &out_branches.Sigma.Px2);
    fOutputTree->Branch((name_v0 + "_SigmaXPy").c_str(), &out_branches.Sigma.XPy);
    fOutputTree->Branch((name_v0 + "_SigmaYPy").c_str(), &out_branches.Sigma.YPy);
    fOutputTree->Branch((name_v0 + "_SigmaZPy").c_str(), &out_branches.Sigma.ZPy);
    fOutputTree->Branch((name_v0 + "_SigmaPxPy").c_str(), &out_branches.Sigma.PxPy);
    fOutputTree->Branch((name_v0 + "_SigmaPy2").c_str(), &out_branches.Sigma.Py2);
    fOutputTree->Branch((name_v0 + "_SigmaXPz").c_str(), &out_branches.Sigma.XPz);
    fOutputTree->Branch((name_v0 + "_SigmaYPz").c_str(), &out_branches.Sigma.YPz);
    fOutputTree->Branch((name_v0 + "_SigmaZPz").c_str(), &out_branches.Sigma.ZPz);
    fOutputTree->Branch((name_v0 + "_SigmaPxPz").c_str(), &out_branches.Sigma.PxPz);
    fOutputTree->Branch((name_v0 + "_SigmaPyPz").c_str(), &out_branches.Sigma.PyPz);
    fOutputTree->Branch((name_v0 + "_SigmaPz2").c_str(), &out_branches.Sigma.Pz2);
    fOutputTree->Branch((name_v0 + "_SigmaXE").c_str(), &out_branches.Sigma.XE);
    fOutputTree->Branch((name_v0 + "_SigmaYE").c_str(), &out_branches.Sigma.YE);
    fOutputTree->Branch((name_v0 + "_SigmaZE").c_str(), &out_branches.Sigma.ZE);
    fOutputTree->Branch((name_v0 + "_SigmaPxE").c_str(), &out_branches.Sigma.PxE);
    fOutputTree->Branch((name_v0 + "_SigmaPyE").c_str(), &out_branches.Sigma.PyE);
    fOutputTree->Branch((name_v0 + "_SigmaPzE").c_str(), &out_branches.Sigma.PzE);
    fOutputTree->Branch((name_v0 + "_SigmaE2").c_str(), &out_branches.Sigma.E2);

    fOutputTree->Branch((name_v0 + "_Neg_Entry").c_str(), &out_branches.Neg.Entry);
    fOutputTree->Branch((name_v0 + "_Neg_X").c_str(), &out_branches.Neg.X);
    fOutputTree->Branch((name_v0 + "_Neg_Y").c_str(), &out_branches.Neg.Y);
    fOutputTree->Branch((name_v0 + "_Neg_Z").c_str(), &out_branches.Neg.Z);
    fOutputTree->Branch((name_v0 + "_Neg_Px").c_str(), &out_branches.Neg.Px);
    fOutputTree->Branch((name_v0 + "_Neg_Py").c_str(), &out_branches.Neg.Py);
    fOutputTree->Branch((name_v0 + "_Neg_Pz").c_str(), &out_branches.Neg.Pz);
    fOutputTree->Branch((name_v0 + "_Neg_E").c_str(), &out_branches.Neg.E);

    fOutputTree->Branch((name_v0 + "_Pos_Entry").c_str(), &out_branches.Pos.Entry);
    fOutputTree->Branch((name_v0 + "_Pos_X").c_str(), &out_branches.Pos.X);
    fOutputTree->Branch((name_v0 + "_Pos_Y").c_str(), &out_branches.Pos.Y);
    fOutputTree->Branch((name_v0 + "_Pos_Z").c_str(), &out_branches.Pos.Z);
    fOutputTree->Branch((name_v0 + "_Pos_Px").c_str(), &out_branches.Pos.Px);
    fOutputTree->Branch((name_v0 + "_Pos_Py").c_str(), &out_branches.Pos.Py);
    fOutputTree->Branch((name_v0 + "_Pos_Pz").c_str(), &out_branches.Pos.Pz);
    fOutputTree->Branch((name_v0 + "_Pos_E").c_str(), &out_branches.Pos.E);

    if (IsMC()) {
        fOutputTree->Branch((name_v0 + "_MC_Entry").c_str(), &out_branches.True.Entry);

        fOutputTree->Branch((name_v0 + "_MC_X").c_str(), &out_branches.True.X);
        fOutputTree->Branch((name_v0 + "_MC_Y").c_str(), &out_branches.True.Y);
        fOutputTree->Branch((name_v0 + "_MC_Z").c_str(), &out_branches.True.Z);
        fOutputTree->Branch((name_v0 + "_MC_Px").c_str(), &out_branches.True.Px);
        fOutputTree->Branch((name_v0 + "_MC_Py").c_str(), &out_branches.True.Py);
        fOutputTree->Branch((name_v0 + "_MC_Pz").c_str(), &out_branches.True.Pz);
        fOutputTree->Branch((name_v0 + "_MC_E").c_str(), &out_branches.True.E);

        fOutputTree->Branch((name_v0 + "_MC_PdgCode").c_str(), &out_branches.True.PdgCode);
        fOutputTree->Branch((name_v0 + "_MC_MotherEntry").c_str(), &out_branches.True.MotherEntry);
        fOutputTree->Branch((name_v0 + "_MC_IsSignal").c_str(), &out_branches.True.IsSignal);
        fOutputTree->Branch((name_v0 + "_MC_IsPrimary").c_str(), &out_branches.True.IsPrimary);
        fOutputTree->Branch((name_v0 + "_MC_IsSecFromMat").c_str(), &out_branches.True.IsSecFromMat);
        fOutputTree->Branch((name_v0 + "_MC_IsSecFromWeak").c_str(), &out_branches.True.IsSecFromWeak);
        fOutputTree->Branch((name_v0 + "_MC_ReactionID").c_str(), &out_branches.True.ReactionID);
        /*
                fOutputTree->Branch((name_v0 + "_MC_Neg_Entry").c_str(), &out_branches.Neg.True.Entry);

                fOutputTree->Branch((name_v0 + "_MC_Neg_X").c_str(), &out_branches.Neg.True.X);
                fOutputTree->Branch((name_v0 + "_MC_Neg_Y").c_str(), &out_branches.Neg.True.Y);
                fOutputTree->Branch((name_v0 + "_MC_Neg_Z").c_str(), &out_branches.Neg.True.Z);
                fOutputTree->Branch((name_v0 + "_MC_Neg_Px").c_str(), &out_branches.Neg.True.Px);
                fOutputTree->Branch((name_v0 + "_MC_Neg_Py").c_str(), &out_branches.Neg.True.Py);
                fOutputTree->Branch((name_v0 + "_MC_Neg_Pz").c_str(), &out_branches.Neg.True.Pz);
                fOutputTree->Branch((name_v0 + "_MC_Neg_E").c_str(), &out_branches.Neg.True.E);

                fOutputTree->Branch((name_v0 + "_MC_Neg_PdgCode").c_str(), &out_branches.Neg.True.PdgCode);
                fOutputTree->Branch((name_v0 + "_MC_Neg_MotherEntry").c_str(), &out_branches.Neg.True.MotherEntry);
                fOutputTree->Branch((name_v0 + "_MC_Neg_IsSignal").c_str(), &out_branches.Neg.True.IsSignal);
                fOutputTree->Branch((name_v0 + "_MC_Neg_IsPrimary").c_str(), &out_branches.Neg.True.IsPrimary);
                fOutputTree->Branch((name_v0 + "_MC_Neg_IsSecFromMat").c_str(), &out_branches.Neg.True.IsSecFromMat);
                fOutputTree->Branch((name_v0 + "_MC_Neg_IsSecFromWeak").c_str(), &out_branches.Neg.True.IsSecFromWeak);
                fOutputTree->Branch((name_v0 + "_MC_Neg_ReactionID").c_str(), &out_branches.Neg.True.ReactionID);

                fOutputTree->Branch((name_v0 + "_MC_Pos_Entry").c_str(), &out_branches.Pos.True.Entry);

                fOutputTree->Branch((name_v0 + "_MC_Pos_Xv").c_str(), &out_branches.Pos.True.Xv);
                fOutputTree->Branch((name_v0 + "_MC_Pos_Yv").c_str(), &out_branches.Pos.True.Yv);
                fOutputTree->Branch((name_v0 + "_MC_Pos_Zv").c_str(), &out_branches.Pos.True.Zv);
                fOutputTree->Branch((name_v0 + "_MC_Pos_Px").c_str(), &out_branches.Pos.True.Px);
                fOutputTree->Branch((name_v0 + "_MC_Pos_Py").c_str(), &out_branches.Pos.True.Py);
                fOutputTree->Branch((name_v0 + "_MC_Pos_Pz").c_str(), &out_branches.Pos.True.Pz);
                fOutputTree->Branch((name_v0 + "_MC_Pos_E").c_str(), &out_branches.Pos.True.E);

                fOutputTree->Branch((name_v0 + "_MC_Pos_PdgCode").c_str(), &out_branches.Pos.True.PdgCode);
                fOutputTree->Branch((name_v0 + "_MC_Pos_MotherEntry").c_str(), &out_branches.Pos.True.MotherEntry);
                fOutputTree->Branch((name_v0 + "_MC_Pos_IsSignal").c_str(), &out_branches.Pos.True.IsSignal);
                fOutputTree->Branch((name_v0 + "_MC_Pos_IsPrimary").c_str(), &out_branches.Pos.True.IsPrimary);
                fOutputTree->Branch((name_v0 + "_MC_Pos_IsSecFromMat").c_str(), &out_branches.Pos.True.IsSecFromMat);
                fOutputTree->Branch((name_v0 + "_MC_Pos_IsSecFromWeak").c_str(), &out_branches.Pos.True.IsSecFromWeak);
                fOutputTree->Branch((name_v0 + "_MC_Pos_ReactionID").c_str(), &out_branches.Pos.True.ReactionID);
                */
    }
#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Packager::CreateOutputBranchesTracks(const std::string& name_part, PackedEvents::Tracks& out_branches) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    fOutputTree->Branch((name_part + "_Entry").c_str(), &out_branches.Entry);

    fOutputTree->Branch((name_part + "_X").c_str(), &out_branches.X);
    fOutputTree->Branch((name_part + "_Y").c_str(), &out_branches.Y);
    fOutputTree->Branch((name_part + "_Z").c_str(), &out_branches.Z);
    fOutputTree->Branch((name_part + "_Px").c_str(), &out_branches.Px);
    fOutputTree->Branch((name_part + "_Py").c_str(), &out_branches.Py);
    fOutputTree->Branch((name_part + "_Pz").c_str(), &out_branches.Pz);
    fOutputTree->Branch((name_part + "_E").c_str(), &out_branches.E);

    fOutputTree->Branch((name_part + "_SigmaX2").c_str(), &out_branches.Sigma.X2);
    fOutputTree->Branch((name_part + "_SigmaXY").c_str(), &out_branches.Sigma.XY);
    fOutputTree->Branch((name_part + "_SigmaY2").c_str(), &out_branches.Sigma.Y2);
    fOutputTree->Branch((name_part + "_SigmaXZ").c_str(), &out_branches.Sigma.XZ);
    fOutputTree->Branch((name_part + "_SigmaYZ").c_str(), &out_branches.Sigma.YZ);
    fOutputTree->Branch((name_part + "_SigmaZ2").c_str(), &out_branches.Sigma.Z2);
    fOutputTree->Branch((name_part + "_SigmaXPx").c_str(), &out_branches.Sigma.XPx);
    fOutputTree->Branch((name_part + "_SigmaYPx").c_str(), &out_branches.Sigma.YPx);
    fOutputTree->Branch((name_part + "_SigmaZPx").c_str(), &out_branches.Sigma.ZPx);
    fOutputTree->Branch((name_part + "_SigmaPx2").c_str(), &out_branches.Sigma.Px2);
    fOutputTree->Branch((name_part + "_SigmaXPy").c_str(), &out_branches.Sigma.XPy);
    fOutputTree->Branch((name_part + "_SigmYaPy").c_str(), &out_branches.Sigma.YPy);
    fOutputTree->Branch((name_part + "_SigmaZPy").c_str(), &out_branches.Sigma.ZPy);
    fOutputTree->Branch((name_part + "_SigmaPxPy").c_str(), &out_branches.Sigma.PxPy);
    fOutputTree->Branch((name_part + "_SigmaPy2").c_str(), &out_branches.Sigma.Py2);
    fOutputTree->Branch((name_part + "_SigmaXPz").c_str(), &out_branches.Sigma.XPz);
    fOutputTree->Branch((name_part + "_SigmaYPz").c_str(), &out_branches.Sigma.YPz);
    fOutputTree->Branch((name_part + "_SigmaZPz").c_str(), &out_branches.Sigma.ZPz);
    fOutputTree->Branch((name_part + "_SigmaPxPz").c_str(), &out_branches.Sigma.PxPz);
    fOutputTree->Branch((name_part + "_SigmaPyPz").c_str(), &out_branches.Sigma.PyPz);
    fOutputTree->Branch((name_part + "_SigmaPz2").c_str(), &out_branches.Sigma.Pz2);
    fOutputTree->Branch((name_part + "_SigmaXE").c_str(), &out_branches.Sigma.XE);
    fOutputTree->Branch((name_part + "_SigmaYE").c_str(), &out_branches.Sigma.YE);
    fOutputTree->Branch((name_part + "_SigmaZE").c_str(), &out_branches.Sigma.ZE);
    fOutputTree->Branch((name_part + "_SigmaPxE").c_str(), &out_branches.Sigma.PxE);
    fOutputTree->Branch((name_part + "_SigmaPyE").c_str(), &out_branches.Sigma.PyE);
    fOutputTree->Branch((name_part + "_SigmaPzE").c_str(), &out_branches.Sigma.PzE);
    fOutputTree->Branch((name_part + "_SigmaE2").c_str(), &out_branches.Sigma.E2);

    if (IsMC()) {
        fOutputTree->Branch((name_part + "_MC_Entry").c_str(), &out_branches.True.Entry);

        fOutputTree->Branch((name_part + "_MC_X").c_str(), &out_branches.True.X);
        fOutputTree->Branch((name_part + "_MC_Y").c_str(), &out_branches.True.Y);
        fOutputTree->Branch((name_part + "_MC_Z").c_str(), &out_branches.True.Z);
        fOutputTree->Branch((name_part + "_MC_Px").c_str(), &out_branches.True.Px);
        fOutputTree->Branch((name_part + "_MC_Py").c_str(), &out_branches.True.Py);
        fOutputTree->Branch((name_part + "_MC_Pz").c_str(), &out_branches.True.Pz);
        fOutputTree->Branch((name_part + "_MC_E").c_str(), &out_branches.True.E);

        fOutputTree->Branch((name_part + "_MC_PdgCode").c_str(), &out_branches.True.PdgCode);
        fOutputTree->Branch((name_part + "_MC_MotherEntry").c_str(), &out_branches.True.MotherEntry);
        fOutputTree->Branch((name_part + "_MC_IsSignal").c_str(), &out_branches.True.IsSignal);
        fOutputTree->Branch((name_part + "_MC_IsPrimary").c_str(), &out_branches.True.IsPrimary);
        fOutputTree->Branch((name_part + "_MC_IsSecFromMat").c_str(), &out_branches.True.IsSecFromMat);
        fOutputTree->Branch((name_part + "_MC_IsSecFromWeak").c_str(), &out_branches.True.IsSecFromWeak);
        fOutputTree->Branch((name_part + "_MC_ReactionID").c_str(), &out_branches.True.ReactionID);
    }
#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Packager::ProcessEvent() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
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
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Packager::ProcessInjected() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif
    fOutput_Injected.ReactionID = fInput_Injected.ReactionID;
    fOutput_Injected.Px = fInput_Injected.Px;
    fOutput_Injected.Py = fInput_Injected.Py;
    fOutput_Injected.Pz = fInput_Injected.Pz;
    fOutput_Injected.Nucleon_Px = fInput_Injected.Nucleon_Px;
    fOutput_Injected.Nucleon_Py = fInput_Injected.Nucleon_Py;
    fOutput_Injected.Nucleon_Pz = fInput_Injected.Nucleon_Pz;
#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

// ## Tracks ZONE ## //

// Store tracks' ESD indices into vectors, according to their respective track PID and charge.
void Packager::ProcessTracks() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif
    const auto n_tracks = NumberTracks();
    if (IsMC()) fTruthHandler.InitMap(n_tracks);
    for (size_t esd_track{0}; esd_track < n_tracks; ++esd_track) {
        // get charge //
        auto charge{fInput_Tracks.Charge->at(esd_track)};
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
        // MC //
        if (IsMC()) fTruthHandler.Link(esd_track, fInput_Tracks.McEntry->at(esd_track));  // PENDING: should i add protection?
    }
#ifdef T2S_DEBUG
    std::cout << "   " << __FUNCTION__ << " :: n_antiprotons = " << fVec_AntiProtons.size() << '\n';
    std::cout << "   " << __FUNCTION__ << " :: n_protons     = " << fVec_Protons.size() << '\n';
    std::cout << "   " << __FUNCTION__ << " :: n_negkaons    = " << fVec_NegKaons.size() << '\n';
    std::cout << "   " << __FUNCTION__ << " :: n_poskaons    = " << fVec_PosKaons.size() << '\n';
    std::cout << "   " << __FUNCTION__ << " :: n_piminus     = " << fVec_PiMinus.size() << '\n';
    std::cout << "   " << __FUNCTION__ << " :: n_piplus      = " << fVec_PiPlus.size() << '\n';
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

// Note: intended for light particles only, i.e., kaons and pions.
void Packager::PackTracks(PdgCode pdg_code) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif
    // determine rules based on particle pdg code //
    const std::vector<size_t>* vec;
    PackedEvents::Tracks* out;
    double mass;
    switch (pdg_code) {
        case PdgCode::NegKaon:
            vec = &fVec_NegKaons;
            out = &fOutput_NegKaons;
            mass = PdgMass::Kaon;
            break;
        case PdgCode::PosKaon:
            vec = &fVec_PosKaons;
            out = &fOutput_PosKaons;
            mass = PdgMass::Kaon;
            break;
        case PdgCode::PiMinus:
            vec = &fVec_PiMinus;
            out = &fOutput_PiMinus;
            mass = PdgMass::Pion;
            break;
        case PdgCode::PiPlus:
            vec = &fVec_PiPlus;
            out = &fOutput_PiPlus;
            mass = PdgMass::Pion;
            break;
        default:
            return;
    }

    // loop over selected tracks //
    for (auto esd_idx : *vec) {

        // prepare kf object //
        std::array<double, 6> neg_kf_params = PackParams_KF(esd_idx);
        std::array<float, 5> neg_alice_params = PackParams_ALICE(esd_idx);
        std::array<float, 15> neg_alice_cov = PackCovMatrix_ALICE(esd_idx);
        auto kf = KF::CreateParticle(neg_kf_params, neg_alice_params, neg_alice_cov, fInput_Tracks.Alpha->at(esd_idx),
                                     fInput_Tracks.Charge->at(esd_idx), mass);

        // store //
        KF::Track kf_track;
        static_cast<KF::Particle&>(kf_track) = kf;
        Store(kf_track, *out);
    }  // end of loop over selected tracks
#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Packager::Store(const KF::Track& track, PackedEvents::Tracks& out_branches) {

    out_branches.Entry->push_back(track.idx);
    out_branches.X->push_back(static_cast<float>(track.GetParameter(0)));
    out_branches.Y->push_back(static_cast<float>(track.GetParameter(1)));
    out_branches.Z->push_back(static_cast<float>(track.GetParameter(2)));
    out_branches.Px->push_back(static_cast<float>(track.GetParameter(3)));
    out_branches.Py->push_back(static_cast<float>(track.GetParameter(4)));
    out_branches.Pz->push_back(static_cast<float>(track.GetParameter(5)));
    out_branches.E->push_back(static_cast<float>(track.GetParameter(6)));

    out_branches.Sigma.X2->push_back(static_cast<float>(track.GetCovariance(0)));
    out_branches.Sigma.XY->push_back(static_cast<float>(track.GetCovariance(1)));
    out_branches.Sigma.Y2->push_back(static_cast<float>(track.GetCovariance(2)));
    out_branches.Sigma.XZ->push_back(static_cast<float>(track.GetCovariance(3)));
    out_branches.Sigma.YZ->push_back(static_cast<float>(track.GetCovariance(4)));
    out_branches.Sigma.Z2->push_back(static_cast<float>(track.GetCovariance(5)));
    out_branches.Sigma.XPx->push_back(static_cast<float>(track.GetCovariance(6)));
    out_branches.Sigma.YPx->push_back(static_cast<float>(track.GetCovariance(7)));
    out_branches.Sigma.ZPx->push_back(static_cast<float>(track.GetCovariance(8)));
    out_branches.Sigma.Px2->push_back(static_cast<float>(track.GetCovariance(9)));
    out_branches.Sigma.XPy->push_back(static_cast<float>(track.GetCovariance(10)));
    out_branches.Sigma.YPy->push_back(static_cast<float>(track.GetCovariance(11)));
    out_branches.Sigma.ZPy->push_back(static_cast<float>(track.GetCovariance(12)));
    out_branches.Sigma.PxPy->push_back(static_cast<float>(track.GetCovariance(13)));
    out_branches.Sigma.Py2->push_back(static_cast<float>(track.GetCovariance(14)));
    out_branches.Sigma.XPz->push_back(static_cast<float>(track.GetCovariance(15)));
    out_branches.Sigma.YPz->push_back(static_cast<float>(track.GetCovariance(16)));
    out_branches.Sigma.ZPz->push_back(static_cast<float>(track.GetCovariance(17)));
    out_branches.Sigma.PxPz->push_back(static_cast<float>(track.GetCovariance(18)));
    out_branches.Sigma.PyPz->push_back(static_cast<float>(track.GetCovariance(19)));
    out_branches.Sigma.Pz2->push_back(static_cast<float>(track.GetCovariance(20)));
    out_branches.Sigma.XE->push_back(static_cast<float>(track.GetCovariance(21)));
    out_branches.Sigma.YE->push_back(static_cast<float>(track.GetCovariance(22)));
    out_branches.Sigma.ZE->push_back(static_cast<float>(track.GetCovariance(23)));
    out_branches.Sigma.PxE->push_back(static_cast<float>(track.GetCovariance(24)));
    out_branches.Sigma.PyE->push_back(static_cast<float>(track.GetCovariance(25)));
    out_branches.Sigma.PzE->push_back(static_cast<float>(track.GetCovariance(26)));
    out_branches.Sigma.E2->push_back(static_cast<float>(track.GetCovariance(27)));

    if (IsMC()) {
        const int mc_entry{fTruthHandler.McEntry(track.idx)};
        out_branches.True.Entry->push_back(mc_entry);
        if (mc_entry == Const::DummyInt) return;
        out_branches.True.X->push_back(fTruthHandler.X(mc_entry));
        out_branches.True.Y->push_back(fTruthHandler.Y(mc_entry));
        out_branches.True.Z->push_back(fTruthHandler.Z(mc_entry));
        out_branches.True.Px->push_back(fTruthHandler.Px(mc_entry));
        out_branches.True.Py->push_back(fTruthHandler.Py(mc_entry));
        out_branches.True.Pz->push_back(fTruthHandler.Pz(mc_entry));
        out_branches.True.E->push_back(fTruthHandler.E(mc_entry));

        out_branches.True.PdgCode->push_back(fTruthHandler.PdgCode(mc_entry));
        out_branches.True.MotherEntry->push_back(fTruthHandler.MotherEntry(mc_entry));
        out_branches.True.IsSignal->push_back(fTruthHandler.IsSignal(mc_entry));
        out_branches.True.IsPrimary->push_back(fTruthHandler.IsPrimary(mc_entry));
        out_branches.True.IsSecFromMat->push_back(fTruthHandler.IsSecFromMat(mc_entry));
        out_branches.True.IsSecFromWeak->push_back(fTruthHandler.IsSecFromWeak(mc_entry));
        out_branches.True.ReactionID->push_back(fTruthHandler.ReactionID(mc_entry));
    }
}

// ## V0s ZONE ## //

void Packager::FindV0s(PdgCode pdg_code_v0) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    // determine rules based on V0 pdg code //
    const std::vector<size_t>* vec_neg;
    const std::vector<size_t>* vec_pos;
    PackedEvents::V0s* out;
    double mass_neg;
    double mass_pos;
    double mass_v0;
    switch (pdg_code_v0) {
        case PdgCode::AntiLambda:
            vec_neg = &fVec_AntiProtons;
            vec_pos = &fVec_PiPlus;
            out = &fOutput_AntiLambdas;
            mass_neg = PdgMass::Proton;
            mass_pos = PdgMass::Pion;
            mass_v0 = PdgMass::Lambda;
            break;
        case PdgCode::Lambda:
            vec_neg = &fVec_PiMinus;
            vec_pos = &fVec_Protons;
            out = &fOutput_Lambdas;
            mass_neg = PdgMass::Pion;
            mass_pos = PdgMass::Proton;
            mass_v0 = PdgMass::Lambda;
            break;
        case PdgCode::KaonZeroShort:
            vec_neg = &fVec_PiMinus;
            vec_pos = &fVec_PiPlus;
            out = &fOutput_KaonsZeroShort;
            mass_neg = PdgMass::Pion;
            mass_pos = PdgMass::Pion;
            mass_v0 = PdgMass::KaonZeroShort;
            break;
        default:
            return;
    }

    // loop over all possible pairs of tracks //
    size_t v0_entry{0};
    for (auto esd_neg : *vec_neg) {
        for (auto esd_pos : *vec_pos) {

            // sanity check //
            if (esd_neg == esd_pos) continue;

            // prepare neg //
            std::array<double, 6> neg_kf_params = PackParams_KF(esd_neg);
            std::array<float, 5> neg_alice_params = PackParams_ALICE(esd_neg);
            std::array<float, 15> neg_alice_cov = PackCovMatrix_ALICE(esd_neg);
            auto neg = KF::CreateParticle(neg_kf_params, neg_alice_params, neg_alice_cov, fInput_Tracks.Alpha->at(esd_neg),
                                          fInput_Tracks.Charge->at(esd_neg), mass_neg);

            // prepare pos //
            std::array<double, 6> pos_kf_params = PackParams_KF(esd_pos);
            std::array<float, 5> pos_alice_params = PackParams_ALICE(esd_pos);
            std::array<float, 15> pos_alice_cov = PackCovMatrix_ALICE(esd_pos);
            auto pos = KF::CreateParticle(pos_kf_params, pos_alice_params, pos_alice_cov, fInput_Tracks.Alpha->at(esd_pos),
                                          fInput_Tracks.Charge->at(esd_pos), mass_pos);

            // build v0 //
            KF::V0 v0;
            v0.AddDaughter(neg, fInput_Event.MagneticField);
            v0.AddDaughter(pos, fInput_Event.MagneticField);
            v0.AddMassConstraint(mass_v0);

            // apply cuts //
            if (!PassesCuts(v0, pdg_code_v0)) continue;

            // add info //
            v0.SetIndices({v0_entry, esd_neg, esd_pos});
            v0.Neg_Energy = neg.E();
            v0.Pos_Energy = pos.E();
            v0.pdg_code_hyp = static_cast<int>(pdg_code_v0);
#ifdef T2S_DEBUG
            std::cout << "   " << __FUNCTION__ << " :: idx,neg,pos=" << v0.idx << "," << v0.idx_neg << "," << v0.idx_pos;
            std::cout << ";x,y,z=" << v0.X() << "," << v0.Y() << "," << v0.Z();
            std::cout << ";x,y,z(neg)=" << v0.PCA_Neg()[0] << "," << v0.PCA_Neg()[1] << "," << v0.PCA_Neg()[2];
            std::cout << ";x,y,z(pos)=" << v0.PCA_Pos()[0] << "," << v0.PCA_Pos()[1] << "," << v0.PCA_Pos()[2];
            std::cout << ";mass=" << v0.Mass();
            std::cout << ";dca_dau=" << v0.DCA_Daughters();
            std::cout << ";radius=" << v0.Radius();
            std::cout << ";dca_neg=" << v0.DCA_Neg_V0();
            std::cout << ";dca_pos=" << v0.DCA_Pos_V0();
            std::cout << ";pt=" << v0.Pt();
            std::cout << ";eta=" << v0.Eta();
            std::cout << ";qt=" << v0.ArmenterosQt();
            std::cout << ";alpha=" << v0.ArmenterosAlpha();
            std::cout << ";cpa_pv=" << v0.CPA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv);
            std::cout << ";dca_pv=" << v0.DCA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv) << '\n';
#endif

            // store //
            Store(v0, *out);
            ++v0_entry;
        }  // end of loop over pos
    }  // end of loop over neg
#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

bool Packager::PassesCuts_Lambda(const KF::V0& v0) const {

    if (v0.Mass() < Cuts::Lambda::Min_Mass || v0.Mass() > Cuts::Lambda::Max_Mass) return false;
    if (v0.DCA_Daughters() > Cuts::Lambda::Max_DCAbtwDau) return false;
    if (v0.AbsZ() > Cuts::Lambda::AbsMax_Zv) return false;
    if (v0.Radius() < Cuts::Lambda::Min_Radius || v0.Radius() > Cuts::Lambda::Max_Radius) return false;
    if (v0.DCA_Neg_V0() > Cuts::Lambda::Max_DCAnegV0) return false;
    if (v0.DCA_Pos_V0() > Cuts::Lambda::Max_DCAposV0) return false;
    if (v0.Pt() < Cuts::Lambda::Min_Pt) return false;
    if (v0.AbsEta() > Cuts::Lambda::AbsMax_Eta) return false;
    if (v0.AbsArmQtOverAlpha() > Cuts::Lambda::AbsMax_ArmQtOverAlpha) return false;
    if (v0.CPA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv) < Cuts::Lambda::Min_CPAwrtPV ||
        v0.CPA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv) > Cuts::Lambda::Max_CPAwrtPV) {
        return false;
    }
    if (v0.DCA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv) < Cuts::Lambda::Min_DCAwrtPV) return false;

    return true;
}

bool Packager::PassesCuts_KaonZeroShort(const KF::V0& v0) const {

    if (v0.DCA_Daughters() > Cuts::KaonZeroShort::Max_DCAbtwDau) return false;
    if (v0.Pt() < Cuts::KaonZeroShort::Min_Pt) return false;
    if (v0.Mass() < Cuts::KaonZeroShort::Min_Mass || v0.Mass() > Cuts::KaonZeroShort::Max_Mass) return false;
    if (v0.AbsEta() > Cuts::KaonZeroShort::AbsMax_Eta) return false;
    if (v0.AbsZ() > Cuts::KaonZeroShort::AbsMax_Zv) return false;
    if (v0.Radius() < Cuts::KaonZeroShort::Min_Radius || v0.Radius() > Cuts::KaonZeroShort::Max_Radius) return false;
    if (v0.DCA_Neg_V0() > Cuts::KaonZeroShort::Max_DCAnegV0) return false;
    if (v0.DCA_Pos_V0() > Cuts::KaonZeroShort::Max_DCAposV0) return false;
    if (v0.CPA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv) < Cuts::KaonZeroShort::Min_CPAwrtPV ||
        v0.CPA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv) > Cuts::KaonZeroShort::Max_CPAwrtPV) {
        return false;
    }
    if (v0.DCA_Point(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv) < Cuts::KaonZeroShort::Min_DCAwrtPV) return false;

    return true;
}

void Packager::Store(const KF::V0& v0, PackedEvents::V0s& out_branches) {

    out_branches.Entry->push_back(v0.idx);
    out_branches.X->push_back(static_cast<float>(v0.X()));
    out_branches.Y->push_back(static_cast<float>(v0.Y()));
    out_branches.Z->push_back(static_cast<float>(v0.Z()));
    out_branches.Px->push_back(static_cast<float>(v0.Px()));
    out_branches.Py->push_back(static_cast<float>(v0.Py()));
    out_branches.Pz->push_back(static_cast<float>(v0.Pz()));
    out_branches.E->push_back(static_cast<float>(v0.E()));

    out_branches.Sigma.X2->push_back(static_cast<float>(v0.GetCovariance(0)));
    out_branches.Sigma.XY->push_back(static_cast<float>(v0.GetCovariance(1)));
    out_branches.Sigma.Y2->push_back(static_cast<float>(v0.GetCovariance(2)));
    out_branches.Sigma.XZ->push_back(static_cast<float>(v0.GetCovariance(3)));
    out_branches.Sigma.YZ->push_back(static_cast<float>(v0.GetCovariance(4)));
    out_branches.Sigma.Z2->push_back(static_cast<float>(v0.GetCovariance(5)));
    out_branches.Sigma.XPx->push_back(static_cast<float>(v0.GetCovariance(6)));
    out_branches.Sigma.YPx->push_back(static_cast<float>(v0.GetCovariance(7)));
    out_branches.Sigma.ZPx->push_back(static_cast<float>(v0.GetCovariance(8)));
    out_branches.Sigma.Px2->push_back(static_cast<float>(v0.GetCovariance(9)));
    out_branches.Sigma.XPy->push_back(static_cast<float>(v0.GetCovariance(10)));
    out_branches.Sigma.YPy->push_back(static_cast<float>(v0.GetCovariance(11)));
    out_branches.Sigma.ZPy->push_back(static_cast<float>(v0.GetCovariance(12)));
    out_branches.Sigma.PxPy->push_back(static_cast<float>(v0.GetCovariance(13)));
    out_branches.Sigma.Py2->push_back(static_cast<float>(v0.GetCovariance(14)));
    out_branches.Sigma.XPz->push_back(static_cast<float>(v0.GetCovariance(15)));
    out_branches.Sigma.YPz->push_back(static_cast<float>(v0.GetCovariance(16)));
    out_branches.Sigma.ZPz->push_back(static_cast<float>(v0.GetCovariance(17)));
    out_branches.Sigma.PxPz->push_back(static_cast<float>(v0.GetCovariance(18)));
    out_branches.Sigma.PyPz->push_back(static_cast<float>(v0.GetCovariance(19)));
    out_branches.Sigma.Pz2->push_back(static_cast<float>(v0.GetCovariance(20)));
    out_branches.Sigma.XE->push_back(static_cast<float>(v0.GetCovariance(21)));
    out_branches.Sigma.YE->push_back(static_cast<float>(v0.GetCovariance(22)));
    out_branches.Sigma.ZE->push_back(static_cast<float>(v0.GetCovariance(23)));
    out_branches.Sigma.PxE->push_back(static_cast<float>(v0.GetCovariance(24)));
    out_branches.Sigma.PyE->push_back(static_cast<float>(v0.GetCovariance(25)));
    out_branches.Sigma.PzE->push_back(static_cast<float>(v0.GetCovariance(26)));
    out_branches.Sigma.E2->push_back(static_cast<float>(v0.GetCovariance(27)));

    out_branches.Neg.Entry->push_back(v0.idx_neg);
    out_branches.Neg.X->push_back(static_cast<float>(v0.PCA_Neg()[0]));
    out_branches.Neg.Y->push_back(static_cast<float>(v0.PCA_Neg()[1]));
    out_branches.Neg.Z->push_back(static_cast<float>(v0.PCA_Neg()[2]));
    out_branches.Neg.Px->push_back(static_cast<float>(v0.Mom_Neg()[0]));
    out_branches.Neg.Py->push_back(static_cast<float>(v0.Mom_Neg()[1]));
    out_branches.Neg.Pz->push_back(static_cast<float>(v0.Mom_Neg()[2]));
    out_branches.Neg.E->push_back(static_cast<float>(v0.Neg_Energy));

    out_branches.Pos.Entry->push_back(v0.idx_pos);
    out_branches.Pos.X->push_back(static_cast<float>(v0.PCA_Pos()[0]));
    out_branches.Pos.Y->push_back(static_cast<float>(v0.PCA_Pos()[1]));
    out_branches.Pos.Z->push_back(static_cast<float>(v0.PCA_Pos()[2]));
    out_branches.Pos.Px->push_back(static_cast<float>(v0.Mom_Pos()[0]));
    out_branches.Pos.Py->push_back(static_cast<float>(v0.Mom_Pos()[1]));
    out_branches.Pos.Pz->push_back(static_cast<float>(v0.Mom_Pos()[2]));
    out_branches.Pos.E->push_back(static_cast<float>(v0.Pos_Energy));

    if (IsMC()) {
        /*
                const int neg_mc_entry{fTruthHandler.McEntry(neg_entry)};
                out_branches.Neg.True.Entry->push_back(neg_mc_entry);
                if (neg_mc_entry > Const::DummyInt) {
                    out_branches.Neg.True.Xv->push_back(fTruthHandler.Xv(neg_mc_entry));
                    out_branches.Neg.True.Yv->push_back(fTruthHandler.Yv(neg_mc_entry));
                    out_branches.Neg.True.Zv->push_back(fTruthHandler.Zv(neg_mc_entry));
                    out_branches.Neg.True.Px->push_back(fTruthHandler.Px(neg_mc_entry));
                    out_branches.Neg.True.Py->push_back(fTruthHandler.Py(neg_mc_entry));
                    out_branches.Neg.True.Pz->push_back(fTruthHandler.Pz(neg_mc_entry));
                    // out_branches.Neg.True.E->push_back(fTruthHandler.E(neg_mc_entry));

                    out_branches.Neg.True.PdgCode->push_back(fTruthHandler.PdgCode(neg_mc_entry));
                    out_branches.Neg.True.MotherEntry->push_back(fTruthHandler.MotherEntry(neg_mc_entry));
                    out_branches.Neg.True.IsSignal->push_back(fTruthHandler.IsSignal(neg_mc_entry));
                    out_branches.Neg.True.IsPrimary->push_back(fTruthHandler.IsPrimary(neg_mc_entry));
                    out_branches.Neg.True.IsSecFromMat->push_back(fTruthHandler.IsSecFromMat(neg_mc_entry));
                    out_branches.Neg.True.IsSecFromWeak->push_back(fTruthHandler.IsSecFromWeak(neg_mc_entry));
                    out_branches.Neg.True.ReactionID->push_back(fTruthHandler.ReactionID(neg_mc_entry));
                }

                const int pos_mc_entry{fTruthHandler.McEntry(pos_entry)};
                out_branches.Pos.True.Entry->push_back(pos_mc_entry);
                if (pos_mc_entry > Const::DummyInt) {
                    out_branches.Pos.True.Xv->push_back(fTruthHandler.Xv(pos_mc_entry));
                    out_branches.Pos.True.Yv->push_back(fTruthHandler.Yv(pos_mc_entry));
                    out_branches.Pos.True.Zv->push_back(fTruthHandler.Zv(pos_mc_entry));
                    out_branches.Pos.True.Px->push_back(fTruthHandler.Px(pos_mc_entry));
                    out_branches.Pos.True.Py->push_back(fTruthHandler.Py(pos_mc_entry));
                    out_branches.Pos.True.Pz->push_back(fTruthHandler.Pz(pos_mc_entry));
                    // out_branches.Pos.True.E->push_back(fTruthHandler.E(pos_mc_entry));

                    out_branches.Pos.True.PdgCode->push_back(fTruthHandler.PdgCode(pos_mc_entry));
                    out_branches.Pos.True.MotherEntry->push_back(fTruthHandler.MotherEntry(pos_mc_entry));
                    out_branches.Pos.True.IsSignal->push_back(fTruthHandler.IsSignal(pos_mc_entry));
                    out_branches.Pos.True.IsPrimary->push_back(fTruthHandler.IsPrimary(pos_mc_entry));
                    out_branches.Pos.True.IsSecFromMat->push_back(fTruthHandler.IsSecFromMat(pos_mc_entry));
                    out_branches.Pos.True.IsSecFromWeak->push_back(fTruthHandler.IsSecFromWeak(pos_mc_entry));
                    out_branches.Pos.True.ReactionID->push_back(fTruthHandler.ReactionID(pos_mc_entry));
                }
         */
        const int v0_mc_entry{fTruthHandler.SameMother(fTruthHandler.McEntry(v0.idx_neg), fTruthHandler.McEntry(v0.idx_pos))};
        out_branches.True.Entry->push_back(v0_mc_entry);
        if (v0_mc_entry > Const::DummyInt) {
            out_branches.True.X->push_back(fTruthHandler.X(v0_mc_entry));
            out_branches.True.Y->push_back(fTruthHandler.Y(v0_mc_entry));
            out_branches.True.Z->push_back(fTruthHandler.Z(v0_mc_entry));
            out_branches.True.Px->push_back(fTruthHandler.Px(v0_mc_entry));
            out_branches.True.Py->push_back(fTruthHandler.Py(v0_mc_entry));
            out_branches.True.Pz->push_back(fTruthHandler.Pz(v0_mc_entry));
            out_branches.True.E->push_back(fTruthHandler.E(v0_mc_entry));

            out_branches.True.PdgCode->push_back(fTruthHandler.PdgCode(v0_mc_entry));
            out_branches.True.MotherEntry->push_back(fTruthHandler.MotherEntry(v0_mc_entry));
            out_branches.True.IsSignal->push_back(fTruthHandler.IsSignal(v0_mc_entry, v0.pdg_code_hyp));
            out_branches.True.IsPrimary->push_back(fTruthHandler.IsPrimary(v0_mc_entry));
            out_branches.True.IsSecFromMat->push_back(fTruthHandler.IsSecFromMat(v0_mc_entry));
            out_branches.True.IsSecFromWeak->push_back(fTruthHandler.IsSecFromWeak(v0_mc_entry));
            out_branches.True.ReactionID->push_back(fTruthHandler.ReactionID(v0_mc_entry));
        } else {
            out_branches.True.X->push_back(Const::DummyFloat);
            out_branches.True.Y->push_back(Const::DummyFloat);
            out_branches.True.Z->push_back(Const::DummyFloat);
            out_branches.True.Px->push_back(Const::DummyFloat);
            out_branches.True.Py->push_back(Const::DummyFloat);
            out_branches.True.Pz->push_back(Const::DummyFloat);
            out_branches.True.E->push_back(Const::DummyFloat);

            out_branches.True.PdgCode->push_back(Const::DummyInt);
            out_branches.True.MotherEntry->push_back(Const::DummyInt);
            out_branches.True.IsSignal->push_back(false);
            out_branches.True.IsPrimary->push_back(false);
            out_branches.True.IsSecFromMat->push_back(false);
            out_branches.True.IsSecFromWeak->push_back(false);
            out_branches.True.ReactionID->push_back(Const::DummyInt);
        }
    }
}

// ## END OF CYCLES ## //

void Packager::EndOfEvent() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    // fill tree //
    fOutputTree->Fill();

    // clear temporary containers //
    if (IsMC()) fTruthHandler.Clear();
    fVec_AntiProtons.clear();
    fVec_Protons.clear();
    fVec_NegKaons.clear();
    fVec_PosKaons.clear();
    fVec_PiMinus.clear();
    fVec_PiPlus.clear();

    // clear output branches //
    switch (GetReactionChannel()) {
        // standard channels //
        case ReactionChannel::A:
            fOutput_AntiLambdas.Clear(IsMC());
            fOutput_KaonsZeroShort.Clear(IsMC());
            break;
        case ReactionChannel::D:
            fOutput_AntiLambdas.Clear(IsMC());
            fOutput_PosKaons.Clear(IsMC());
            break;
        case ReactionChannel::E:
            fOutput_AntiLambdas.Clear(IsMC());
            fOutput_PosKaons.Clear(IsMC());
            fOutput_PiMinus.Clear(IsMC());
            fOutput_PiPlus.Clear(IsMC());
            break;
        case ReactionChannel::H:
            fOutput_PosKaons.Clear(IsMC());
            break;
        // anti-channels //
        case ReactionChannel::AntiA:
            fOutput_Lambdas.Clear(IsMC());
            fOutput_KaonsZeroShort.Clear(IsMC());
            break;
        case ReactionChannel::AntiD:
            fOutput_Lambdas.Clear(IsMC());
            fOutput_NegKaons.Clear(IsMC());
            break;
        case ReactionChannel::AntiE:
            fOutput_Lambdas.Clear(IsMC());
            fOutput_NegKaons.Clear(IsMC());
            fOutput_PiMinus.Clear(IsMC());
            fOutput_PiPlus.Clear(IsMC());
            break;
        case ReactionChannel::AntiH:
            fOutput_NegKaons.Clear(IsMC());
            break;
        // for data //
        case ReactionChannel::All:
            fOutput_AntiLambdas.Clear(IsMC());
            fOutput_Lambdas.Clear(IsMC());
            fOutput_KaonsZeroShort.Clear(IsMC());
            fOutput_NegKaons.Clear(IsMC());
            fOutput_PosKaons.Clear(IsMC());
            fOutput_PiMinus.Clear(IsMC());
            fOutput_PiPlus.Clear(IsMC());
            break;
    }

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Packager::EndOfAnalysis() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    fOutputTree->Write();
    std::cout << "TTree \"" << fOutputTree->GetName() << "\" has been written onto TFile \"" << fSettings.PathOutputFile << "\"" << '\n';

    fEventsTree->ResetBranchAddresses();
    fOutputTree->ResetBranchAddresses();

    std::cout << "Done." << '\n';
#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

}  // namespace Tree2Secondaries
