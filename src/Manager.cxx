#include <iostream>
#include <memory>
#include <string_view>

#include "Analysis/Cuts.hxx"
#include "Analysis/InputFormat.hxx"
#include "Analysis/Manager.hxx"
#include "Math/Constants.hxx"
#include "Math/KFWrapper.hxx"

namespace Tree2Secondaries::Analysis {

bool Manager::Initialize() {
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
    std::cout << "TChain \"Events\" loaded successfully with " << fEventsTree->GetNtrees() << " trees and " << fEventsTree->GetEntries()
              << " total entries." << '\n';

    ConnectInputBranches();

    if (!PrepareOutputFile()) return false;
    if (!PrepareOutputTree()) return false;
    CreateOutputBranches();

    std::cout << "Manager initialized successfully" << '\n';

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
    return true;
}

// ## INPUT ZONE ## //

void Manager::ConnectInputBranches() {
    fEventsTree->SetBranchStatus("*", false);

    ConnectBranchesEvents();
    if (IsMC()) {
        ConnectBranchesMC();
        if (IsSignalMC()) ConnectBranchesInjected();
    }
    ConnectBranchesTracks();
}

void Manager::ConnectBranchesEvents() {
    fEventsTree->SetBranchStatus("RunNumber", true);
    fEventsTree->SetBranchStatus("DirNumber", true);
    fEventsTree->SetBranchStatus("EventNumber", true);
    fEventsTree->SetBranchStatus("Centrality", true);
    fEventsTree->SetBranchStatus("MagneticField", true);
    fEventsTree->SetBranchStatus("PV_Xv", true);
    fEventsTree->SetBranchStatus("PV_Yv", true);
    fEventsTree->SetBranchStatus("PV_Zv", true);

    fEventsTree->SetBranchAddress("RunNumber", &fInput_Event.RunNumber);
    fEventsTree->SetBranchAddress("DirNumber", &fInput_Event.DirNumber);
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

        fEventsTree->SetBranchAddress("MC_PV_Xv", &fInput_Event.PV_TrueXv);
        fEventsTree->SetBranchAddress("MC_PV_Yv", &fInput_Event.PV_TrueYv);
        fEventsTree->SetBranchAddress("MC_PV_Zv", &fInput_Event.PV_TrueZv);
    }
}

void Manager::ConnectBranchesInjected() {
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

void Manager::ConnectBranchesMC() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif
    fEventsTree->SetBranchStatus("MC_Xv", true);
    fEventsTree->SetBranchStatus("MC_Yv", true);
    fEventsTree->SetBranchStatus("MC_Zv", true);
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

    fEventsTree->SetBranchAddress("MC_Xv", &fTruthHandler.fInput_MC.Xv);
    fEventsTree->SetBranchAddress("MC_Yv", &fTruthHandler.fInput_MC.Yv);
    fEventsTree->SetBranchAddress("MC_Zv", &fTruthHandler.fInput_MC.Zv);
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

void Manager::ConnectBranchesTracks() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif
    fEventsTree->SetBranchStatus("Track_Px", true);
    fEventsTree->SetBranchStatus("Track_Py", true);
    fEventsTree->SetBranchStatus("Track_Pz", true);
    fEventsTree->SetBranchStatus("Track_X", true);
    fEventsTree->SetBranchStatus("Track_Y", true);
    fEventsTree->SetBranchStatus("Track_Z", true);
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

    fEventsTree->SetBranchAddress("Track_Px", &fInput_Tracks.Px);
    fEventsTree->SetBranchAddress("Track_Py", &fInput_Tracks.Py);
    fEventsTree->SetBranchAddress("Track_Pz", &fInput_Tracks.Pz);
    fEventsTree->SetBranchAddress("Track_X", &fInput_Tracks.X);
    fEventsTree->SetBranchAddress("Track_Y", &fInput_Tracks.Y);
    fEventsTree->SetBranchAddress("Track_Z", &fInput_Tracks.Z);
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

bool Manager::PrepareOutputFile() {
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

bool Manager::PrepareOutputTree() {
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

void Manager::CreateOutputBranches() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    CreateOutputBranchesEvents();
    if (IsSignalMC()) CreateOutputBranchesInjected();

    switch (GetReactionChannel()) {
        // standard channels //
        case ReactionChannel::A:
            CreateOutputBranchesV0s(Acronym::AntiLambda, fOutput_AntiLambdas);
            CreateOutputBranchesV0s(Acronym::KaonZeroShort, fOutput_KaonsZeroShort);
            break;
        case ReactionChannel::D:
            CreateOutputBranchesV0s(Acronym::AntiLambda, fOutput_AntiLambdas);
            CreateOutputBranchesTracks(Acronym::PosKaon, fOutput_PosKaons);
            break;
        case ReactionChannel::E:
            CreateOutputBranchesV0s(Acronym::AntiLambda, fOutput_AntiLambdas);
            CreateOutputBranchesTracks(Acronym::PosKaon, fOutput_PosKaons);
            CreateOutputBranchesTracks(Acronym::PiMinus, fOutput_PiMinus);
            CreateOutputBranchesTracks(Acronym::PiPlus, fOutput_PiPlus);
            break;
        case ReactionChannel::H:
            CreateOutputBranchesTracks(Acronym::PosKaon, fOutput_PosKaons);
            break;
        // anti-channels //
        case ReactionChannel::AntiA:
            CreateOutputBranchesV0s(Acronym::Lambda, fOutput_Lambdas);
            CreateOutputBranchesV0s(Acronym::KaonZeroShort, fOutput_KaonsZeroShort);
            break;
        case ReactionChannel::AntiD:
            CreateOutputBranchesV0s(Acronym::Lambda, fOutput_Lambdas);
            CreateOutputBranchesTracks(Acronym::NegKaon, fOutput_NegKaons);
            break;
        case ReactionChannel::AntiE:
            CreateOutputBranchesV0s(Acronym::Lambda, fOutput_Lambdas);
            CreateOutputBranchesTracks(Acronym::NegKaon, fOutput_NegKaons);
            CreateOutputBranchesTracks(Acronym::PiMinus, fOutput_PiMinus);
            CreateOutputBranchesTracks(Acronym::PiPlus, fOutput_PiPlus);
            break;
        case ReactionChannel::AntiH:
            CreateOutputBranchesTracks(Acronym::NegKaon, fOutput_NegKaons);
            break;
        // for data //
        case ReactionChannel::All:
            CreateOutputBranchesV0s(Acronym::AntiLambda, fOutput_AntiLambdas);
            CreateOutputBranchesV0s(Acronym::Lambda, fOutput_Lambdas);
            CreateOutputBranchesV0s(Acronym::KaonZeroShort, fOutput_KaonsZeroShort);
            CreateOutputBranchesTracks(Acronym::NegKaon, fOutput_NegKaons);
            CreateOutputBranchesTracks(Acronym::PosKaon, fOutput_PosKaons);
            CreateOutputBranchesTracks(Acronym::PiMinus, fOutput_PiMinus);
            CreateOutputBranchesTracks(Acronym::PiPlus, fOutput_PiPlus);
            break;
    }  // end of switch statement
#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Manager::CreateOutputBranchesEvents() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    fOutputTree->Branch("RunNumber", &fOutput_Event.RunNumber);
    fOutputTree->Branch("DirNumber", &fOutput_Event.DirNumber);
    fOutputTree->Branch("EventNumber", &fOutput_Event.EventNumber);
    fOutputTree->Branch("Centrality", &fOutput_Event.Centrality);
    fOutputTree->Branch("MagneticField", &fOutput_Event.MagneticField);
    fOutputTree->Branch("PV_Xv", &fOutput_Event.PV_Xv);
    fOutputTree->Branch("PV_Yv", &fOutput_Event.PV_Yv);
    fOutputTree->Branch("PV_Zv", &fOutput_Event.PV_Zv);

    if (IsMC()) {
        fOutputTree->Branch("MC_PV_Xv", &fOutput_Event.PV_TrueXv);
        fOutputTree->Branch("MC_PV_Yv", &fOutput_Event.PV_TrueYv);
        fOutputTree->Branch("MC_PV_Zv", &fOutput_Event.PV_TrueZv);
    }
#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Manager::CreateOutputBranchesInjected() {
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

void Manager::CreateOutputBranchesV0s(std::string_view v0_sv, OutputSOA::V0s& out_branches) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    const std::string v0_name{v0_sv};

    fOutputTree->Branch((v0_name + "_Entry").c_str(), &out_branches.Entry);
    fOutputTree->Branch((v0_name + "_Xv").c_str(), &out_branches.Xv);
    fOutputTree->Branch((v0_name + "_Yv").c_str(), &out_branches.Yv);
    fOutputTree->Branch((v0_name + "_Zv").c_str(), &out_branches.Zv);
    fOutputTree->Branch((v0_name + "_Px").c_str(), &out_branches.Px);
    fOutputTree->Branch((v0_name + "_Py").c_str(), &out_branches.Py);
    fOutputTree->Branch((v0_name + "_Pz").c_str(), &out_branches.Pz);
    fOutputTree->Branch((v0_name + "_E").c_str(), &out_branches.E);

    fOutputTree->Branch((v0_name + "_SigmaX2").c_str(), &out_branches.SigmaX2);
    fOutputTree->Branch((v0_name + "_SigmaXY").c_str(), &out_branches.SigmaXY);
    fOutputTree->Branch((v0_name + "_SigmaY2").c_str(), &out_branches.SigmaY2);
    fOutputTree->Branch((v0_name + "_SigmaXZ").c_str(), &out_branches.SigmaXZ);
    fOutputTree->Branch((v0_name + "_SigmaYZ").c_str(), &out_branches.SigmaYZ);
    fOutputTree->Branch((v0_name + "_SigmaZ2").c_str(), &out_branches.SigmaZ2);
    fOutputTree->Branch((v0_name + "_SigmaXPx").c_str(), &out_branches.SigmaXPx);
    fOutputTree->Branch((v0_name + "_SigmaYPx").c_str(), &out_branches.SigmaYPx);
    fOutputTree->Branch((v0_name + "_SigmaZPx").c_str(), &out_branches.SigmaZPx);
    fOutputTree->Branch((v0_name + "_SigmaPx2").c_str(), &out_branches.SigmaPx2);
    fOutputTree->Branch((v0_name + "_SigmaXPy").c_str(), &out_branches.SigmaXPy);
    fOutputTree->Branch((v0_name + "_SigmYaPy").c_str(), &out_branches.SigmaYPy);
    fOutputTree->Branch((v0_name + "_SigmaZPy").c_str(), &out_branches.SigmaZPy);
    fOutputTree->Branch((v0_name + "_SigmaPxPy").c_str(), &out_branches.SigmaPxPy);
    fOutputTree->Branch((v0_name + "_SigmaPy2").c_str(), &out_branches.SigmaPy2);
    fOutputTree->Branch((v0_name + "_SigmaXPz").c_str(), &out_branches.SigmaXPz);
    fOutputTree->Branch((v0_name + "_SigmaYPz").c_str(), &out_branches.SigmaYPz);
    fOutputTree->Branch((v0_name + "_SigmaZPz").c_str(), &out_branches.SigmaZPz);
    fOutputTree->Branch((v0_name + "_SigmaPxPz").c_str(), &out_branches.SigmaPxPz);
    fOutputTree->Branch((v0_name + "_SigmaPyPz").c_str(), &out_branches.SigmaPyPz);
    fOutputTree->Branch((v0_name + "_SigmaPz2").c_str(), &out_branches.SigmaPz2);
    fOutputTree->Branch((v0_name + "_SigmaXE").c_str(), &out_branches.SigmaXE);
    fOutputTree->Branch((v0_name + "_SigmaYE").c_str(), &out_branches.SigmaYE);
    fOutputTree->Branch((v0_name + "_SigmaZE").c_str(), &out_branches.SigmaZE);
    fOutputTree->Branch((v0_name + "_SigmaPxE").c_str(), &out_branches.SigmaPxE);
    fOutputTree->Branch((v0_name + "_SigmaPyE").c_str(), &out_branches.SigmaPyE);
    fOutputTree->Branch((v0_name + "_SigmaPzE").c_str(), &out_branches.SigmaPzE);
    fOutputTree->Branch((v0_name + "_SigmaE2").c_str(), &out_branches.SigmaE2);

    fOutputTree->Branch((v0_name + "_Neg_Entry").c_str(), &out_branches.Neg_Entry);
    fOutputTree->Branch((v0_name + "_Neg_Xv").c_str(), &out_branches.Neg_Xv);
    fOutputTree->Branch((v0_name + "_Neg_Yv").c_str(), &out_branches.Neg_Yv);
    fOutputTree->Branch((v0_name + "_Neg_Zv").c_str(), &out_branches.Neg_Zv);
    fOutputTree->Branch((v0_name + "_Neg_Px").c_str(), &out_branches.Neg_Px);
    fOutputTree->Branch((v0_name + "_Neg_Py").c_str(), &out_branches.Neg_Py);
    fOutputTree->Branch((v0_name + "_Neg_Pz").c_str(), &out_branches.Neg_Pz);

    fOutputTree->Branch((v0_name + "_Pos_Entry").c_str(), &out_branches.Pos_Entry);
    fOutputTree->Branch((v0_name + "_Pos_Xv").c_str(), &out_branches.Pos_Xv);
    fOutputTree->Branch((v0_name + "_Pos_Yv").c_str(), &out_branches.Pos_Yv);
    fOutputTree->Branch((v0_name + "_Pos_Zv").c_str(), &out_branches.Pos_Zv);
    fOutputTree->Branch((v0_name + "_Pos_Px").c_str(), &out_branches.Pos_Px);
    fOutputTree->Branch((v0_name + "_Pos_Py").c_str(), &out_branches.Pos_Py);
    fOutputTree->Branch((v0_name + "_Pos_Pz").c_str(), &out_branches.Pos_Pz);

    if (IsMC()) {
        fOutputTree->Branch((v0_name + "_MC_Entry").c_str(), &out_branches.True.Entry);

        fOutputTree->Branch((v0_name + "_MC_Xv").c_str(), &out_branches.True.Xv);
        fOutputTree->Branch((v0_name + "_MC_Yv").c_str(), &out_branches.True.Yv);
        fOutputTree->Branch((v0_name + "_MC_Zv").c_str(), &out_branches.True.Zv);
        fOutputTree->Branch((v0_name + "_MC_Px").c_str(), &out_branches.True.Px);
        fOutputTree->Branch((v0_name + "_MC_Py").c_str(), &out_branches.True.Py);
        fOutputTree->Branch((v0_name + "_MC_Pz").c_str(), &out_branches.True.Pz);
        fOutputTree->Branch((v0_name + "_MC_E").c_str(), &out_branches.True.E);

        fOutputTree->Branch((v0_name + "_MC_PdgCode").c_str(), &out_branches.True.PdgCode);
        fOutputTree->Branch((v0_name + "_MC_MotherEntry").c_str(), &out_branches.True.MotherEntry);
        fOutputTree->Branch((v0_name + "_MC_IsSignal").c_str(), &out_branches.True.IsSignal);
        fOutputTree->Branch((v0_name + "_MC_IsPrimary").c_str(), &out_branches.True.IsPrimary);
        fOutputTree->Branch((v0_name + "_MC_IsSecFromMat").c_str(), &out_branches.True.IsSecFromMat);
        fOutputTree->Branch((v0_name + "_MC_IsSecFromWeak").c_str(), &out_branches.True.IsSecFromWeak);
        fOutputTree->Branch((v0_name + "_MC_ReactionID").c_str(), &out_branches.True.ReactionID);
        /*
                fOutputTree->Branch((v0_name + "_MC_Neg_Entry").c_str(), &out_branches.Neg.True.Entry);

                fOutputTree->Branch((v0_name + "_MC_Neg_Xv").c_str(), &out_branches.Neg.True.Xv);
                fOutputTree->Branch((v0_name + "_MC_Neg_Yv").c_str(), &out_branches.Neg.True.Yv);
                fOutputTree->Branch((v0_name + "_MC_Neg_Zv").c_str(), &out_branches.Neg.True.Zv);
                fOutputTree->Branch((v0_name + "_MC_Neg_Px").c_str(), &out_branches.Neg.True.Px);
                fOutputTree->Branch((v0_name + "_MC_Neg_Py").c_str(), &out_branches.Neg.True.Py);
                fOutputTree->Branch((v0_name + "_MC_Neg_Pz").c_str(), &out_branches.Neg.True.Pz);
                fOutputTree->Branch((v0_name + "_MC_Neg_E").c_str(), &out_branches.Neg.True.E);

                fOutputTree->Branch((v0_name + "_MC_Neg_PdgCode").c_str(), &out_branches.Neg.True.PdgCode);
                fOutputTree->Branch((v0_name + "_MC_Neg_MotherEntry").c_str(), &out_branches.Neg.True.MotherEntry);
                fOutputTree->Branch((v0_name + "_MC_Neg_IsSignal").c_str(), &out_branches.Neg.True.IsSignal);
                fOutputTree->Branch((v0_name + "_MC_Neg_IsPrimary").c_str(), &out_branches.Neg.True.IsPrimary);
                fOutputTree->Branch((v0_name + "_MC_Neg_IsSecFromMat").c_str(), &out_branches.Neg.True.IsSecFromMat);
                fOutputTree->Branch((v0_name + "_MC_Neg_IsSecFromWeak").c_str(), &out_branches.Neg.True.IsSecFromWeak);
                fOutputTree->Branch((v0_name + "_MC_Neg_ReactionID").c_str(), &out_branches.Neg.True.ReactionID);

                fOutputTree->Branch((v0_name + "_MC_Pos_Entry").c_str(), &out_branches.Pos.True.Entry);

                fOutputTree->Branch((v0_name + "_MC_Pos_Xv").c_str(), &out_branches.Pos.True.Xv);
                fOutputTree->Branch((v0_name + "_MC_Pos_Yv").c_str(), &out_branches.Pos.True.Yv);
                fOutputTree->Branch((v0_name + "_MC_Pos_Zv").c_str(), &out_branches.Pos.True.Zv);
                fOutputTree->Branch((v0_name + "_MC_Pos_Px").c_str(), &out_branches.Pos.True.Px);
                fOutputTree->Branch((v0_name + "_MC_Pos_Py").c_str(), &out_branches.Pos.True.Py);
                fOutputTree->Branch((v0_name + "_MC_Pos_Pz").c_str(), &out_branches.Pos.True.Pz);
                fOutputTree->Branch((v0_name + "_MC_Pos_E").c_str(), &out_branches.Pos.True.E);

                fOutputTree->Branch((v0_name + "_MC_Pos_PdgCode").c_str(), &out_branches.Pos.True.PdgCode);
                fOutputTree->Branch((v0_name + "_MC_Pos_MotherEntry").c_str(), &out_branches.Pos.True.MotherEntry);
                fOutputTree->Branch((v0_name + "_MC_Pos_IsSignal").c_str(), &out_branches.Pos.True.IsSignal);
                fOutputTree->Branch((v0_name + "_MC_Pos_IsPrimary").c_str(), &out_branches.Pos.True.IsPrimary);
                fOutputTree->Branch((v0_name + "_MC_Pos_IsSecFromMat").c_str(), &out_branches.Pos.True.IsSecFromMat);
                fOutputTree->Branch((v0_name + "_MC_Pos_IsSecFromWeak").c_str(), &out_branches.Pos.True.IsSecFromWeak);
                fOutputTree->Branch((v0_name + "_MC_Pos_ReactionID").c_str(), &out_branches.Pos.True.ReactionID);
                */
    }
#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Manager::CreateOutputBranchesTracks(std::string_view charged_sv, OutputSOA::Tracks& out_branches) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    const std::string charged_name{charged_sv};

    fOutputTree->Branch((charged_name + "_Entry").c_str(), &out_branches.Entry);

    fOutputTree->Branch((charged_name + "_Xv").c_str(), &out_branches.Xv);
    fOutputTree->Branch((charged_name + "_Yv").c_str(), &out_branches.Yv);
    fOutputTree->Branch((charged_name + "_Zv").c_str(), &out_branches.Zv);
    fOutputTree->Branch((charged_name + "_Px").c_str(), &out_branches.Px);
    fOutputTree->Branch((charged_name + "_Py").c_str(), &out_branches.Py);
    fOutputTree->Branch((charged_name + "_Pz").c_str(), &out_branches.Pz);
    fOutputTree->Branch((charged_name + "_E").c_str(), &out_branches.E);

    fOutputTree->Branch((charged_name + "_SigmaX2").c_str(), &out_branches.SigmaX2);
    fOutputTree->Branch((charged_name + "_SigmaXY").c_str(), &out_branches.SigmaXY);
    fOutputTree->Branch((charged_name + "_SigmaY2").c_str(), &out_branches.SigmaY2);
    fOutputTree->Branch((charged_name + "_SigmaXZ").c_str(), &out_branches.SigmaXZ);
    fOutputTree->Branch((charged_name + "_SigmaYZ").c_str(), &out_branches.SigmaYZ);
    fOutputTree->Branch((charged_name + "_SigmaZ2").c_str(), &out_branches.SigmaZ2);
    fOutputTree->Branch((charged_name + "_SigmaXPx").c_str(), &out_branches.SigmaXPx);
    fOutputTree->Branch((charged_name + "_SigmaYPx").c_str(), &out_branches.SigmaYPx);
    fOutputTree->Branch((charged_name + "_SigmaZPx").c_str(), &out_branches.SigmaZPx);
    fOutputTree->Branch((charged_name + "_SigmaPx2").c_str(), &out_branches.SigmaPx2);
    fOutputTree->Branch((charged_name + "_SigmaXPy").c_str(), &out_branches.SigmaXPy);
    fOutputTree->Branch((charged_name + "_SigmYaPy").c_str(), &out_branches.SigmaYPy);
    fOutputTree->Branch((charged_name + "_SigmaZPy").c_str(), &out_branches.SigmaZPy);
    fOutputTree->Branch((charged_name + "_SigmaPxPy").c_str(), &out_branches.SigmaPxPy);
    fOutputTree->Branch((charged_name + "_SigmaPy2").c_str(), &out_branches.SigmaPy2);
    fOutputTree->Branch((charged_name + "_SigmaXPz").c_str(), &out_branches.SigmaXPz);
    fOutputTree->Branch((charged_name + "_SigmaYPz").c_str(), &out_branches.SigmaYPz);
    fOutputTree->Branch((charged_name + "_SigmaZPz").c_str(), &out_branches.SigmaZPz);
    fOutputTree->Branch((charged_name + "_SigmaPxPz").c_str(), &out_branches.SigmaPxPz);
    fOutputTree->Branch((charged_name + "_SigmaPyPz").c_str(), &out_branches.SigmaPyPz);
    fOutputTree->Branch((charged_name + "_SigmaPz2").c_str(), &out_branches.SigmaPz2);
    fOutputTree->Branch((charged_name + "_SigmaXE").c_str(), &out_branches.SigmaXE);
    fOutputTree->Branch((charged_name + "_SigmaYE").c_str(), &out_branches.SigmaYE);
    fOutputTree->Branch((charged_name + "_SigmaZE").c_str(), &out_branches.SigmaZE);
    fOutputTree->Branch((charged_name + "_SigmaPxE").c_str(), &out_branches.SigmaPxE);
    fOutputTree->Branch((charged_name + "_SigmaPyE").c_str(), &out_branches.SigmaPyE);
    fOutputTree->Branch((charged_name + "_SigmaPzE").c_str(), &out_branches.SigmaPzE);
    fOutputTree->Branch((charged_name + "_SigmaE2").c_str(), &out_branches.SigmaE2);

    if (IsMC()) {
        fOutputTree->Branch((charged_name + "_MC_Entry").c_str(), &out_branches.True.Entry);

        fOutputTree->Branch((charged_name + "_MC_Xv").c_str(), &out_branches.True.Xv);
        fOutputTree->Branch((charged_name + "_MC_Yv").c_str(), &out_branches.True.Yv);
        fOutputTree->Branch((charged_name + "_MC_Zv").c_str(), &out_branches.True.Zv);
        fOutputTree->Branch((charged_name + "_MC_Px").c_str(), &out_branches.True.Px);
        fOutputTree->Branch((charged_name + "_MC_Py").c_str(), &out_branches.True.Py);
        fOutputTree->Branch((charged_name + "_MC_Pz").c_str(), &out_branches.True.Pz);
        fOutputTree->Branch((charged_name + "_MC_E").c_str(), &out_branches.True.E);

        fOutputTree->Branch((charged_name + "_MC_PdgCode").c_str(), &out_branches.True.PdgCode);
        fOutputTree->Branch((charged_name + "_MC_MotherEntry").c_str(), &out_branches.True.MotherEntry);
        fOutputTree->Branch((charged_name + "_MC_IsSignal").c_str(), &out_branches.True.IsSignal);
        fOutputTree->Branch((charged_name + "_MC_IsPrimary").c_str(), &out_branches.True.IsPrimary);
        fOutputTree->Branch((charged_name + "_MC_IsSecFromMat").c_str(), &out_branches.True.IsSecFromMat);
        fOutputTree->Branch((charged_name + "_MC_IsSecFromWeak").c_str(), &out_branches.True.IsSecFromWeak);
        fOutputTree->Branch((charged_name + "_MC_ReactionID").c_str(), &out_branches.True.ReactionID);
    }
#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Manager::ProcessEvent() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif
    fOutput_Event.RunNumber = fInput_Event.RunNumber;
    fOutput_Event.DirNumber = fInput_Event.DirNumber;
    fOutput_Event.EventNumber = fInput_Event.EventNumber;
    fOutput_Event.Centrality = fInput_Event.Centrality;
    fOutput_Event.MagneticField = fInput_Event.MagneticField;
    fOutput_Event.PV_Xv = fInput_Event.PV_Xv;
    fOutput_Event.PV_Yv = fInput_Event.PV_Yv;
    fOutput_Event.PV_Zv = fInput_Event.PV_Zv;

    if (IsMC()) {
        fOutput_Event.PV_TrueXv = fInput_Event.PV_TrueXv;
        fOutput_Event.PV_TrueYv = fInput_Event.PV_TrueYv;
        fOutput_Event.PV_TrueZv = fInput_Event.PV_TrueZv;
    }
#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Manager::ProcessInjected() {
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
void Manager::ProcessTracks() {
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
void Manager::StoreTracks(PdgCode pdg_code) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif
    // determine rules based on particle pdg code //
    const std::vector<size_t>* vec;
    OutputSOA::Tracks* out;
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

void Manager::Store(const KF::Track& track, OutputSOA::Tracks& out_branches) {

    out_branches.Entry->push_back(track.idx);

    out_branches.Xv->push_back(static_cast<float>(track.GetParameter(0)));
    out_branches.Yv->push_back(static_cast<float>(track.GetParameter(1)));
    out_branches.Zv->push_back(static_cast<float>(track.GetParameter(2)));
    out_branches.Px->push_back(static_cast<float>(track.GetParameter(3)));
    out_branches.Py->push_back(static_cast<float>(track.GetParameter(4)));
    out_branches.Pz->push_back(static_cast<float>(track.GetParameter(5)));
    out_branches.E->push_back(static_cast<float>(track.GetParameter(6)));

    out_branches.SigmaX2->push_back(static_cast<float>(track.GetCovariance(0)));
    out_branches.SigmaXY->push_back(static_cast<float>(track.GetCovariance(1)));
    out_branches.SigmaY2->push_back(static_cast<float>(track.GetCovariance(2)));
    out_branches.SigmaXZ->push_back(static_cast<float>(track.GetCovariance(3)));
    out_branches.SigmaYZ->push_back(static_cast<float>(track.GetCovariance(4)));
    out_branches.SigmaZ2->push_back(static_cast<float>(track.GetCovariance(5)));
    out_branches.SigmaXPx->push_back(static_cast<float>(track.GetCovariance(6)));
    out_branches.SigmaYPx->push_back(static_cast<float>(track.GetCovariance(7)));
    out_branches.SigmaZPx->push_back(static_cast<float>(track.GetCovariance(8)));
    out_branches.SigmaPx2->push_back(static_cast<float>(track.GetCovariance(9)));
    out_branches.SigmaXPy->push_back(static_cast<float>(track.GetCovariance(10)));
    out_branches.SigmaYPy->push_back(static_cast<float>(track.GetCovariance(11)));
    out_branches.SigmaZPy->push_back(static_cast<float>(track.GetCovariance(12)));
    out_branches.SigmaPxPy->push_back(static_cast<float>(track.GetCovariance(13)));
    out_branches.SigmaPy2->push_back(static_cast<float>(track.GetCovariance(14)));
    out_branches.SigmaXPz->push_back(static_cast<float>(track.GetCovariance(15)));
    out_branches.SigmaYPz->push_back(static_cast<float>(track.GetCovariance(16)));
    out_branches.SigmaZPz->push_back(static_cast<float>(track.GetCovariance(17)));
    out_branches.SigmaPxPz->push_back(static_cast<float>(track.GetCovariance(18)));
    out_branches.SigmaPyPz->push_back(static_cast<float>(track.GetCovariance(19)));
    out_branches.SigmaPz2->push_back(static_cast<float>(track.GetCovariance(20)));
    out_branches.SigmaXE->push_back(static_cast<float>(track.GetCovariance(21)));
    out_branches.SigmaYE->push_back(static_cast<float>(track.GetCovariance(22)));
    out_branches.SigmaZE->push_back(static_cast<float>(track.GetCovariance(23)));
    out_branches.SigmaPxE->push_back(static_cast<float>(track.GetCovariance(24)));
    out_branches.SigmaPyE->push_back(static_cast<float>(track.GetCovariance(25)));
    out_branches.SigmaPzE->push_back(static_cast<float>(track.GetCovariance(26)));
    out_branches.SigmaE2->push_back(static_cast<float>(track.GetCovariance(27)));

    if (IsMC()) {
        const int mc_entry{fTruthHandler.McEntry(track.idx)};
        out_branches.True.Entry->push_back(mc_entry);
        if (mc_entry == Const::DummyInt) return;
        out_branches.True.Xv->push_back(fTruthHandler.Xv(mc_entry));
        out_branches.True.Yv->push_back(fTruthHandler.Yv(mc_entry));
        out_branches.True.Zv->push_back(fTruthHandler.Zv(mc_entry));
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

void Manager::FindV0s(PdgCode pdg_code_v0) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    // determine rules based on V0 pdg code //
    const std::vector<size_t>* vec_neg;
    const std::vector<size_t>* vec_pos;
    OutputSOA::V0s* out;
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
            v0.idx = v0_entry;
            v0.idx_neg = esd_neg;
            v0.idx_pos = esd_pos;
            v0.pdg_code_hyp = static_cast<int>(pdg_code_v0);
            v0.AddDaughter(neg, fInput_Event.MagneticField);
            v0.AddDaughter(pos, fInput_Event.MagneticField);
            v0.AddMassConstraint(mass_v0);

            // apply cuts //
            if (!PassesV0Cuts(v0, pdg_code_v0)) continue;
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

bool Manager::PassesLambdaCuts(const KF::V0& v0) const {
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

bool Manager::PassesKaonZeroCuts(const KF::V0& v0) const {
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

void Manager::Store(const KF::V0& v0, OutputSOA::V0s& out_branches) {

    out_branches.Entry->push_back(v0.idx);

    out_branches.Xv->push_back(static_cast<float>(v0.X()));
    out_branches.Yv->push_back(static_cast<float>(v0.Y()));
    out_branches.Zv->push_back(static_cast<float>(v0.Z()));
    out_branches.Px->push_back(static_cast<float>(v0.Px()));
    out_branches.Py->push_back(static_cast<float>(v0.Py()));
    out_branches.Pz->push_back(static_cast<float>(v0.Pz()));
    out_branches.E->push_back(static_cast<float>(v0.E()));

    out_branches.SigmaX2->push_back(static_cast<float>(v0.GetCovariance(0)));
    out_branches.SigmaXY->push_back(static_cast<float>(v0.GetCovariance(1)));
    out_branches.SigmaY2->push_back(static_cast<float>(v0.GetCovariance(2)));
    out_branches.SigmaXZ->push_back(static_cast<float>(v0.GetCovariance(3)));
    out_branches.SigmaYZ->push_back(static_cast<float>(v0.GetCovariance(4)));
    out_branches.SigmaZ2->push_back(static_cast<float>(v0.GetCovariance(5)));
    out_branches.SigmaXPx->push_back(static_cast<float>(v0.GetCovariance(6)));
    out_branches.SigmaYPx->push_back(static_cast<float>(v0.GetCovariance(7)));
    out_branches.SigmaZPx->push_back(static_cast<float>(v0.GetCovariance(8)));
    out_branches.SigmaPx2->push_back(static_cast<float>(v0.GetCovariance(9)));
    out_branches.SigmaXPy->push_back(static_cast<float>(v0.GetCovariance(10)));
    out_branches.SigmaYPy->push_back(static_cast<float>(v0.GetCovariance(11)));
    out_branches.SigmaZPy->push_back(static_cast<float>(v0.GetCovariance(12)));
    out_branches.SigmaPxPy->push_back(static_cast<float>(v0.GetCovariance(13)));
    out_branches.SigmaPy2->push_back(static_cast<float>(v0.GetCovariance(14)));
    out_branches.SigmaXPz->push_back(static_cast<float>(v0.GetCovariance(15)));
    out_branches.SigmaYPz->push_back(static_cast<float>(v0.GetCovariance(16)));
    out_branches.SigmaZPz->push_back(static_cast<float>(v0.GetCovariance(17)));
    out_branches.SigmaPxPz->push_back(static_cast<float>(v0.GetCovariance(18)));
    out_branches.SigmaPyPz->push_back(static_cast<float>(v0.GetCovariance(19)));
    out_branches.SigmaPz2->push_back(static_cast<float>(v0.GetCovariance(20)));
    out_branches.SigmaXE->push_back(static_cast<float>(v0.GetCovariance(21)));
    out_branches.SigmaYE->push_back(static_cast<float>(v0.GetCovariance(22)));
    out_branches.SigmaZE->push_back(static_cast<float>(v0.GetCovariance(23)));
    out_branches.SigmaPxE->push_back(static_cast<float>(v0.GetCovariance(24)));
    out_branches.SigmaPyE->push_back(static_cast<float>(v0.GetCovariance(25)));
    out_branches.SigmaPzE->push_back(static_cast<float>(v0.GetCovariance(26)));
    out_branches.SigmaE2->push_back(static_cast<float>(v0.GetCovariance(27)));

    out_branches.Neg_Entry->push_back(v0.idx_neg);
    out_branches.Neg_Xv->push_back(static_cast<float>(v0.PCA_Neg()[0]));
    out_branches.Neg_Yv->push_back(static_cast<float>(v0.PCA_Neg()[1]));
    out_branches.Neg_Zv->push_back(static_cast<float>(v0.PCA_Neg()[2]));
    out_branches.Neg_Px->push_back(static_cast<float>(v0.Mom_Neg()[0]));
    out_branches.Neg_Py->push_back(static_cast<float>(v0.Mom_Neg()[1]));
    out_branches.Neg_Pz->push_back(static_cast<float>(v0.Mom_Neg()[2]));

    out_branches.Pos_Entry->push_back(v0.idx_pos);
    out_branches.Pos_Xv->push_back(static_cast<float>(v0.PCA_Pos()[0]));
    out_branches.Pos_Yv->push_back(static_cast<float>(v0.PCA_Pos()[1]));
    out_branches.Pos_Zv->push_back(static_cast<float>(v0.PCA_Pos()[2]));
    out_branches.Pos_Px->push_back(static_cast<float>(v0.Mom_Pos()[0]));
    out_branches.Pos_Py->push_back(static_cast<float>(v0.Mom_Pos()[1]));
    out_branches.Pos_Pz->push_back(static_cast<float>(v0.Mom_Pos()[2]));

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
            out_branches.True.Xv->push_back(fTruthHandler.Xv(v0_mc_entry));
            out_branches.True.Yv->push_back(fTruthHandler.Yv(v0_mc_entry));
            out_branches.True.Zv->push_back(fTruthHandler.Zv(v0_mc_entry));
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
            out_branches.True.Xv->push_back(Const::DummyFloat);
            out_branches.True.Yv->push_back(Const::DummyFloat);
            out_branches.True.Zv->push_back(Const::DummyFloat);
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

void Manager::EndOfEvent() {
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

void Manager::EndOfAnalysis() {
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

}  // namespace Tree2Secondaries::Analysis
