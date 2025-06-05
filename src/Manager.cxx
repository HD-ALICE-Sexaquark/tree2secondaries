#include <memory>
#include <string_view>

#include "Math/Point3D.h"
#include "Math/Vector4D.h"

#include "Analysis/Cuts.hxx"
#include "Analysis/InputFormat.hxx"
#include "Analysis/Manager.hxx"
#include "Math/Constants.hxx"
#include "Math/Vertexer.hxx"
#include "Utilities/Logger.hxx"

using XYZPoint = ROOT::Math::XYZPoint;
using PxPyPzMVector = ROOT::Math::PxPyPzMVector;
using PxPyPzEVector = ROOT::Math::PxPyPzEVector;

namespace Tree2Secondaries::Analysis {

bool Manager::Initialize() {

    fEventsTree = std::make_unique<TChain>("Events");
    for (const auto& path : fSettings.PathInputFiles) {
        if (fEventsTree->Add(path.c_str()) == 0) {
            WARNING("Couldn't add TFile \"%s\"", path.c_str());
        }
    }
    if (!fEventsTree->GetEntries()) {
        ERROR("Couldn't manage to read any entry.");
        return false;
    }
    INFO("TChain \"Events\" loaded successfully with %i trees and %lld total entries.", fEventsTree->GetNtrees(), fEventsTree->GetEntries());

    ConnectInputBranches();

    if (!PrepareOutputFile()) return false;
    if (!PrepareOutputTree()) return false;
    CreateOutputBranches();

    INFO("Manager initialized successfully");
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
}

void Manager::ConnectBranchesMC() {
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
}

void Manager::ConnectBranchesTracks() {
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

    if (IsMC()) {
        fEventsTree->SetBranchStatus("Track_McEntry", true);

        fEventsTree->SetBranchAddress("Track_McEntry", &fInput_Tracks.McEntry);
    }
}

// ## OUTPUT ZONE ## //

bool Manager::PrepareOutputFile() {

    fOutputFile = std::unique_ptr<TFile>(TFile::Open(fSettings.PathOutputFile.c_str(), "RECREATE"));
    if (!fOutputFile) {
        ERROR("TFile \"%s\" couldn't be created", fSettings.PathOutputFile.c_str());
        return false;
    }

    return true;
}

bool Manager::PrepareOutputTree() {

    fOutputTree = std::make_unique<TTree>("ProEvents", "Processed Events");
    if (!fOutputTree) {
        ERROR("TTree \"ProEvents\" couldn't be created");
        return false;
    }

    return true;
}

void Manager::CreateOutputBranches() {

    CreateOutputBranchesEvents();
    if (IsSignalMC()) CreateOutputBranchesInjected();

    switch (GetReactionChannel()) {
        // standard channels //
        case ReactionChannel::A:
            CreateOutputBranchesV0s(Acronym::AntiLambda, fOutput_AntiLambdas);
            CreateOutputBranchesV0s(Acronym::KaonZeroShort, fOutput_NeutralKaons);
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
            CreateOutputBranchesV0s(Acronym::KaonZeroShort, fOutput_NeutralKaons);
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
            CreateOutputBranchesV0s(Acronym::KaonZeroShort, fOutput_NeutralKaons);
            CreateOutputBranchesTracks(Acronym::NegKaon, fOutput_NegKaons);
            CreateOutputBranchesTracks(Acronym::PosKaon, fOutput_PosKaons);
            CreateOutputBranchesTracks(Acronym::PiMinus, fOutput_PiMinus);
            CreateOutputBranchesTracks(Acronym::PiPlus, fOutput_PiPlus);
            break;
    }  // end of switch statement
}

void Manager::CreateOutputBranchesEvents() {
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
}

void Manager::CreateOutputBranchesInjected() {
    fOutputTree->Branch("ReactionID", &fOutput_Injected.ReactionID);
    fOutputTree->Branch("Sexaquark_Px", &fOutput_Injected.Px);
    fOutputTree->Branch("Sexaquark_Py", &fOutput_Injected.Py);
    fOutputTree->Branch("Sexaquark_Pz", &fOutput_Injected.Pz);
    fOutputTree->Branch("Nucleon_Px", &fOutput_Injected.Nucleon_Px);
    fOutputTree->Branch("Nucleon_Py", &fOutput_Injected.Nucleon_Py);
    fOutputTree->Branch("Nucleon_Pz", &fOutput_Injected.Nucleon_Pz);
}

void Manager::CreateOutputBranchesV0s(std::string_view v0_sv, OutputSOA::V0s& out_branches) {
    const std::string v0_name{v0_sv};

    fOutputTree->Branch((v0_name + "_Entry").c_str(), &out_branches.Entry);

    fOutputTree->Branch((v0_name + "_Xv").c_str(), &out_branches.Xv);
    fOutputTree->Branch((v0_name + "_Yv").c_str(), &out_branches.Yv);
    fOutputTree->Branch((v0_name + "_Zv").c_str(), &out_branches.Zv);
    fOutputTree->Branch((v0_name + "_Px").c_str(), &out_branches.Px);
    fOutputTree->Branch((v0_name + "_Py").c_str(), &out_branches.Py);
    fOutputTree->Branch((v0_name + "_Pz").c_str(), &out_branches.Pz);
    fOutputTree->Branch((v0_name + "_E").c_str(), &out_branches.E);

    fOutputTree->Branch((v0_name + "_Neg_Entry").c_str(), &out_branches.Neg.Entry);

    fOutputTree->Branch((v0_name + "_Neg_Xv").c_str(), &out_branches.Neg.Xv);
    fOutputTree->Branch((v0_name + "_Neg_Yv").c_str(), &out_branches.Neg.Yv);
    fOutputTree->Branch((v0_name + "_Neg_Zv").c_str(), &out_branches.Neg.Zv);
    fOutputTree->Branch((v0_name + "_Neg_Px").c_str(), &out_branches.Neg.Px);
    fOutputTree->Branch((v0_name + "_Neg_Py").c_str(), &out_branches.Neg.Py);
    fOutputTree->Branch((v0_name + "_Neg_Pz").c_str(), &out_branches.Neg.Pz);
    fOutputTree->Branch((v0_name + "_Neg_E").c_str(), &out_branches.Neg.E);

    fOutputTree->Branch((v0_name + "_Pos_Entry").c_str(), &out_branches.Pos.Entry);

    fOutputTree->Branch((v0_name + "_Pos_Xv").c_str(), &out_branches.Pos.Xv);
    fOutputTree->Branch((v0_name + "_Pos_Yv").c_str(), &out_branches.Pos.Yv);
    fOutputTree->Branch((v0_name + "_Pos_Zv").c_str(), &out_branches.Pos.Zv);
    fOutputTree->Branch((v0_name + "_Pos_Px").c_str(), &out_branches.Pos.Px);
    fOutputTree->Branch((v0_name + "_Pos_Py").c_str(), &out_branches.Pos.Py);
    fOutputTree->Branch((v0_name + "_Pos_Pz").c_str(), &out_branches.Pos.Pz);
    fOutputTree->Branch((v0_name + "_Pos_E").c_str(), &out_branches.Pos.E);

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
    }
}

void Manager::CreateOutputBranchesTracks(std::string_view charged_sv, OutputSOA::Tracks& out_branches) {
    const std::string charged_name{charged_sv};

    fOutputTree->Branch((charged_name + "_Entry").c_str(), &out_branches.Entry);

    fOutputTree->Branch((charged_name + "_Xv").c_str(), &out_branches.Xv);
    fOutputTree->Branch((charged_name + "_Yv").c_str(), &out_branches.Yv);
    fOutputTree->Branch((charged_name + "_Zv").c_str(), &out_branches.Zv);
    fOutputTree->Branch((charged_name + "_Px").c_str(), &out_branches.Px);
    fOutputTree->Branch((charged_name + "_Py").c_str(), &out_branches.Py);
    fOutputTree->Branch((charged_name + "_Pz").c_str(), &out_branches.Pz);
    fOutputTree->Branch((charged_name + "_E").c_str(), &out_branches.E);

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
}

void Manager::ProcessEvent() {
    fPropagator.SetBz(fInput_Event.MagneticField);
    fPropagator.SetPrimaryVertex(fInput_Event.PV_Xv, fInput_Event.PV_Yv, fInput_Event.PV_Zv);

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
}

void Manager::ProcessInjected() {
    fOutput_Injected.ReactionID = fInput_Injected.ReactionID;
    fOutput_Injected.Px = fInput_Injected.Px;
    fOutput_Injected.Py = fInput_Injected.Py;
    fOutput_Injected.Pz = fInput_Injected.Pz;
    fOutput_Injected.Nucleon_Px = fInput_Injected.Nucleon_Px;
    fOutput_Injected.Nucleon_Py = fInput_Injected.Nucleon_Py;
    fOutput_Injected.Nucleon_Pz = fInput_Injected.Nucleon_Pz;
}

// ## Tracks ZONE ## //

void Manager::ProcessTracks() {
    const int n_tracks{NumberTracks()};
    if (IsMC()) fTruthHandler.InitMap(n_tracks);
    for (auto track_entry{0}; track_entry < n_tracks; ++track_entry) {
        // get properties //
        auto charge{fInput_Tracks.Charge->at(track_entry)};
        XYZPoint xyz0{fInput_Tracks.X->at(track_entry), fInput_Tracks.Y->at(track_entry), fInput_Tracks.Z->at(track_entry)};
        auto px0{fInput_Tracks.Px->at(track_entry)};
        auto py0{fInput_Tracks.Py->at(track_entry)};
        auto pz0{fInput_Tracks.Pz->at(track_entry)};
        // PID //
        if (std::abs(fInput_Tracks.NSigmaProton->at(track_entry)) < Cuts::Track::AbsMax_PID_NSigma) {
            auto proton = std::make_shared<Charged>(track_entry, charge, xyz0, PxPyPzMVector{px0, py0, pz0, PdgMass::Proton});
            if (charge < 0) fVec_AntiProtons.push_back(proton);
            if (charge > 0) fVec_Protons.push_back(proton);
        }
        if (std::abs(fInput_Tracks.NSigmaKaon->at(track_entry)) < Cuts::Track::AbsMax_PID_NSigma) {
            auto kaon = std::make_shared<Charged>(track_entry, charge, xyz0, PxPyPzMVector{px0, py0, pz0, PdgMass::Kaon});
            if (charge < 0) fVec_NegKaons.push_back(kaon);
            if (charge > 0) fVec_PosKaons.push_back(kaon);
        }
        if (std::abs(fInput_Tracks.NSigmaPion->at(track_entry)) < Cuts::Track::AbsMax_PID_NSigma) {
            auto pion = std::make_shared<Charged>(track_entry, charge, xyz0, PxPyPzMVector{px0, py0, pz0, PdgMass::Pion});
            if (charge < 0) fVec_PiMinus.push_back(pion);
            if (charge > 0) fVec_PiPlus.push_back(pion);
        }
        // MC //
        if (IsMC()) fTruthHandler.Link(track_entry, fInput_Tracks.McEntry->at(track_entry));
    }
#ifdef T2S_DEBUG
    INFO("finished ProcessTracks() with n_protons = %zu, n_poskaons = %zu, n_piplus = %zu", fProtons.size(), fPosKaons.size(), fPiPlus.size());
    INFO("finished ProcessTracks() with n_antiprotons = %zu, n_negkaons = %zu, n_piminus = %zu", fAntiProtons.size(), fNegKaons.size(),
         fPiMinus.size());
#endif
}

void Manager::StoreTracks(const std::vector<std::shared_ptr<Charged>>& charged_vec, OutputSOA::Tracks& out_branch) {
    for (const auto& track : charged_vec) {
        const int track_entry{track->Entry()};
        out_branch.Entry->push_back(track_entry);

        out_branch.Xv->push_back(static_cast<float>(track->X0()));
        out_branch.Yv->push_back(static_cast<float>(track->Y0()));
        out_branch.Zv->push_back(static_cast<float>(track->Z0()));
        out_branch.Px->push_back(static_cast<float>(track->Px0()));
        out_branch.Py->push_back(static_cast<float>(track->Py0()));
        out_branch.Pz->push_back(static_cast<float>(track->Pz0()));
        out_branch.E->push_back(static_cast<float>(track->Energy()));

        if (IsMC()) {
            const int mc_entry{fTruthHandler.McEntry(track_entry)};
            out_branch.True.Entry->push_back(mc_entry);
            if (mc_entry > Const::DummyInt) {
                out_branch.True.Xv->push_back(fTruthHandler.Xv(mc_entry));
                out_branch.True.Yv->push_back(fTruthHandler.Yv(mc_entry));
                out_branch.True.Zv->push_back(fTruthHandler.Zv(mc_entry));
                out_branch.True.Px->push_back(fTruthHandler.Px(mc_entry));
                out_branch.True.Py->push_back(fTruthHandler.Py(mc_entry));
                out_branch.True.Pz->push_back(fTruthHandler.Pz(mc_entry));
                out_branch.True.E->push_back(fTruthHandler.E(mc_entry));

                out_branch.True.PdgCode->push_back(fTruthHandler.PdgCode(mc_entry));
                out_branch.True.MotherEntry->push_back(fTruthHandler.MotherEntry(mc_entry));
                out_branch.True.IsSignal->push_back(fTruthHandler.IsSignal(mc_entry));
                out_branch.True.IsPrimary->push_back(fTruthHandler.IsPrimary(mc_entry));
                out_branch.True.IsSecFromMat->push_back(fTruthHandler.IsSecFromMat(mc_entry));
                out_branch.True.IsSecFromWeak->push_back(fTruthHandler.IsSecFromWeak(mc_entry));
                out_branch.True.ReactionID->push_back(fTruthHandler.ReactionID(mc_entry));
            } else {
                out_branch.True.Xv->push_back(Const::DummyFloat);
                out_branch.True.Yv->push_back(Const::DummyFloat);
                out_branch.True.Zv->push_back(Const::DummyFloat);
                out_branch.True.Px->push_back(Const::DummyFloat);
                out_branch.True.Py->push_back(Const::DummyFloat);
                out_branch.True.Pz->push_back(Const::DummyFloat);
                out_branch.True.E->push_back(Const::DummyFloat);

                out_branch.True.PdgCode->push_back(Const::DummyInt);
                out_branch.True.MotherEntry->push_back(Const::DummyInt);
                out_branch.True.IsSignal->push_back(false);
                out_branch.True.IsPrimary->push_back(false);
                out_branch.True.IsSecFromMat->push_back(false);
                out_branch.True.IsSecFromWeak->push_back(false);
                out_branch.True.ReactionID->push_back(Const::DummyInt);
            }
        }
    }
}

// ## V0s ZONE ## //

void Manager::FindV0s(int pdg_code_v0, int pdg_code_neg, int pdg_code_pos) {
    // choose tracks species to loop over //
    const auto& neg_vec{pdg_code_neg == PdgCode::AntiProton ? fVec_AntiProtons : fVec_PiMinus};
    const auto& pos_vec{pdg_code_pos == PdgCode::Proton ? fVec_Protons : fVec_PiPlus};
    const double dca_threshold{std::abs(pdg_code_v0) == PdgCode::Lambda ? Cuts::Lambda::Max_DCAbtwDau : Cuts::KaonZeroShort::Max_DCAbtwDau};
    // loop over all possible pairs of tracks //
    int v0_entry{0};
    for (const auto& neg : neg_vec) {
        for (const auto& pos : pos_vec) {
            // sanity check //
            if (neg->Entry() == pos->Entry()) continue;
            // fit //
            Particle::Pair res{Vertexer::MinimizeDistanceHelixHelix(*neg, *pos, fPropagator, dca_threshold)};
            auto V0 = std::make_shared<Neutral>(v0_entry, neg->Entry(), pos->Entry(), pdg_code_v0, res);
            // apply cuts //
            if (!PassesV0Cuts(V0, pdg_code_v0)) continue;
#ifdef T2S_DEBUG
            INFO("neg%i pos%i m%f dbd%f z%f r%f dn%f dp%f pt%f et%f qt%f a%f cpv%f dpv%f", fInput_EsdIdx->at(neg->Entry()),
                 fInput_EsdIdx->at(pos->Entry()), V0->Mass(), V0->DCAbtwDaughters(), V0->DecayZ(), V0->DecayRadius(), V0->DCANegWrtV0(),
                 V0->DCAPosWrtV0(), V0->Pt(), V0->Eta(), V0->ArmenterosQt(), V0->ArmenterosAlpha(), V0->CPAwrt(fPropagator.PrimaryVertex()),
                 V0->DCAwrt(fPropagator.PrimaryVertex()));
#endif
            // store //
            Store(V0, pdg_code_v0);
            ++v0_entry;
        }  // end of loop over pos
    }  // end of loop over neg
}

bool Manager::PassesLambdaCuts(const std::shared_ptr<Neutral>& v0) const {
#ifdef T2S_DEBUG
    INFO("m%f dbd%f z%f r%f dn%f dp%f pt%f et%f qt%f a%f cpv%f dpv%f", v0->Mass(), v0->DCAbtwDaughters(), v0->DecayZ(), v0->DecayRadius(),
         v0->DCANegWrtV0(), v0->DCAPosWrtV0(), v0->Pt(), v0->Eta(), v0->ArmenterosQt(), v0->ArmenterosAlpha(),
         v0->CPAwrt(fPropagator.PrimaryVertex()), v0->DCAwrt(fPropagator.PrimaryVertex()));
#endif
    if (v0->Mass() < Cuts::Lambda::Min_Mass || v0->Mass() > Cuts::Lambda::Max_Mass) return false;
    if (v0->DCAbtwDaughters() > Cuts::Lambda::Max_DCAbtwDau) return false;
    if (std::abs(v0->DecayZ()) > Cuts::Lambda::AbsMax_Zv) return false;
    if (v0->DecayRadius() < Cuts::Lambda::Min_Radius || v0->DecayRadius() > Cuts::Lambda::Max_Radius) return false;
    if (v0->DCANegWrtV0() > Cuts::Lambda::Max_DCAnegV0) return false;
    if (v0->DCAPosWrtV0() > Cuts::Lambda::Max_DCAposV0) return false;
    if (v0->Pt() < Cuts::Lambda::Min_Pt) return false;
    if (std::abs(v0->Eta()) > Cuts::Lambda::AbsMax_Eta) return false;
    if (v0->ArmenterosQt() / std::abs(v0->ArmenterosAlpha()) > Cuts::Lambda::AbsMax_ArmQtOverAlpha) return false;
    if (v0->CPAwrt(fPropagator.PrimaryVertex()) < Cuts::Lambda::Min_CPAwrtPV ||
        v0->CPAwrt(fPropagator.PrimaryVertex()) > Cuts::Lambda::Max_CPAwrtPV) {
        return false;
    }
    if (v0->DCAwrt(fPropagator.PrimaryVertex()) < Cuts::Lambda::Min_DCAwrtPV) return false;

    return true;
}

bool Manager::PassesKaonZeroCuts(const std::shared_ptr<Neutral>& v0) const {
#ifdef T2S_DEBUG
    INFO("m%f dbd%f z%f r%f dn%f dp%f pt%f et%f qt%f a%f cpv%f dpv%f", v0->Mass(), v0->DCAbtwDaughters(), v0->DecayZ(), v0->DecayRadius(),
         v0->DCANegWrtV0(), v0->DCAPosWrtV0(), v0->Pt(), v0->Eta(), v0->ArmenterosQt(), v0->ArmenterosAlpha(),
         v0->CPAwrt(fPropagator.PrimaryVertex()), v0->DCAwrt(fPropagator.PrimaryVertex()));
#endif
    if (v0->DCAbtwDaughters() > Cuts::KaonZeroShort::Max_DCAbtwDau) return false;
    if (v0->Pt() < Cuts::KaonZeroShort::Min_Pt) return false;
    if (v0->Mass() < Cuts::KaonZeroShort::Min_Mass || v0->Mass() > Cuts::KaonZeroShort::Max_Mass) return false;
    if (std::abs(v0->Eta()) > Cuts::KaonZeroShort::AbsMax_Eta) return false;
    if (std::abs(v0->DecayZ()) > Cuts::KaonZeroShort::AbsMax_Zv) return false;
    if (v0->DecayRadius() < Cuts::KaonZeroShort::Min_Radius || v0->DecayRadius() > Cuts::KaonZeroShort::Max_Radius) return false;
    if (v0->DCANegWrtV0() > Cuts::KaonZeroShort::Max_DCAnegV0) return false;
    if (v0->DCAPosWrtV0() > Cuts::KaonZeroShort::Max_DCAposV0) return false;
    if (v0->CPAwrt(fPropagator.PrimaryVertex()) < Cuts::KaonZeroShort::Min_CPAwrtPV ||
        v0->CPAwrt(fPropagator.PrimaryVertex()) > Cuts::KaonZeroShort::Max_CPAwrtPV) {
        return false;
    }
    if (v0->DCAwrt(fPropagator.PrimaryVertex()) < Cuts::KaonZeroShort::Min_DCAwrtPV) return false;

    return true;
}

void Manager::Store(const std::shared_ptr<Neutral>& v0, OutputSOA::V0s& out_branches) {

    out_branches.Entry->push_back(v0->Entry());

    out_branches.Xv->push_back(v0->DecayX());
    out_branches.Yv->push_back(v0->DecayY());
    out_branches.Zv->push_back(v0->DecayZ());
    out_branches.Px->push_back(v0->Px());
    out_branches.Py->push_back(v0->Py());
    out_branches.Pz->push_back(v0->Pz());
    out_branches.E->push_back(v0->Energy());

    const int neg_entry{v0->NegEntry()};
    out_branches.Neg.Entry->push_back(v0->NegEntry());

    out_branches.Neg.Xv->push_back(v0->NegVertex().X());
    out_branches.Neg.Yv->push_back(v0->NegVertex().Y());
    out_branches.Neg.Zv->push_back(v0->NegVertex().Z());
    out_branches.Neg.Px->push_back(v0->NegPxPyPzE().Px());
    out_branches.Neg.Py->push_back(v0->NegPxPyPzE().Py());
    out_branches.Neg.Pz->push_back(v0->NegPxPyPzE().Pz());
    out_branches.Neg.E->push_back(v0->NegPxPyPzE().E());

    const int pos_entry{v0->PosEntry()};
    out_branches.Pos.Entry->push_back(v0->PosEntry());

    out_branches.Pos.Xv->push_back(v0->PosVertex().X());
    out_branches.Pos.Yv->push_back(v0->PosVertex().Y());
    out_branches.Pos.Zv->push_back(v0->PosVertex().Z());
    out_branches.Pos.Px->push_back(v0->PosPxPyPzE().Px());
    out_branches.Pos.Py->push_back(v0->PosPxPyPzE().Py());
    out_branches.Pos.Pz->push_back(v0->PosPxPyPzE().Pz());
    out_branches.Pos.E->push_back(v0->PosPxPyPzE().E());

    if (IsMC()) {
        const int neg_mc_entry{fTruthHandler.McEntry(neg_entry)};
        out_branches.Neg.True.Entry->push_back(neg_mc_entry);
        if (neg_mc_entry > Const::DummyInt) {
            out_branches.Neg.True.Xv->push_back(fTruthHandler.Xv(neg_mc_entry));
            out_branches.Neg.True.Yv->push_back(fTruthHandler.Yv(neg_mc_entry));
            out_branches.Neg.True.Zv->push_back(fTruthHandler.Zv(neg_mc_entry));
            out_branches.Neg.True.Px->push_back(fTruthHandler.Px(neg_mc_entry));
            out_branches.Neg.True.Py->push_back(fTruthHandler.Py(neg_mc_entry));
            out_branches.Neg.True.Pz->push_back(fTruthHandler.Pz(neg_mc_entry));
            out_branches.Neg.True.E->push_back(fTruthHandler.E(neg_mc_entry));

            out_branches.Neg.True.PdgCode->push_back(fTruthHandler.PdgCode(neg_mc_entry));
            out_branches.Neg.True.MotherEntry->push_back(fTruthHandler.MotherEntry(neg_mc_entry));
            out_branches.Neg.True.IsSignal->push_back(fTruthHandler.IsSignal(neg_mc_entry));
            out_branches.Neg.True.IsPrimary->push_back(fTruthHandler.IsPrimary(neg_mc_entry));
            out_branches.Neg.True.IsSecFromMat->push_back(fTruthHandler.IsSecFromMat(neg_mc_entry));
            out_branches.Neg.True.IsSecFromWeak->push_back(fTruthHandler.IsSecFromWeak(neg_mc_entry));
            out_branches.Neg.True.ReactionID->push_back(fTruthHandler.ReactionID(neg_mc_entry));
        } else {
            out_branches.Neg.True.Xv->push_back(Const::DummyFloat);
            out_branches.Neg.True.Yv->push_back(Const::DummyFloat);
            out_branches.Neg.True.Zv->push_back(Const::DummyFloat);
            out_branches.Neg.True.Px->push_back(Const::DummyFloat);
            out_branches.Neg.True.Py->push_back(Const::DummyFloat);
            out_branches.Neg.True.Pz->push_back(Const::DummyFloat);
            out_branches.Neg.True.E->push_back(Const::DummyFloat);

            out_branches.Neg.True.PdgCode->push_back(Const::DummyInt);
            out_branches.Neg.True.MotherEntry->push_back(Const::DummyInt);
            out_branches.Neg.True.IsSignal->push_back(false);
            out_branches.Neg.True.IsPrimary->push_back(false);
            out_branches.Neg.True.IsSecFromMat->push_back(false);
            out_branches.Neg.True.IsSecFromWeak->push_back(false);
            out_branches.Neg.True.ReactionID->push_back(Const::DummyInt);
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
            out_branches.Pos.True.E->push_back(fTruthHandler.E(pos_mc_entry));

            out_branches.Pos.True.PdgCode->push_back(fTruthHandler.PdgCode(pos_mc_entry));
            out_branches.Pos.True.MotherEntry->push_back(fTruthHandler.MotherEntry(pos_mc_entry));
            out_branches.Pos.True.IsSignal->push_back(fTruthHandler.IsSignal(pos_mc_entry));
            out_branches.Pos.True.IsPrimary->push_back(fTruthHandler.IsPrimary(pos_mc_entry));
            out_branches.Pos.True.IsSecFromMat->push_back(fTruthHandler.IsSecFromMat(pos_mc_entry));
            out_branches.Pos.True.IsSecFromWeak->push_back(fTruthHandler.IsSecFromWeak(pos_mc_entry));
            out_branches.Pos.True.ReactionID->push_back(fTruthHandler.ReactionID(pos_mc_entry));
        } else {
            out_branches.Pos.True.Xv->push_back(Const::DummyFloat);
            out_branches.Pos.True.Yv->push_back(Const::DummyFloat);
            out_branches.Pos.True.Zv->push_back(Const::DummyFloat);
            out_branches.Pos.True.Px->push_back(Const::DummyFloat);
            out_branches.Pos.True.Py->push_back(Const::DummyFloat);
            out_branches.Pos.True.Pz->push_back(Const::DummyFloat);
            out_branches.Pos.True.E->push_back(Const::DummyFloat);

            out_branches.Pos.True.PdgCode->push_back(Const::DummyInt);
            out_branches.Pos.True.MotherEntry->push_back(Const::DummyInt);
            out_branches.Pos.True.IsSignal->push_back(false);
            out_branches.Pos.True.IsPrimary->push_back(false);
            out_branches.Pos.True.IsSecFromMat->push_back(false);
            out_branches.Pos.True.IsSecFromWeak->push_back(false);
            out_branches.Pos.True.ReactionID->push_back(Const::DummyInt);
        }

        const int v0_mc_entry{fTruthHandler.SameMother(neg_mc_entry, pos_mc_entry)};
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
            out_branches.True.IsSignal->push_back(fTruthHandler.IsSignal(v0_mc_entry, v0->HypothesisPID()));
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
            fOutput_NeutralKaons.Clear(IsMC());
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
            fOutput_NeutralKaons.Clear(IsMC());
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
            fOutput_NeutralKaons.Clear(IsMC());
            fOutput_NegKaons.Clear(IsMC());
            fOutput_PosKaons.Clear(IsMC());
            fOutput_PiMinus.Clear(IsMC());
            fOutput_PiPlus.Clear(IsMC());
            break;
    }
}

void Manager::EndOfAnalysis() {
    fOutputTree->Write();
    INFO("TTree \"%s\" has been written onto TFile \"%s\".", fOutputTree->GetName(), fSettings.PathOutputFile.c_str());

    fEventsTree->ResetBranchAddresses();
    fOutputTree->ResetBranchAddresses();

    INFO("Done.");
}

}  // namespace Tree2Secondaries::Analysis
