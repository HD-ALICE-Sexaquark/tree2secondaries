#include <iostream>
#include <memory>

#include "App/Utilities.hxx"
#include "KF/Utilities.hxx"
#include "Math/Constants.hxx"
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
            std::cerr << "   " << __FUNCTION__ << " :: Couldn't add TFile \"" << path << "\"" << '\n';
        }
    }
    if (!fEventsTree->GetEntries()) {
        std::cerr << "   " << __FUNCTION__ << " :: Couldn't manage to read any entry." << '\n';
        return false;
    }
    std::cout << "   " << __FUNCTION__ << " :: TChain \"" << fEventsTree->GetName() << "\" loaded successfully with " << fEventsTree->GetNtrees()
              << " trees and " << fEventsTree->GetEntries() << " total entries." << '\n';

    ConnectInputBranches();

    if (!PrepareOutputFile()) return false;
    if (!PrepareOutputTree()) return false;
    CreateOutputBranches();

    std::cout << "   " << __FUNCTION__ << " :: Packager initialized successfully" << '\n';

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
    return true;
}

// ## INPUT ZONE ## //

void Packager::ConnectInputBranches() {
    fEventsTree->SetBranchStatus("*", false);

    ConnectBranches_Events();
    if (IsMC()) {
        ConnectBranches_MC();
        if (IsSignalMC()) ConnectBranches_Injected();
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
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    Utils::ConnectBranch(fEventsTree.get(), "ReactionID", &fInput_Injected.ReactionID);
    Utils::ConnectBranch(fEventsTree.get(), "Sexaquark_Px", &fInput_Injected.Px);
    Utils::ConnectBranch(fEventsTree.get(), "Sexaquark_Py", &fInput_Injected.Py);
    Utils::ConnectBranch(fEventsTree.get(), "Sexaquark_Pz", &fInput_Injected.Pz);
    Utils::ConnectBranch(fEventsTree.get(), "Nucleon_Px", &fInput_Injected.Nucleon_Px);
    Utils::ConnectBranch(fEventsTree.get(), "Nucleon_Py", &fInput_Injected.Nucleon_Py);
    Utils::ConnectBranch(fEventsTree.get(), "Nucleon_Pz", &fInput_Injected.Nucleon_Pz);

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Packager::ConnectBranches_MC() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
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
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Packager::ConnectBranches_Tracks() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
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
        std::cerr << "   " << __FUNCTION__ << " :: TFile \"" << fSettings.PathOutputFile << "\" couldn't be created" << '\n';
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

    CreateOutputBranches_Events();
    if (IsSignalMC()) CreateOutputBranches_Injected();

    switch (GetReactionChannel()) {
        // standard channels //
        case ReactionChannel::A:
            CreateOutputBranches_V0s(static_cast<std::string>(Acronym::AntiLambda), fOutput_AntiLambdas);
            CreateOutputBranches_V0s(static_cast<std::string>(Acronym::KaonZeroShort), fOutput_KaonsZeroShort);
            if (IsMC()) {
                CreateOutputBranches_MC_V0s(static_cast<std::string>(Acronym::AntiLambda), fOutput_MC_AntiLambdas);
                CreateOutputBranches_MC_V0s(static_cast<std::string>(Acronym::KaonZeroShort), fOutput_MC_KaonsZeroShort);
            }
            break;
        case ReactionChannel::D:
            CreateOutputBranches_V0s(static_cast<std::string>(Acronym::AntiLambda), fOutput_AntiLambdas);
            CreateOutputBranches_Tracks(static_cast<std::string>(Acronym::PosKaon), fOutput_PosKaons);
            if (IsMC()) {
                CreateOutputBranches_MC_V0s(static_cast<std::string>(Acronym::AntiLambda), fOutput_MC_AntiLambdas);
                CreateOutputBranches_MC_Tracks(static_cast<std::string>(Acronym::PosKaon), fOutput_MC_PosKaons);
            }
            break;
        case ReactionChannel::E:
            CreateOutputBranches_V0s(static_cast<std::string>(Acronym::AntiLambda), fOutput_AntiLambdas);
            CreateOutputBranches_Tracks(static_cast<std::string>(Acronym::PosKaon), fOutput_PosKaons);
            CreateOutputBranches_Tracks(static_cast<std::string>(Acronym::PiMinus), fOutput_PiMinus);
            CreateOutputBranches_Tracks(static_cast<std::string>(Acronym::PiPlus), fOutput_PiPlus);
            if (IsMC()) {
                CreateOutputBranches_MC_V0s(static_cast<std::string>(Acronym::AntiLambda), fOutput_MC_AntiLambdas);
                CreateOutputBranches_MC_Tracks(static_cast<std::string>(Acronym::PosKaon), fOutput_MC_PosKaons);
                CreateOutputBranches_MC_Tracks(static_cast<std::string>(Acronym::PiMinus), fOutput_MC_PiMinus);
                CreateOutputBranches_MC_Tracks(static_cast<std::string>(Acronym::PiPlus), fOutput_MC_PiPlus);
            }
            break;
        case ReactionChannel::H:
            CreateOutputBranches_Tracks(static_cast<std::string>(Acronym::PosKaon), fOutput_PosKaons);
            if (IsMC()) CreateOutputBranches_MC_Tracks(static_cast<std::string>(Acronym::PosKaon), fOutput_MC_PosKaons);
            break;
        // anti-channels //
        case ReactionChannel::AntiA:
            CreateOutputBranches_V0s(static_cast<std::string>(Acronym::Lambda), fOutput_Lambdas);
            CreateOutputBranches_V0s(static_cast<std::string>(Acronym::KaonZeroShort), fOutput_KaonsZeroShort);
            if (IsMC()) {
                CreateOutputBranches_V0s(static_cast<std::string>(Acronym::Lambda), fOutput_Lambdas);
                CreateOutputBranches_V0s(static_cast<std::string>(Acronym::KaonZeroShort), fOutput_KaonsZeroShort);
            }
            break;
        case ReactionChannel::AntiD:
            CreateOutputBranches_V0s(static_cast<std::string>(Acronym::Lambda), fOutput_Lambdas);
            CreateOutputBranches_Tracks(static_cast<std::string>(Acronym::NegKaon), fOutput_NegKaons);
            if (IsMC()) {
                CreateOutputBranches_MC_V0s(static_cast<std::string>(Acronym::Lambda), fOutput_MC_Lambdas);
                CreateOutputBranches_MC_Tracks(static_cast<std::string>(Acronym::NegKaon), fOutput_MC_NegKaons);
            }
            break;
        case ReactionChannel::AntiE:
            CreateOutputBranches_V0s(static_cast<std::string>(Acronym::Lambda), fOutput_Lambdas);
            CreateOutputBranches_Tracks(static_cast<std::string>(Acronym::NegKaon), fOutput_NegKaons);
            CreateOutputBranches_Tracks(static_cast<std::string>(Acronym::PiMinus), fOutput_PiMinus);
            CreateOutputBranches_Tracks(static_cast<std::string>(Acronym::PiPlus), fOutput_PiPlus);
            if (IsMC()) {
                CreateOutputBranches_MC_V0s(static_cast<std::string>(Acronym::Lambda), fOutput_MC_Lambdas);
                CreateOutputBranches_MC_Tracks(static_cast<std::string>(Acronym::NegKaon), fOutput_MC_NegKaons);
                CreateOutputBranches_MC_Tracks(static_cast<std::string>(Acronym::PiMinus), fOutput_MC_PiMinus);
                CreateOutputBranches_MC_Tracks(static_cast<std::string>(Acronym::PiPlus), fOutput_MC_PiPlus);
            }
            break;
        case ReactionChannel::AntiH:
            CreateOutputBranches_Tracks(static_cast<std::string>(Acronym::NegKaon), fOutput_NegKaons);
            if (IsMC()) CreateOutputBranches_MC_Tracks(static_cast<std::string>(Acronym::NegKaon), fOutput_MC_NegKaons);
            break;
        // for data //
        case ReactionChannel::All:
            CreateOutputBranches_V0s(static_cast<std::string>(Acronym::AntiLambda), fOutput_AntiLambdas);
            CreateOutputBranches_V0s(static_cast<std::string>(Acronym::Lambda), fOutput_Lambdas);
            CreateOutputBranches_V0s(static_cast<std::string>(Acronym::KaonZeroShort), fOutput_KaonsZeroShort);
            CreateOutputBranches_Tracks(static_cast<std::string>(Acronym::NegKaon), fOutput_NegKaons);
            CreateOutputBranches_Tracks(static_cast<std::string>(Acronym::PosKaon), fOutput_PosKaons);
            CreateOutputBranches_Tracks(static_cast<std::string>(Acronym::PiMinus), fOutput_PiMinus);
            CreateOutputBranches_Tracks(static_cast<std::string>(Acronym::PiPlus), fOutput_PiPlus);
            if (IsMC()) {
                CreateOutputBranches_MC_V0s(static_cast<std::string>(Acronym::AntiLambda), fOutput_MC_AntiLambdas);
                CreateOutputBranches_MC_V0s(static_cast<std::string>(Acronym::Lambda), fOutput_MC_Lambdas);
                CreateOutputBranches_MC_V0s(static_cast<std::string>(Acronym::KaonZeroShort), fOutput_MC_KaonsZeroShort);
                CreateOutputBranches_MC_Tracks(static_cast<std::string>(Acronym::NegKaon), fOutput_MC_NegKaons);
                CreateOutputBranches_MC_Tracks(static_cast<std::string>(Acronym::PosKaon), fOutput_MC_PosKaons);
                CreateOutputBranches_MC_Tracks(static_cast<std::string>(Acronym::PiMinus), fOutput_MC_PiMinus);
                CreateOutputBranches_MC_Tracks(static_cast<std::string>(Acronym::PiPlus), fOutput_MC_PiPlus);
            }
            break;
    }  // end of switch statement
#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Packager::CreateOutputBranches_Events() {
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

void Packager::CreateOutputBranches_Injected() {
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

void Packager::CreateOutputBranches_V0s(const std::string& name_v0, PackedEvents::V0s& sov) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    fOutputTree->Branch((name_v0 + "_Entry").c_str(), &sov.Entry);
    fOutputTree->Branch((name_v0 + "_X").c_str(), &sov.X);
    fOutputTree->Branch((name_v0 + "_Y").c_str(), &sov.Y);
    fOutputTree->Branch((name_v0 + "_Z").c_str(), &sov.Z);
    fOutputTree->Branch((name_v0 + "_Px").c_str(), &sov.Px);
    fOutputTree->Branch((name_v0 + "_Py").c_str(), &sov.Py);
    fOutputTree->Branch((name_v0 + "_Pz").c_str(), &sov.Pz);
    fOutputTree->Branch((name_v0 + "_E").c_str(), &sov.E);

    fOutputTree->Branch((name_v0 + "_SigmaX2").c_str(), &sov.Sigma.X2);
    fOutputTree->Branch((name_v0 + "_SigmaXY").c_str(), &sov.Sigma.XY);
    fOutputTree->Branch((name_v0 + "_SigmaY2").c_str(), &sov.Sigma.Y2);
    fOutputTree->Branch((name_v0 + "_SigmaXZ").c_str(), &sov.Sigma.XZ);
    fOutputTree->Branch((name_v0 + "_SigmaYZ").c_str(), &sov.Sigma.YZ);
    fOutputTree->Branch((name_v0 + "_SigmaZ2").c_str(), &sov.Sigma.Z2);
    fOutputTree->Branch((name_v0 + "_SigmaXPx").c_str(), &sov.Sigma.XPx);
    fOutputTree->Branch((name_v0 + "_SigmaYPx").c_str(), &sov.Sigma.YPx);
    fOutputTree->Branch((name_v0 + "_SigmaZPx").c_str(), &sov.Sigma.ZPx);
    fOutputTree->Branch((name_v0 + "_SigmaPx2").c_str(), &sov.Sigma.Px2);
    fOutputTree->Branch((name_v0 + "_SigmaXPy").c_str(), &sov.Sigma.XPy);
    fOutputTree->Branch((name_v0 + "_SigmaYPy").c_str(), &sov.Sigma.YPy);
    fOutputTree->Branch((name_v0 + "_SigmaZPy").c_str(), &sov.Sigma.ZPy);
    fOutputTree->Branch((name_v0 + "_SigmaPxPy").c_str(), &sov.Sigma.PxPy);
    fOutputTree->Branch((name_v0 + "_SigmaPy2").c_str(), &sov.Sigma.Py2);
    fOutputTree->Branch((name_v0 + "_SigmaXPz").c_str(), &sov.Sigma.XPz);
    fOutputTree->Branch((name_v0 + "_SigmaYPz").c_str(), &sov.Sigma.YPz);
    fOutputTree->Branch((name_v0 + "_SigmaZPz").c_str(), &sov.Sigma.ZPz);
    fOutputTree->Branch((name_v0 + "_SigmaPxPz").c_str(), &sov.Sigma.PxPz);
    fOutputTree->Branch((name_v0 + "_SigmaPyPz").c_str(), &sov.Sigma.PyPz);
    fOutputTree->Branch((name_v0 + "_SigmaPz2").c_str(), &sov.Sigma.Pz2);
    fOutputTree->Branch((name_v0 + "_SigmaXE").c_str(), &sov.Sigma.XE);
    fOutputTree->Branch((name_v0 + "_SigmaYE").c_str(), &sov.Sigma.YE);
    fOutputTree->Branch((name_v0 + "_SigmaZE").c_str(), &sov.Sigma.ZE);
    fOutputTree->Branch((name_v0 + "_SigmaPxE").c_str(), &sov.Sigma.PxE);
    fOutputTree->Branch((name_v0 + "_SigmaPyE").c_str(), &sov.Sigma.PyE);
    fOutputTree->Branch((name_v0 + "_SigmaPzE").c_str(), &sov.Sigma.PzE);
    fOutputTree->Branch((name_v0 + "_SigmaE2").c_str(), &sov.Sigma.E2);

    fOutputTree->Branch((name_v0 + "_Neg_Entry").c_str(), &sov.Neg.Entry);
    fOutputTree->Branch((name_v0 + "_Neg_X").c_str(), &sov.Neg.X);
    fOutputTree->Branch((name_v0 + "_Neg_Y").c_str(), &sov.Neg.Y);
    fOutputTree->Branch((name_v0 + "_Neg_Z").c_str(), &sov.Neg.Z);
    fOutputTree->Branch((name_v0 + "_Neg_Px").c_str(), &sov.Neg.Px);
    fOutputTree->Branch((name_v0 + "_Neg_Py").c_str(), &sov.Neg.Py);
    fOutputTree->Branch((name_v0 + "_Neg_Pz").c_str(), &sov.Neg.Pz);
    fOutputTree->Branch((name_v0 + "_Neg_E").c_str(), &sov.Neg.E);

    fOutputTree->Branch((name_v0 + "_Pos_Entry").c_str(), &sov.Pos.Entry);
    fOutputTree->Branch((name_v0 + "_Pos_X").c_str(), &sov.Pos.X);
    fOutputTree->Branch((name_v0 + "_Pos_Y").c_str(), &sov.Pos.Y);
    fOutputTree->Branch((name_v0 + "_Pos_Z").c_str(), &sov.Pos.Z);
    fOutputTree->Branch((name_v0 + "_Pos_Px").c_str(), &sov.Pos.Px);
    fOutputTree->Branch((name_v0 + "_Pos_Py").c_str(), &sov.Pos.Py);
    fOutputTree->Branch((name_v0 + "_Pos_Pz").c_str(), &sov.Pos.Pz);
    fOutputTree->Branch((name_v0 + "_Pos_E").c_str(), &sov.Pos.E);

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Packager::CreateOutputBranches_Tracks(const std::string& name_part, PackedEvents::Tracks& sov) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    fOutputTree->Branch((name_part + "_Entry").c_str(), &sov.Entry);

    fOutputTree->Branch((name_part + "_X").c_str(), &sov.X);
    fOutputTree->Branch((name_part + "_Y").c_str(), &sov.Y);
    fOutputTree->Branch((name_part + "_Z").c_str(), &sov.Z);
    fOutputTree->Branch((name_part + "_Px").c_str(), &sov.Px);
    fOutputTree->Branch((name_part + "_Py").c_str(), &sov.Py);
    fOutputTree->Branch((name_part + "_Pz").c_str(), &sov.Pz);
    fOutputTree->Branch((name_part + "_E").c_str(), &sov.E);

    fOutputTree->Branch((name_part + "_SigmaX2").c_str(), &sov.Sigma.X2);
    fOutputTree->Branch((name_part + "_SigmaXY").c_str(), &sov.Sigma.XY);
    fOutputTree->Branch((name_part + "_SigmaY2").c_str(), &sov.Sigma.Y2);
    fOutputTree->Branch((name_part + "_SigmaXZ").c_str(), &sov.Sigma.XZ);
    fOutputTree->Branch((name_part + "_SigmaYZ").c_str(), &sov.Sigma.YZ);
    fOutputTree->Branch((name_part + "_SigmaZ2").c_str(), &sov.Sigma.Z2);
    fOutputTree->Branch((name_part + "_SigmaXPx").c_str(), &sov.Sigma.XPx);
    fOutputTree->Branch((name_part + "_SigmaYPx").c_str(), &sov.Sigma.YPx);
    fOutputTree->Branch((name_part + "_SigmaZPx").c_str(), &sov.Sigma.ZPx);
    fOutputTree->Branch((name_part + "_SigmaPx2").c_str(), &sov.Sigma.Px2);
    fOutputTree->Branch((name_part + "_SigmaXPy").c_str(), &sov.Sigma.XPy);
    fOutputTree->Branch((name_part + "_SigmaYPy").c_str(), &sov.Sigma.YPy);
    fOutputTree->Branch((name_part + "_SigmaZPy").c_str(), &sov.Sigma.ZPy);
    fOutputTree->Branch((name_part + "_SigmaPxPy").c_str(), &sov.Sigma.PxPy);
    fOutputTree->Branch((name_part + "_SigmaPy2").c_str(), &sov.Sigma.Py2);
    fOutputTree->Branch((name_part + "_SigmaXPz").c_str(), &sov.Sigma.XPz);
    fOutputTree->Branch((name_part + "_SigmaYPz").c_str(), &sov.Sigma.YPz);
    fOutputTree->Branch((name_part + "_SigmaZPz").c_str(), &sov.Sigma.ZPz);
    fOutputTree->Branch((name_part + "_SigmaPxPz").c_str(), &sov.Sigma.PxPz);
    fOutputTree->Branch((name_part + "_SigmaPyPz").c_str(), &sov.Sigma.PyPz);
    fOutputTree->Branch((name_part + "_SigmaPz2").c_str(), &sov.Sigma.Pz2);
    fOutputTree->Branch((name_part + "_SigmaXE").c_str(), &sov.Sigma.XE);
    fOutputTree->Branch((name_part + "_SigmaYE").c_str(), &sov.Sigma.YE);
    fOutputTree->Branch((name_part + "_SigmaZE").c_str(), &sov.Sigma.ZE);
    fOutputTree->Branch((name_part + "_SigmaPxE").c_str(), &sov.Sigma.PxE);
    fOutputTree->Branch((name_part + "_SigmaPyE").c_str(), &sov.Sigma.PyE);
    fOutputTree->Branch((name_part + "_SigmaPzE").c_str(), &sov.Sigma.PzE);
    fOutputTree->Branch((name_part + "_SigmaE2").c_str(), &sov.Sigma.E2);

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Packager::CreateOutputBranches_MC_V0s(const std::string& name_v0, PackedEvents::MC_V0s& sov) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    fOutputTree->Branch((name_v0 + "_MC_Entry").c_str(), &sov.Entry);
    fOutputTree->Branch((name_v0 + "_MC_X").c_str(), &sov.X);
    fOutputTree->Branch((name_v0 + "_MC_Y").c_str(), &sov.Y);
    fOutputTree->Branch((name_v0 + "_MC_Z").c_str(), &sov.Z);
    fOutputTree->Branch((name_v0 + "_MC_Px").c_str(), &sov.Px);
    fOutputTree->Branch((name_v0 + "_MC_Py").c_str(), &sov.Py);
    fOutputTree->Branch((name_v0 + "_MC_Pz").c_str(), &sov.Pz);
    fOutputTree->Branch((name_v0 + "_MC_E").c_str(), &sov.E);

    fOutputTree->Branch((name_v0 + "_MC_DecayX").c_str(), &sov.DecayX);
    fOutputTree->Branch((name_v0 + "_MC_DecayY").c_str(), &sov.DecayY);
    fOutputTree->Branch((name_v0 + "_MC_DecayZ").c_str(), &sov.DecayZ);

    fOutputTree->Branch((name_v0 + "_MC_PdgCode").c_str(), &sov.PdgCode);
    fOutputTree->Branch((name_v0 + "_MC_MotherEntry").c_str(), &sov.MotherEntry);
    fOutputTree->Branch((name_v0 + "_MC_PdgCode_Mother").c_str(), &sov.PdgCode_Mother);
    fOutputTree->Branch((name_v0 + "_MC_IsTrue").c_str(), &sov.IsTrue);
    fOutputTree->Branch((name_v0 + "_MC_IsSignal").c_str(), &sov.IsSignal);
    fOutputTree->Branch((name_v0 + "_MC_IsSecondary").c_str(), &sov.IsSecondary);
    fOutputTree->Branch((name_v0 + "_MC_ReactionID").c_str(), &sov.ReactionID);
    fOutputTree->Branch((name_v0 + "_MC_IsHybrid").c_str(), &sov.IsHybrid);

    fOutputTree->Branch((name_v0 + "_MC_Neg_Entry").c_str(), &sov.Neg_Entry);
    fOutputTree->Branch((name_v0 + "_MC_Neg_Px").c_str(), &sov.Neg_Px);
    fOutputTree->Branch((name_v0 + "_MC_Neg_Py").c_str(), &sov.Neg_Py);
    fOutputTree->Branch((name_v0 + "_MC_Neg_Pz").c_str(), &sov.Neg_Pz);
    fOutputTree->Branch((name_v0 + "_MC_Neg_PdgCode").c_str(), &sov.Neg_PdgCode);
    fOutputTree->Branch((name_v0 + "_MC_Neg_IsTrue").c_str(), &sov.Neg_IsTrue);
    fOutputTree->Branch((name_v0 + "_MC_Neg_IsSignal").c_str(), &sov.Neg_IsSignal);
    fOutputTree->Branch((name_v0 + "_MC_Neg_IsSecondary").c_str(), &sov.Neg_IsSecondary);
    fOutputTree->Branch((name_v0 + "_MC_Neg_ReactionID").c_str(), &sov.Neg_ReactionID);

    fOutputTree->Branch((name_v0 + "_MC_Pos_Entry").c_str(), &sov.Pos_Entry);
    fOutputTree->Branch((name_v0 + "_MC_Pos_Px").c_str(), &sov.Pos_Px);
    fOutputTree->Branch((name_v0 + "_MC_Pos_Py").c_str(), &sov.Pos_Py);
    fOutputTree->Branch((name_v0 + "_MC_Pos_Pz").c_str(), &sov.Pos_Pz);
    fOutputTree->Branch((name_v0 + "_MC_Pos_PdgCode").c_str(), &sov.Pos_PdgCode);
    fOutputTree->Branch((name_v0 + "_MC_Pos_IsTrue").c_str(), &sov.Pos_IsTrue);
    fOutputTree->Branch((name_v0 + "_MC_Pos_IsSignal").c_str(), &sov.Pos_IsSignal);
    fOutputTree->Branch((name_v0 + "_MC_Pos_IsSecondary").c_str(), &sov.Pos_IsSecondary);
    fOutputTree->Branch((name_v0 + "_MC_Pos_ReactionID").c_str(), &sov.Pos_ReactionID);

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
#endif
}

void Packager::CreateOutputBranches_MC_Tracks(const std::string& name_part, PackedEvents::MC_Tracks& sov) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    fOutputTree->Branch((name_part + "_MC_Entry").c_str(), &sov.Entry);
    fOutputTree->Branch((name_part + "_MC_X").c_str(), &sov.X);
    fOutputTree->Branch((name_part + "_MC_Y").c_str(), &sov.Y);
    fOutputTree->Branch((name_part + "_MC_Z").c_str(), &sov.Z);
    fOutputTree->Branch((name_part + "_MC_Px").c_str(), &sov.Px);
    fOutputTree->Branch((name_part + "_MC_Py").c_str(), &sov.Py);
    fOutputTree->Branch((name_part + "_MC_Pz").c_str(), &sov.Pz);
    fOutputTree->Branch((name_part + "_MC_E").c_str(), &sov.E);

    fOutputTree->Branch((name_part + "_MC_MotherEntry").c_str(), &sov.MotherEntry);
    fOutputTree->Branch((name_part + "_MC_GrandMotherEntry").c_str(), &sov.MotherEntry);
    fOutputTree->Branch((name_part + "_MC_PdgCode").c_str(), &sov.PdgCode);
    fOutputTree->Branch((name_part + "_MC_PdgCode_Mother").c_str(), &sov.PdgCode_Mother);
    fOutputTree->Branch((name_part + "_MC_PdgCode_GrandMother").c_str(), &sov.PdgCode_GrandMother);
    fOutputTree->Branch((name_part + "_MC_IsTrue").c_str(), &sov.IsTrue);
    fOutputTree->Branch((name_part + "_MC_IsSignal").c_str(), &sov.IsSignal);
    fOutputTree->Branch((name_part + "_MC_IsSecondary").c_str(), &sov.IsSecondary);
    fOutputTree->Branch((name_part + "_MC_ReactionID").c_str(), &sov.ReactionID);

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

    for (size_t esd_track{0}; esd_track < NumberTracks(); ++esd_track) {
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
    PackedEvents::MC_Tracks* mc_out{nullptr};
    double mass;
    switch (pdg_code) {
        case PdgCode::NegKaon:
            vec = &fVec_NegKaons;
            out = &fOutput_NegKaons;
            if (IsMC()) mc_out = &fOutput_MC_NegKaons;
            mass = PdgMass::Kaon;
            break;
        case PdgCode::PosKaon:
            vec = &fVec_PosKaons;
            out = &fOutput_PosKaons;
            if (IsMC()) mc_out = &fOutput_MC_PosKaons;
            mass = PdgMass::Kaon;
            break;
        case PdgCode::PiMinus:
            vec = &fVec_PiMinus;
            out = &fOutput_PiMinus;
            if (IsMC()) mc_out = &fOutput_MC_PiMinus;
            mass = PdgMass::Pion;
            break;
        case PdgCode::PiPlus:
            vec = &fVec_PiPlus;
            out = &fOutput_PiPlus;
            if (IsMC()) mc_out = &fOutput_MC_PiPlus;
            mass = PdgMass::Pion;
            break;
        default:
            return;
    }

    // loop over selected tracks //
    for (auto esd_idx : *vec) {

        // prepare kf object //
        KF::Vector<6> neg_kf_params = KF::PackParams(fInput_Tracks, esd_idx);
        std::array<float, 5> neg_alice_params = KF::PackParams_ALICE(fInput_Tracks, esd_idx);
        std::array<float, 15> neg_alice_cov = KF::PackCovMatrix_ALICE(fInput_Tracks, esd_idx);
        KF::Track kf_track{KF::CreateParticle(neg_kf_params, neg_alice_params, neg_alice_cov, fInput_Tracks.Alpha->at(esd_idx),
                                              fInput_Tracks.Charge->at(esd_idx), mass),
                           esd_idx};

        // store //
        Store(kf_track, *out);
        if (IsMC()) {
            MC::Track mc_track{fInput_MC, fInput_Tracks.McEntry->at(esd_idx), pdg_code};
            StoreMC(mc_track, *mc_out);
        }
    }  // end of loop over selected tracks

#ifdef T2S_DEBUG
    std::cout << "-- finished (" << __FUNCTION__ << ") --" << '\n';
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

void Packager::StoreMC(const MC::Track& mc, PackedEvents::MC_Tracks& sov) {

    sov.Entry->push_back(mc.Entry);
    sov.X->push_back(mc.X);
    sov.Y->push_back(mc.Y);
    sov.Z->push_back(mc.Z);
    sov.Px->push_back(mc.Px);
    sov.Py->push_back(mc.Py);
    sov.Pz->push_back(mc.Pz);
    sov.E->push_back(mc.Energy);

    sov.PdgCode->push_back(mc.PdgCode);
    sov.MotherEntry->push_back(mc.MotherEntry);
    sov.PdgCode_Mother->push_back(mc.PdgCode_Mother);
    sov.GrandMotherEntry->push_back(mc.GrandMotherEntry);
    sov.PdgCode_GrandMother->push_back(mc.PdgCode_GrandMother);
    sov.IsTrue->push_back(mc.IsTrue);
    sov.IsSignal->push_back(mc.IsSignal);
    sov.IsSecondary->push_back(mc.IsSecondary);
    sov.ReactionID->push_back(mc.ReactionID);
}

// ## V0s ZONE ## //

void Packager::FindV0s(PdgCode pdg_code) {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    // determine rules based on V0 pdg code //
    const std::vector<size_t>* vec_neg;
    const std::vector<size_t>* vec_pos;
    PackedEvents::V0s* out;
    PackedEvents::MC_V0s* mc_out{nullptr};
    double mass_neg;
    double mass_pos;
    double mass_v0;
    PdgCode pdg_code_neg;
    PdgCode pdg_code_pos;
    switch (pdg_code) {
        case PdgCode::AntiLambda:
            vec_neg = &fVec_AntiProtons;
            vec_pos = &fVec_PiPlus;
            out = &fOutput_AntiLambdas;
            if (IsMC()) mc_out = &fOutput_MC_AntiLambdas;
            mass_neg = PdgMass::Proton;
            mass_pos = PdgMass::Pion;
            mass_v0 = PdgMass::Lambda;
            pdg_code_neg = PdgCode::AntiProton;
            pdg_code_pos = PdgCode::PiPlus;
            break;
        case PdgCode::Lambda:
            vec_neg = &fVec_PiMinus;
            vec_pos = &fVec_Protons;
            out = &fOutput_Lambdas;
            if (IsMC()) mc_out = &fOutput_MC_Lambdas;
            mass_neg = PdgMass::Pion;
            mass_pos = PdgMass::Proton;
            mass_v0 = PdgMass::Lambda;
            pdg_code_neg = PdgCode::PiMinus;
            pdg_code_pos = PdgCode::Proton;
            break;
        case PdgCode::KaonZeroShort:
            vec_neg = &fVec_PiMinus;
            vec_pos = &fVec_PiPlus;
            out = &fOutput_KaonsZeroShort;
            if (IsMC()) mc_out = &fOutput_MC_KaonsZeroShort;
            mass_neg = PdgMass::Pion;
            mass_pos = PdgMass::Pion;
            mass_v0 = PdgMass::KaonZeroShort;
            pdg_code_neg = PdgCode::PiMinus;
            pdg_code_pos = PdgCode::PiPlus;
            break;
        default:
            std::cerr << "   " << __FUNCTION__ << " :: Invalid PDG Code " << int(pdg_code) << " for a V0." << '\n';
            return;
    }

    // loop over all possible pairs of tracks //
    size_t v0_entry{0};
    for (auto esd_neg : *vec_neg) {
        for (auto esd_pos : *vec_pos) {

            // sanity check //
            if (esd_neg == esd_pos) continue;

            // prepare neg //
            std::array<double, 6> neg_kf_params = KF::PackParams(fInput_Tracks, esd_neg);
            std::array<float, 5> neg_alice_params = KF::PackParams_ALICE(fInput_Tracks, esd_neg);
            std::array<float, 15> neg_alice_cov = KF::PackCovMatrix_ALICE(fInput_Tracks, esd_neg);
            auto neg = KF::CreateParticle(neg_kf_params, neg_alice_params, neg_alice_cov, fInput_Tracks.Alpha->at(esd_neg),
                                          fInput_Tracks.Charge->at(esd_neg), mass_neg);

            // prepare pos //
            std::array<double, 6> pos_kf_params = KF::PackParams(fInput_Tracks, esd_pos);
            std::array<float, 5> pos_alice_params = KF::PackParams_ALICE(fInput_Tracks, esd_pos);
            std::array<float, 15> pos_alice_cov = KF::PackCovMatrix_ALICE(fInput_Tracks, esd_pos);
            auto pos = KF::CreateParticle(pos_kf_params, pos_alice_params, pos_alice_cov, fInput_Tracks.Alpha->at(esd_pos),
                                          fInput_Tracks.Charge->at(esd_pos), mass_pos);

            // build v0 //
            KF::V0 v0;
            v0.AddDaughter(neg, fInput_Event.MagneticField);
            v0.AddDaughter(pos, fInput_Event.MagneticField);
            v0.AddMassConstraint(mass_v0);

            // apply cuts //
            if (!PassesCuts(v0, pdg_code)) continue;

            // add info //
            v0.SetIndices({v0_entry, esd_neg, esd_pos});
            v0.Neg_Energy = neg.E();
            v0.Pos_Energy = pos.E();
            v0.pdg_code_hyp = static_cast<int>(pdg_code);
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
            if (IsMC()) {
                MC::V0 mc_v0{fInput_MC, fInput_Tracks.McEntry->at(esd_neg), fInput_Tracks.McEntry->at(esd_pos), pdg_code, pdg_code_neg, pdg_code_pos};
                StoreMC(mc_v0, *mc_out);
            }
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

    sov.Neg.Entry->push_back(v0.idx_neg);
    sov.Neg.X->push_back(static_cast<float>(v0.PCA_Neg()[0]));
    sov.Neg.Y->push_back(static_cast<float>(v0.PCA_Neg()[1]));
    sov.Neg.Z->push_back(static_cast<float>(v0.PCA_Neg()[2]));
    sov.Neg.Px->push_back(static_cast<float>(v0.Mom_Neg()[0]));
    sov.Neg.Py->push_back(static_cast<float>(v0.Mom_Neg()[1]));
    sov.Neg.Pz->push_back(static_cast<float>(v0.Mom_Neg()[2]));
    sov.Neg.E->push_back(static_cast<float>(v0.Neg_Energy));

    sov.Pos.Entry->push_back(v0.idx_pos);
    sov.Pos.X->push_back(static_cast<float>(v0.PCA_Pos()[0]));
    sov.Pos.Y->push_back(static_cast<float>(v0.PCA_Pos()[1]));
    sov.Pos.Z->push_back(static_cast<float>(v0.PCA_Pos()[2]));
    sov.Pos.Px->push_back(static_cast<float>(v0.Mom_Pos()[0]));
    sov.Pos.Py->push_back(static_cast<float>(v0.Mom_Pos()[1]));
    sov.Pos.Pz->push_back(static_cast<float>(v0.Mom_Pos()[2]));
    sov.Pos.E->push_back(static_cast<float>(v0.Pos_Energy));
}

void Packager::StoreMC(const MC::V0& v0, PackedEvents::MC_V0s& sov) {

    sov.Entry->push_back(v0.Entry);
    sov.X->push_back(v0.X);
    sov.Y->push_back(v0.Y);
    sov.Z->push_back(v0.Z);
    sov.Px->push_back(v0.Px);
    sov.Py->push_back(v0.Py);
    sov.Pz->push_back(v0.Pz);
    sov.E->push_back(v0.Energy);

    sov.DecayX->push_back(v0.DecayX());
    sov.DecayY->push_back(v0.DecayY());
    sov.DecayZ->push_back(v0.DecayZ());

    sov.PdgCode->push_back(v0.PdgCode);
    sov.MotherEntry->push_back(v0.MotherEntry);
    sov.PdgCode_Mother->push_back(v0.PdgCode_Mother);
    sov.IsTrue->push_back(v0.IsTrue);
    sov.IsSignal->push_back(v0.IsSignal);
    sov.IsSecondary->push_back(v0.IsSecondary);
    sov.ReactionID->push_back(v0.ReactionID);
    sov.IsHybrid->push_back(v0.IsHybrid);

    // neg //
    sov.Neg_Entry->push_back(v0.neg.Entry);
    sov.Neg_Px->push_back(v0.neg.Px);
    sov.Neg_Py->push_back(v0.neg.Py);
    sov.Neg_Pz->push_back(v0.neg.Pz);
    sov.Neg_PdgCode->push_back(v0.neg.PdgCode);
    sov.Neg_IsTrue->push_back(v0.neg.IsTrue);
    sov.Neg_IsSignal->push_back(v0.neg.IsSignal);
    sov.Neg_IsSecondary->push_back(v0.neg.IsSecondary);
    sov.Neg_ReactionID->push_back(v0.neg.ReactionID);

    // pos //
    sov.Pos_Entry->push_back(v0.pos.Entry);
    sov.Pos_Px->push_back(v0.pos.Px);
    sov.Pos_Py->push_back(v0.pos.Py);
    sov.Pos_Pz->push_back(v0.pos.Pz);
    sov.Pos_PdgCode->push_back(v0.pos.PdgCode);
    sov.Pos_IsTrue->push_back(v0.pos.IsTrue);
    sov.Pos_IsSignal->push_back(v0.pos.IsSignal);
    sov.Pos_IsSecondary->push_back(v0.pos.IsSecondary);
    sov.Pos_ReactionID->push_back(v0.pos.ReactionID);
}

// ## END OF CYCLES ## //

void Packager::EndOfEvent() {
#ifdef T2S_DEBUG
    std::cout << "-- starting (" << __FUNCTION__ << ") --" << '\n';
#endif

    // fill tree //
    fOutputTree->Fill();

    // clear temporary containers //
    fVec_AntiProtons.clear();
    fVec_Protons.clear();
    fVec_NegKaons.clear();
    fVec_PosKaons.clear();
    fVec_PiMinus.clear();
    fVec_PiPlus.clear();

    // clear output branches //
    if (IsSignalMC()) fOutput_Injected.Clear();
    switch (GetReactionChannel()) {
        // standard channels //
        case ReactionChannel::A:
            fOutput_AntiLambdas.Clear();
            fOutput_KaonsZeroShort.Clear();
            if (IsMC()) {
                fOutput_AntiLambdas.Clear();
                fOutput_KaonsZeroShort.Clear();
            }
            break;
        case ReactionChannel::D:
            fOutput_AntiLambdas.Clear();
            fOutput_PosKaons.Clear();
            if (IsMC()) {
                fOutput_AntiLambdas.Clear();
                fOutput_PosKaons.Clear();
            }
            break;
        case ReactionChannel::E:
            fOutput_AntiLambdas.Clear();
            fOutput_PosKaons.Clear();
            fOutput_PiMinus.Clear();
            fOutput_PiPlus.Clear();
            if (IsMC()) {
                fOutput_AntiLambdas.Clear();
                fOutput_PosKaons.Clear();
                fOutput_PiMinus.Clear();
                fOutput_PiPlus.Clear();
            }
            break;
        case ReactionChannel::H:
            fOutput_PosKaons.Clear();
            if (IsMC()) fOutput_PosKaons.Clear();
            break;
        // anti-channels //
        case ReactionChannel::AntiA:
            fOutput_Lambdas.Clear();
            fOutput_KaonsZeroShort.Clear();
            if (IsMC()) {
                fOutput_Lambdas.Clear();
                fOutput_KaonsZeroShort.Clear();
            }
            break;
        case ReactionChannel::AntiD:
            fOutput_Lambdas.Clear();
            fOutput_NegKaons.Clear();
            if (IsMC()) {
                fOutput_Lambdas.Clear();
                fOutput_NegKaons.Clear();
            }
            break;
        case ReactionChannel::AntiE:
            fOutput_Lambdas.Clear();
            fOutput_NegKaons.Clear();
            fOutput_PiMinus.Clear();
            fOutput_PiPlus.Clear();
            if (IsMC()) {
                fOutput_Lambdas.Clear();
                fOutput_NegKaons.Clear();
                fOutput_PiMinus.Clear();
                fOutput_PiPlus.Clear();
            }
            break;
        case ReactionChannel::AntiH:
            fOutput_NegKaons.Clear();
            if (IsMC()) fOutput_NegKaons.Clear();
            break;
        // for data //
        case ReactionChannel::All:
            fOutput_AntiLambdas.Clear();
            fOutput_Lambdas.Clear();
            fOutput_KaonsZeroShort.Clear();
            fOutput_NegKaons.Clear();
            fOutput_PosKaons.Clear();
            fOutput_PiMinus.Clear();
            fOutput_PiPlus.Clear();
            if (IsMC()) {
                fOutput_AntiLambdas.Clear();
                fOutput_Lambdas.Clear();
                fOutput_KaonsZeroShort.Clear();
                fOutput_NegKaons.Clear();
                fOutput_PosKaons.Clear();
                fOutput_PiMinus.Clear();
                fOutput_PiPlus.Clear();
            }
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
