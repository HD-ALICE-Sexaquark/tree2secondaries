#include <memory>
#include <set>

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

namespace Tree2Secondaries::Analysis {

// Initialize the manager.
bool Manager::Initialize() {

    if (!OpenInputFile()) return false;
    if (!LoadInputTree()) return false;
    ConnectInputBranches();

    if (!PrepareOutputFile()) return false;
    if (!PrepareOutputTree()) return false;
    PrepareOutputBranches();

    INFO("Manager initialized successfully");
    return true;
}

// ## INPUT ZONE ## //

bool Manager::OpenInputFile() {

    fInputFile = std::unique_ptr<TFile>(TFile::Open(fSettings.PathInputFiles.c_str(), "READ"));
    if (!fInputFile || fInputFile->IsZombie()) {
        ERROR("TFile \"%s\" couldn't be opened", fSettings.PathInputFiles.c_str());
        return false;
    }
    INFO("TFile \"%s\" opened successfully", fSettings.PathInputFiles.c_str());

    return true;
}

// Load the event tree from the input file.
bool Manager::LoadInputTree() {

    fEventsTree = std::unique_ptr<TTree>(fInputFile->Get<TTree>("Events"));
    if (!fEventsTree) {
        ERROR("TTree \"Events\" couldn't be found in TFile \"%s\"", fSettings.PathInputFiles.c_str());
        return false;
    }
    INFO("TTree \"Events\" found in TFile \"%s\"", fSettings.PathInputFiles.c_str());
    INFO("TTree \"Events\" has %lld entries", fEventsTree->GetEntries());

    return true;
}

// Connect branches to the event tree.
void Manager::ConnectInputBranches() {
    fEventsTree->SetBranchStatus("*", false);

    ConnectBranchesEvents();
    if (IsMC()) {
        ConnectBranchesMC();
        if (IsSignalMC()) ConnectBranchesInjected();
    }
    ConnectBranchesTracks();
}

// Connect branches for event properties.
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

// Connect branches for injected interactions.
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

// Connect branches for MC data.
void Manager::ConnectBranchesMC() {
    fEventsTree->SetBranchStatus("MC_PdgCode", true);
    fEventsTree->SetBranchStatus("MC_Mother_McEntry", true);
    fEventsTree->SetBranchStatus("MC_Status", true);
    fEventsTree->SetBranchStatus("MC_Generator", true);

    fEventsTree->SetBranchAddress("MC_PdgCode", &fInput_MC.PdgCode);
    fEventsTree->SetBranchAddress("MC_Mother_McEntry", &fInput_MC.Mother_McEntry);
    fEventsTree->SetBranchAddress("MC_Status", &fInput_MC.Status);
    fEventsTree->SetBranchAddress("MC_Generator", &fInput_MC.Generator);
}

// Connect branches for reconstructed tracks.
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

//
bool Manager::PrepareOutputFile() {

    fOutputFile = std::unique_ptr<TFile>(TFile::Open(fSettings.PathOutputFile.c_str(), "RECREATE"));
    if (!fOutputFile) {
        ERROR("TFile \"%s\" couldn't be created", fSettings.PathOutputFile.c_str());
        return false;
    }
    INFO("TFile \"%s\" (re)created successfully", fSettings.PathOutputFile.c_str());

    return true;
}

// Create output tree.
bool Manager::PrepareOutputTree() {

    fOutputTree = std::make_unique<TTree>("ProEvents", "Processed Events");
    if (!fOutputTree) {
        ERROR("TTree \"ProEvents\" couldn't be created");
        return false;
    }

    return true;
}

void Manager::PrepareOutputBranches() {
    PrepareBranchesEvents();
    PrepareBranchesV0s();
}

void Manager::PrepareBranchesEvents() {
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

//
void Manager::PrepareBranchesV0s() {
    fOutputTree->Branch("V0_Index", &fOutput_V0s.Index);
    fOutputTree->Branch("V0_PID", &fOutput_V0s.PID);
    fOutputTree->Branch("V0_Neg_Entry", &fOutput_V0s.Neg_Entry);
    fOutputTree->Branch("V0_Pos_Entry", &fOutput_V0s.Pos_Entry);
    fOutputTree->Branch("V0_Xv", &fOutput_V0s.Xv);
    fOutputTree->Branch("V0_Yv", &fOutput_V0s.Yv);
    fOutputTree->Branch("V0_Zv", &fOutput_V0s.Zv);
    fOutputTree->Branch("V0_Px", &fOutput_V0s.Px);
    fOutputTree->Branch("V0_Py", &fOutput_V0s.Py);
    fOutputTree->Branch("V0_Pz", &fOutput_V0s.Pz);
    fOutputTree->Branch("V0_E", &fOutput_V0s.E);
    fOutputTree->Branch("V0_Neg_Xv", &fOutput_V0s.Neg_Xv);
    fOutputTree->Branch("V0_Neg_Yv", &fOutput_V0s.Neg_Yv);
    fOutputTree->Branch("V0_Neg_Zv", &fOutput_V0s.Neg_Zv);
    fOutputTree->Branch("V0_Neg_Px", &fOutput_V0s.Neg_Px);
    fOutputTree->Branch("V0_Neg_Py", &fOutput_V0s.Neg_Py);
    fOutputTree->Branch("V0_Neg_Pz", &fOutput_V0s.Neg_Pz);
    fOutputTree->Branch("V0_Pos_Xv", &fOutput_V0s.Pos_Xv);
    fOutputTree->Branch("V0_Pos_Yv", &fOutput_V0s.Pos_Yv);
    fOutputTree->Branch("V0_Pos_Zv", &fOutput_V0s.Pos_Zv);
    fOutputTree->Branch("V0_Pos_Px", &fOutput_V0s.Pos_Px);
    fOutputTree->Branch("V0_Pos_Py", &fOutput_V0s.Pos_Py);
    fOutputTree->Branch("V0_Pos_Pz", &fOutput_V0s.Pos_Pz);
    if (IsMC()) {
        // PENDING
    }
}

// Read and copy the event data.
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

// Process injected interactions.
void Manager::ProcessInjected() {}

// Process the MC data.
/*
void Manager::ProcessMC() {
    auto n_mc{NumberMC()};
    fMCParticles.reserve(n_mc);
    for (auto mc_entry{0}; mc_entry < n_mc; mc_entry++) {
        // Get properties //
        auto pdg_code{fInput_MC.PdgCode->at(mc_entry)};
        auto mother_mc_entry{fInput_MC.Mother_McEntry->at(mc_entry)};
        auto is_signal{fInput_MC.Generator->at(mc_entry) == 2};
        // auto reaction_id{mother_mc_entry != -1 ? fMCParticles[mother_mc_entry]->ReactionID() : fInput_MC.Status->at(mc_entry)};
        auto reaction_id{is_signal ? fInput_MC.Status->at(mc_entry) : -1};
        // if (mother_mc_entry != -1) {
        // reaction_id = fMCParticles[mother_mc_entry]->ReactionID();
        // }
        // Store //
        fMCParticles.emplace_back(std::make_shared<True>(mc_entry, mother_mc_entry, pdg_code, is_signal, reaction_id));
    }
    // DEBUG MC //
    // for (const auto& mc : fMCParticles) {
    // if (mc->ReactionID() == -1) continue;
    // INFO("MC particle %d, pdg code: %d, mother: %d, reaction: %d", mc->Index(), mc->PdgCode(), mc->MotherMcEntry(), mc->ReactionID());
    // }
}
*/

// Process selected reconstructed tracks.
void Manager::ProcessTracks() {

    for (auto track_entry{0}; track_entry < NumberTracks(); track_entry++) {
        // Get properties //
        auto charge{fInput_Tracks.Charge->at(track_entry)};
        XYZPoint xyz0{fInput_Tracks.X->at(track_entry), fInput_Tracks.Y->at(track_entry), fInput_Tracks.Z->at(track_entry)};
        auto px0{fInput_Tracks.Px->at(track_entry)};
        auto py0{fInput_Tracks.Py->at(track_entry)};
        auto pz0{fInput_Tracks.Pz->at(track_entry)};
        // PID //
        if (std::abs(fInput_Tracks.NSigmaProton->at(track_entry)) < Cuts::Track::AbsMax_PID_NSigma) {
            auto proton = std::make_shared<Charged>(track_entry, charge, xyz0, PxPyPzMVector{px0, py0, pz0, Const::MassProton});
            // if (IsMC()) proton->SetLinkedMc(fMCParticles[fInput_Tracks.McEntry->at(track_entry)]);
            if (charge < 0) fAntiProtons.push_back(proton);
            if (charge > 0) fProtons.push_back(proton);
        }
        if (std::abs(fInput_Tracks.NSigmaKaon->at(track_entry)) < Cuts::Track::AbsMax_PID_NSigma) {
            auto kaon = std::make_shared<Charged>(track_entry, charge, xyz0, PxPyPzMVector{px0, py0, pz0, Const::MassKaon});
            // if (IsMC()) kaon->SetLinkedMc(fMCParticles[fInput_Tracks.McEntry->at(track_entry)]);
            if (charge < 0) fNegKaons.push_back(kaon);
            if (charge > 0) fPosKaons.push_back(kaon);
        }
        if (std::abs(fInput_Tracks.NSigmaPion->at(track_entry)) < Cuts::Track::AbsMax_PID_NSigma) {
            auto pion = std::make_shared<Charged>(track_entry, charge, xyz0, PxPyPzMVector{px0, py0, pz0, Const::MassPion});
            // if (IsMC()) pion->SetLinkedMc(fMCParticles[fInput_Tracks.McEntry->at(track_entry)]);
            if (charge < 0) fPiMinus.push_back(pion);
            if (charge > 0) fPiPlus.push_back(pion);
        }
    }
    // INFO("finished ProcessTracks() with n_protons = %zu, n_poskaons = %zu, n_piplus = %zu", fProtons.size(), fPosKaons.size(), fPiPlus.size());
    // INFO("finished ProcessTracks() with n_antiprotons = %zu, n_negkaons = %zu, n_piminus = %zu", fAntiProtons.size(), fNegKaons.size(),
    //  fPiMinus.size());
}

// ## V0s ZONE ## //

void Manager::FindV0s(int pdg_code_v0, int pdg_code_neg, int pdg_code_pos) {
    // choose tracks species to loop over //
    const auto& neg_vec{pdg_code_neg == PdgCode::AntiProton ? fAntiProtons : fPiMinus};
    const auto& pos_vec{pdg_code_pos == PdgCode::Proton ? fProtons : fPiPlus};
    // loop over all possible pairs of tracks //
    for (const auto& neg : neg_vec) {
        for (const auto& pos : pos_vec) {
            // sanity check //
            if (neg->Entry() == pos->Entry()) continue;
            // fit //
            Particle::Pair res{Vertexer::MinimizeDistanceHelixHelix(*neg, *pos, fPropagator)};
            auto V0 = std::make_shared<Neutral>(neg->Entry(), pos->Entry(), res);  // PENDING: indexing not correct..
            // apply cuts //
            if (!PassesV0Cuts(V0, pdg_code_v0)) continue;
            // store //
            Store(V0, pdg_code_v0);
        }  // end of loop over pos
    }  // end of loop over neg
}

bool Manager::PassesLambdaCuts(const std::shared_ptr<Neutral>& v0) const {
    /*
        INFO("m%f dbd%f z%f r%f dn%f dp%f pt%f et%f qt%f a%f cpv%f dpv%f", v0->Mass(), v0->DCAbtwDaughters(), v0->DecayZ(), v0->DecayRadius(),
             v0->DCANegWrtV0(), v0->DCAPosWrtV0(), v0->Pt(), v0->Eta(), v0->ArmenterosQt(), v0->ArmenterosAlpha(),
             v0->CPAwrt(fPropagator.PrimaryVertex()), v0->DCAwrt(fPropagator.PrimaryVertex()));
    */
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
    /*
        INFO("m%f dbd%f z%f r%f dn%f dp%f pt%f et%f qt%f a%f cpv%f dpv%f", v0->Mass(), v0->DCAbtwDaughters(), v0->DecayZ(), v0->DecayRadius(),
             v0->DCANegWrtV0(), v0->DCAPosWrtV0(), v0->Pt(), v0->Eta(), v0->ArmenterosQt(), v0->ArmenterosAlpha(),
             v0->CPAwrt(fPropagator.PrimaryVertex()), v0->DCAwrt(fPropagator.PrimaryVertex()));
    */
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

//
void Manager::Store(const std::shared_ptr<Neutral>& v0, int pdg_code_v0) {
    // fill temporary container //
    if (pdg_code_v0 == PdgCode::AntiLambda)
        fAntiLambdas.push_back(v0);
    else if (pdg_code_v0 == PdgCode::Lambda)
        fLambdas.push_back(v0);
    else if (pdg_code_v0 == PdgCode::KaonZeroShort)
        fNeutralKaons.push_back(v0);
    else {
        WARNING("Unknown V0 type");
        return;
    }
    // fill branches //
    fOutput_V0s.Index->push_back(0);  // INDEXING PENDING
    fOutput_V0s.PID->push_back(pdg_code_v0);
    fOutput_V0s.Neg_Entry->push_back(v0->NegIndex());
    fOutput_V0s.Pos_Entry->push_back(v0->PosIndex());
    fOutput_V0s.Xv->push_back(v0->DecayX());
    fOutput_V0s.Yv->push_back(v0->DecayY());
    fOutput_V0s.Zv->push_back(v0->DecayZ());
    fOutput_V0s.Px->push_back(v0->Px());
    fOutput_V0s.Py->push_back(v0->Py());
    fOutput_V0s.Pz->push_back(v0->Pz());
    fOutput_V0s.E->push_back(v0->Energy());
    fOutput_V0s.Neg_Xv->push_back(v0->NegVertex().X());
    fOutput_V0s.Neg_Yv->push_back(v0->NegVertex().Y());
    fOutput_V0s.Neg_Zv->push_back(v0->NegVertex().Z());
    fOutput_V0s.Neg_Px->push_back(v0->NegMomentum().X());
    fOutput_V0s.Neg_Py->push_back(v0->NegMomentum().Y());
    fOutput_V0s.Neg_Pz->push_back(v0->NegMomentum().Z());
    fOutput_V0s.Pos_Xv->push_back(v0->PosVertex().X());
    fOutput_V0s.Pos_Yv->push_back(v0->PosVertex().Y());
    fOutput_V0s.Pos_Zv->push_back(v0->PosVertex().Z());
    fOutput_V0s.Pos_Px->push_back(v0->PosMomentum().X());
    fOutput_V0s.Pos_Py->push_back(v0->PosMomentum().Y());
    fOutput_V0s.Pos_Pz->push_back(v0->PosMomentum().Z());
}

// ## END OF CYCLES ## //

// Fill output tree and clean transitory containers.
void Manager::EndOfEvent() {
    // fill tree //
    fOutputTree->Fill();
    // clear temporary containers //
    // if (IsMC()) fMCParticles.clear();
    fAntiProtons.clear();
    fProtons.clear();
    fNegKaons.clear();
    fPosKaons.clear();
    fPiMinus.clear();
    fPiPlus.clear();
    //
    fAntiLambdas.clear();
    fLambdas.clear();
    fNeutralKaons.clear();
    // clear output branches //
    fOutput_V0s.Clear();
}

//
void Manager::EndOfAnalysis() {
    fOutputTree->Write();

    fEventsTree->ResetBranchAddresses();
    fOutputTree->ResetBranchAddresses();
}

}  // namespace Tree2Secondaries::Analysis
