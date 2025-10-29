#include <filesystem>
#include <memory>

#include "App/Logger.hxx"
#include "DataFormats/PackedEvents.hxx"
#include "Finder/Cuts.hxx"
#include "Finder/Finder.hxx"
#include "Fit/ChannelA.hxx"
#include "Fit/ChannelD.hxx"
#include "Fit/Track.hxx"
#include "Fit/Utilities.hxx"
#include "Math/Constants.hxx"
#include "References/References2.hxx"

namespace Tree2Secondaries {

bool Finder::Initialize() {

    fInputChain_PackedEvents = std::make_unique<TChain>(Const::TreeName_PackedEvents.c_str());
    for (const auto& path : fSettings.PathInputFiles) {
        if (fInputChain_PackedEvents->Add(path.c_str()) == 0) {
            Logger::Error(__FUNCTION__, "Couldn't add TFile {}", path);
        }
    }
    if (!fInputChain_PackedEvents->GetEntries()) {
        Logger::Error(__FUNCTION__, "Couldn't manage to read any entry.");
        return false;
    }
    Logger::Info(__FUNCTION__, "TChain \"{}\" loaded successfully with {} trees and {} total entries.", fInputChain_PackedEvents->GetName(),
                 fInputChain_PackedEvents->GetNtrees(), fInputChain_PackedEvents->GetEntries());

    ReadInputBranches();

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

void Finder::ReadInputBranches() {
    // by default, turn off all branches
    fInputChain_PackedEvents->SetBranchStatus("*", false);
    // connect input branches to memory
    ReadBranches_Events();
    if (IsMC()) ReadBranches_Injected();
    // -- depending on reaction channels
    switch (GetReactionChannel()) {
        case EReactionChannel::A:
            ReadBranches_V0s(EParticle::AntiLambda, fInput_AntiLambdas);
            ReadBranches_V0s(EParticle::Lambda, fInput_Lambdas);
            ReadBranches_V0s(EParticle::KaonZeroShort, fInput_KaonsZeroShort);
            if (IsMC()) {
                ReadBranches_LinkedV0s(EParticle::AntiLambda, fInput_Linked_AntiLambdas);
                ReadBranches_LinkedV0s(EParticle::Lambda, fInput_Linked_Lambdas);
                ReadBranches_LinkedV0s(EParticle::KaonZeroShort, fInput_Linked_KaonsZeroShort);
            }
            break;
        case EReactionChannel::D:
            ReadBranches_V0s(EParticle::AntiLambda, fInput_AntiLambdas);
            ReadBranches_V0s(EParticle::Lambda, fInput_Lambdas);
            ReadBranches_Tracks(EParticle::NegKaon, fInput_NegKaons);
            ReadBranches_Tracks(EParticle::PosKaon, fInput_PosKaons);
            if (IsMC()) {
                ReadBranches_LinkedV0s(EParticle::AntiLambda, fInput_Linked_AntiLambdas);
                ReadBranches_LinkedV0s(EParticle::Lambda, fInput_Linked_Lambdas);
                ReadBranches_LinkedTracks(EParticle::NegKaon, fInput_Linked_NegKaons);
                ReadBranches_LinkedTracks(EParticle::PosKaon, fInput_Linked_PosKaons);
            }
            break;
        case EReactionChannel::E:
            ReadBranches_V0s(EParticle::AntiLambda, fInput_AntiLambdas);
            ReadBranches_V0s(EParticle::Lambda, fInput_Lambdas);
            ReadBranches_Tracks(EParticle::NegKaon, fInput_NegKaons);
            ReadBranches_Tracks(EParticle::PosKaon, fInput_PosKaons);
            ReadBranches_Tracks(EParticle::PiMinus, fInput_PiMinus);
            ReadBranches_Tracks(EParticle::PiPlus, fInput_PiPlus);
            if (IsMC()) {
                ReadBranches_LinkedV0s(EParticle::AntiLambda, fInput_Linked_AntiLambdas);
                ReadBranches_LinkedV0s(EParticle::Lambda, fInput_Linked_Lambdas);
                ReadBranches_LinkedTracks(EParticle::NegKaon, fInput_Linked_NegKaons);
                ReadBranches_LinkedTracks(EParticle::PosKaon, fInput_Linked_PosKaons);
                ReadBranches_LinkedTracks(EParticle::PiMinus, fInput_Linked_PiMinus);
                ReadBranches_LinkedTracks(EParticle::PiPlus, fInput_Linked_PiPlus);
            }
            break;
        case EReactionChannel::H:
            ReadBranches_Tracks(EParticle::PosKaon, fInput_PosKaons);
            ReadBranches_Tracks(EParticle::NegKaon, fInput_NegKaons);
            if (IsMC()) {
                ReadBranches_LinkedTracks(EParticle::NegKaon, fInput_Linked_NegKaons);
                ReadBranches_LinkedTracks(EParticle::PosKaon, fInput_Linked_PosKaons);
            }
            break;
    }  // end of switch statement
}

void Finder::ReadBranches_Events() { fInput_Event.ReadBranches_Event(fInputChain_PackedEvents.get(), IsMC()); }

void Finder::ReadBranches_Injected() { fInput_Injected.ReadBranches_SOV_Injected(fInputChain_PackedEvents.get(), true); }

void Finder::ReadBranches_V0s(EParticle pid, DF::Packed::V0s& df) {
    df.ReadBranches_PackedV0s(fInputChain_PackedEvents.get(), Const::Particle_Acronym[pid]);
}

void Finder::ReadBranches_Tracks(EParticle pid, DF::Packed::Tracks& df) {
    df.ReadBranches_PackedTracks(fInputChain_PackedEvents.get(), Const::Particle_Acronym[pid]);
}

void Finder::ReadBranches_LinkedV0s(EParticle pid, DF::Packed::LinkedV0s& df) {
    df.ReadBranches_LinkedV0s(fInputChain_PackedEvents.get(), Const::Particle_Acronym[pid]);
}

void Finder::ReadBranches_LinkedTracks(EParticle pid, DF::Packed::LinkedTracks& df) {
    df.ReadBranches_LinkedTracks(fInputChain_PackedEvents.get(), Const::Particle_Acronym[pid]);
}

// ## OUTPUT ZONE ## //

bool Finder::PrepareOutputFile() {

    const std::filesystem::path output_path(fSettings.PathOutputFile);
    if (output_path.has_parent_path()) std::filesystem::create_directories(output_path.parent_path());

    fOutputFile = std::unique_ptr<TFile>(TFile::Open(fSettings.PathOutputFile.c_str(), "RECREATE"));
    if (!fOutputFile) {
        Logger::Error(__FUNCTION__, "Couldn't create TFile {}", fSettings.PathOutputFile);
        return false;
    }

    return true;
}

bool Finder::PrepareOutputTree() {

    std::string tree_name{std::format("FoundChannel{}", static_cast<char>(fSettings.ReactionChannel))};

    fOutputTree = std::make_unique<TTree>(tree_name.c_str(), "");
    if (!fOutputTree) {
        Logger::Error(__FUNCTION__, "Couldn't create TTree \"{}\"", tree_name);
        return false;
    }

    return true;
}

void Finder::CreateOutputBranches(DF::Found::ChannelA& df) { df.CreateBranches_ChannelA(fOutputTree.get(), IsMC()); }

void Finder::CreateOutputBranches(DF::Found::ChannelD& df) { df.CreateBranches_ChannelD(fOutputTree.get(), IsMC()); }

void Finder::CreateOutputBranches(DF::Found::MC_ChannelA& df) { df.CreateBranches_MC_ChannelA(fOutputTree.get()); }

void Finder::CreateOutputBranches(DF::Found::MC_ChannelD& df) { df.CreateBranches_MC_ChannelD(fOutputTree.get()); }

void Finder::CreateCutFlowHistogram() {
    const int x_nbins{20};
    const float x_min{0.};
    const float x_max{20.};
    std::string hist_title{";Cut N;N Passed Cut"};
    fCutFlowHist = std::make_unique<TH1D>("CutFlow", hist_title.c_str(), x_nbins, x_min, x_max);
    fCutFlowHist_Anti = std::make_unique<TH1D>("CutFlow_Anti", hist_title.c_str(), x_nbins, x_min, x_max);
}

// ## OUTPUT / Injected ZONE ## //

bool Finder::Injected_PrepareOutputTree() {

    std::string tree_name{Const::TreeName_Injected};

    fOutputTree_Injected = std::make_unique<TTree>(tree_name.c_str(), "");
    if (!fOutputTree_Injected) {
        Logger::Error(__FUNCTION__, "Couldn't create TTree \"{}\"", tree_name);
        return false;
    }

    return true;
}

void Finder::Injected_CreateOutputBranches() { fOutput_Injected.CreateBranches_Flat_Injected(fOutputTree_Injected.get()); };

void Finder::Injected_FlattenAndStore() {

    double nucleon_mass{Const::Particle_Mass[Const::ReactionNucleonPID[fSettings.ReactionChannel]]};

    auto n_injected = static_cast<int>(fInput_Injected.ReactionID->size());
    for (int idx_inj{0}; idx_inj < n_injected; ++idx_inj) {
        // `Flat::State`
        float sexa_px{fInput_Injected.Px->at(idx_inj)};
        float sexa_py{fInput_Injected.Py->at(idx_inj)};
        float sexa_pz{fInput_Injected.Pz->at(idx_inj)};
        float sexa_energy{static_cast<float>(
            std::sqrt(fSettings.SexaquarkMass * fSettings.SexaquarkMass + sexa_px * sexa_px + sexa_py * sexa_py + sexa_pz * sexa_pz))};
        fOutput_Injected.Fill_Coordinates(fInput_Injected.X->at(idx_inj), fInput_Injected.Y->at(idx_inj), fInput_Injected.Z->at(idx_inj));
        fOutput_Injected.Fill_LorentzVector(sexa_px, sexa_py, sexa_pz, sexa_energy);
        // `Flat::Injected`
        // -- nucleon (`Flat::LorentzVector`)
        float nucleon_px{fInput_Injected.Nucleon.Px->at(idx_inj)};
        float nucleon_py{fInput_Injected.Nucleon.Py->at(idx_inj)};
        float nucleon_pz{fInput_Injected.Nucleon.Pz->at(idx_inj)};
        float nucleon_energy{static_cast<float>(  //
            std::sqrt(nucleon_mass * nucleon_mass + nucleon_px * nucleon_px + nucleon_py * nucleon_py + nucleon_pz * nucleon_pz))};
        fOutput_Injected.Nucleon.Fill_LorentzVector(nucleon_px, nucleon_py, nucleon_pz, nucleon_energy);
        // -- event properties
        fOutput_Injected.RunNumber = fInput_Event.RunNumber;
        fOutput_Injected.DirNumber = fInput_Event.DirNumber;
        fOutput_Injected.EventNumber = fInput_Event.EventNumber;
        // -- reaction id
        fOutput_Injected.ReactionID = fInput_Injected.ReactionID->at(idx_inj);
        fOutputTree_Injected->Fill();
    }
}

// ## Channel A ZONE ## //

void Finder::FindSexaquarks_ChannelA(bool anti_channel) {

    // determine properties based on anti-channel or not
    // -- v0a
    const DF::Packed::V0s* Packed_V0A{&fInput_AntiLambdas};
    const DF::Packed::LinkedV0s* MC_V0A{&fInput_Linked_AntiLambdas};
    EParticle pid_v0a_neg{EParticle::AntiProton};
    EParticle pid_v0a_pos{EParticle::PiPlus};
    if (anti_channel) {
        Packed_V0A = &fInput_Lambdas;
        MC_V0A = &fInput_Linked_Lambdas;
        pid_v0a_neg = EParticle::Proton;
        pid_v0a_pos = EParticle::PiMinus;
    }
    int charge_v0a_neg{Const::Particle_Charge[pid_v0a_neg]};
    int charge_v0a_pos{Const::Particle_Charge[pid_v0a_pos]};
    double mass_v0a_neg{Const::Particle_Mass[pid_v0a_neg]};
    double mass_v0a_pos{Const::Particle_Mass[pid_v0a_pos]};
    // -- v0b
    const DF::Packed::V0s* Packed_V0B{&fInput_KaonsZeroShort};
    const DF::Packed::LinkedV0s* MC_V0B{&fInput_Linked_KaonsZeroShort};
    EParticle pid_v0b_neg{EParticle::PiMinus};
    EParticle pid_v0b_pos{EParticle::PiPlus};
    int charge_v0b_neg{Const::Particle_Charge[pid_v0b_neg]};
    int charge_v0b_pos{Const::Particle_Charge[pid_v0b_pos]};
    double mass_v0b_neg{Const::Particle_Mass[pid_v0b_neg]};
    double mass_v0b_pos{Const::Particle_Mass[pid_v0b_pos]};
    // -- cut flow hist
    TH1D* hist{anti_channel ? fCutFlowHist_Anti.get() : fCutFlowHist.get()};

    // loop over all possible pairs of (anti)lambda + K0S //
    auto n_v0a = static_cast<int>(Packed_V0A->Entry->size());
    auto n_v0b = static_cast<int>(Packed_V0B->Entry->size());
    for (int v0a_entry{0}; v0a_entry < n_v0a; ++v0a_entry) {

        // unpack (anti)lambda //
        Fit::Track v0a_neg{Fit::UnpackTrack(Packed_V0A->Neg, v0a_entry, charge_v0a_neg, mass_v0a_neg)};
        Fit::Track v0a_pos{Fit::UnpackTrack(Packed_V0A->Pos, v0a_entry, charge_v0a_pos, mass_v0a_pos)};
        Fit::V0 v0a{Fit::UnpackV0(*Packed_V0A, v0a_entry, v0a_neg, v0a_pos)};

        for (int v0b_entry{0}; v0b_entry < n_v0b; ++v0b_entry) {

            // unpack K0S //
            Fit::Track v0b_neg{Fit::UnpackTrack(Packed_V0B->Neg, v0b_entry, charge_v0b_neg, mass_v0b_neg)};
            Fit::Track v0b_pos{Fit::UnpackTrack(Packed_V0B->Pos, v0b_entry, charge_v0b_pos, mass_v0b_pos)};
            Fit::V0 v0b{Fit::UnpackV0(*Packed_V0B, v0b_entry, v0b_neg, v0b_pos)};

            // sanity check //
            if (v0a.Neg.Index == v0b.Neg.Index || v0a.Neg.Index == v0b.Pos.Index || v0a.Pos.Index == v0b.Neg.Index ||
                v0a.Pos.Index == v0b.Pos.Index) {
                continue;
            }

            // fit //
            Fit::ChannelA sexa{v0a, v0b};
            sexa.DoFit(fInput_Event.MagneticField);

#ifdef T2S_DEBUG
            Logger::Debug(__FUNCTION__, "idx(v0a,neg,pos)={},{},{}", v0a.Entry, v0a.Neg.Index, v0a.Pos.Index);
            Logger::Debug(__FUNCTION__, ";x,y,z={},{},{}", v0a.X(), v0a.Y(), v0a.Z());
            Logger::Debug(__FUNCTION__, ";px,py,pz={},{},{}", v0a.Px(), v0a.Py(), v0a.Pz());
            Logger::Debug(__FUNCTION__, ";mass={}", v0a.Mass());

            Logger::Debug(__FUNCTION__, "idx(v0b,neg,pos)={},{},{}", v0b.Entry, v0b.Neg.Index, v0b.Pos.Index);
            Logger::Debug(__FUNCTION__, ";x,y,z={},{},{}", v0b.X(), v0b.Y(), v0b.Z());
            Logger::Debug(__FUNCTION__, ";px,py,pz={},{},{}", v0b.Px(), v0b.Py(), v0b.Pz());
            Logger::Debug(__FUNCTION__, ";mass={}", v0b.Mass());

            Logger::Debug(__FUNCTION__, "x,y,z={},{},{}", sexa.X(), sexa.Y(), sexa.Z());
            // Logger::Debug(__FUNCTION__, ";x,y,z(v0a)={},{},{}", sexa.V0A_PCA_XYZ()[0], sexa.V0A_PCA_XYZ()[1], sexa.V0A_PCA_XYZ()[2]);
            // Logger::Debug(__FUNCTION__, ";x,y,z(v0b)={},{},{}", sexa.V0B_PCA_XYZ()[0], sexa.V0B_PCA_XYZ()[1], sexa.V0B_PCA_XYZ()[2]);
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
            Logger::Debug(__FUNCTION__, ";cpa_pv={}", sexa.CPA_Point(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z));

            // PENDING: to rewrite... chore
            // if (IsMC()) {
            // Ref::PackedV0 mc_v0a_debug{*MC_V0A, v0a.Entry};
            // Ref::PackedV0 mc_v0b_debug{*MC_V0B, v0b.Entry};
            // Link::ChannelA mc_sexa_debug{fInput_Injected, fSettings.SexaquarkMass, mc_v0a_debug, mc_v0b_debug};
            // Logger::Debug(__FUNCTION__, "sexa(is_signal,reaction_id,is_hybrid)={},{},{}", mc_sexa_debug.IsSignal, mc_sexa_debug.ReactionID,
            //   mc_sexa_debug.IsHybrid);
            // Logger::Debug(__FUNCTION__, ";v0a(entry,pdg_code)={},{}", mc_sexa_debug.V0A.Entry, mc_sexa_debug.V0A.PdgCode);
            // Logger::Debug(__FUNCTION__, ";v0a(is_signal,reaction_id,is_hybrid)={},{},{}", mc_sexa_debug.V0A.IsSignal,
            //   mc_sexa_debug.V0A.ReactionID, mc_sexa_debug.V0A.V0_IsHybrid());
            // Logger::Debug(__FUNCTION__, ";v0b(entry,pdg_code)={},{}", mc_sexa_debug.V0B.Entry, mc_sexa_debug.V0B.PdgCode);
            // Logger::Debug(__FUNCTION__, ";v0b(is_signal,reaction_id,is_hybrid)={},{},{}", mc_sexa_debug.V0B.IsSignal,
            //   mc_sexa_debug.V0B.ReactionID, mc_sexa_debug.V0B.V0_IsHybrid());
            // }
#endif

            // apply cuts //
            if (!PassesCuts(sexa, hist)) continue;

            // store //
            Store(sexa, anti_channel);

            if (IsMC()) {
                Ref::PackedV0 mc_v0a{MC_V0A, v0a.Entry};
                Ref::PackedV0 mc_v0b{MC_V0B, v0b.Entry};
                Ref::Injected mc_injected{
                    .source = &fInput_Injected, .mass = fSettings.SexaquarkMass, .nucleon_pid = Const::ReactionNucleonPID[fSettings.ReactionChannel]};
                Ref::ChannelA mc_sexa{mc_injected, mc_v0a, mc_v0b};
                StoreMC(mc_sexa);
            }
            fOutputTree->Fill();
        }
    }
}

bool Finder::PassesCuts(const Fit::ChannelA& sexa, TH1D* cut_flow_hist) const {

    cut_flow_hist->Fill(0.);
    if (sexa.Radius2D() < Cuts::ChannelA::Min_Radius2D || sexa.Radius2D() > Cuts::ChannelA::Max_Radius2D) return false;
    cut_flow_hist->Fill(1.);
    if (sexa.DecayLength_V0A() > Cuts::ChannelA::Max_DecayLengthLa) return false;
    cut_flow_hist->Fill(2.);
    if (sexa.DecayLength_V0B() > Cuts::ChannelA::Max_DecayLengthK0) return false;
    cut_flow_hist->Fill(3.);
    if (sexa.AbsRapidity_MinusNucleon() > Cuts::ChannelA::AbsMax_Rapidity) return false;  // PENDING: kinematic cut, affected by Fermi motion
    cut_flow_hist->Fill(4.);
    if (sexa.Mass_MinusNucleon() < Cuts::ChannelA::Min_MassMinusNucleon || sexa.Mass_MinusNucleon() > Cuts::ChannelA::Max_MassMinusNucleon) {
        return false;  // PENDING: kinematic cut, affected by Fermi motion
    }
    cut_flow_hist->Fill(5.);
    if (sexa.CPA_Point(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) < Cuts::ChannelA::Min_CPAwrtPV ||
        sexa.CPA_Point(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) > Cuts::ChannelA::Max_CPAwrtPV) {
        return false;  // PENDING: kinematic cut, affected by Fermi motion
    }
    cut_flow_hist->Fill(6.);
    if (sexa.DCA_V0ANeg_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelA::Max_DCALaNegSV) return false;
    cut_flow_hist->Fill(7.);
    if (sexa.DCA_V0APos_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelA::Max_DCALaPosSV) return false;
    cut_flow_hist->Fill(8.);
    if (sexa.DCA_V0BNeg_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelA::Max_DCAK0NegSV) return false;
    cut_flow_hist->Fill(9.);
    if (sexa.DCA_V0BPos_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelA::Max_DCAK0PosSV) return false;
    cut_flow_hist->Fill(10.);
    if (sexa.DCA_V0A_wrt_SV() > Cuts::ChannelA::Max_DCALaSV) return false;
    cut_flow_hist->Fill(11.);
    if (sexa.DCA_V0B_wrt_SV() > Cuts::ChannelA::Max_DCAK0SV) return false;
    cut_flow_hist->Fill(12.);
    if (sexa.DCA_btw_V0s() > Cuts::ChannelA::Max_DCAbtwV0s) return false;
    cut_flow_hist->Fill(13.);

    return true;
}

void Finder::Store(const Fit::ChannelA& sexa, bool anti_channel) {
    // `DF::Sexaquark` //
    // -- `DF::Flat::State`
    fOutput_ChannelA.Fill_State(static_cast<float>(sexa.X()), static_cast<float>(sexa.Y()), static_cast<float>(sexa.Z()),
                                static_cast<float>(sexa.Px()), static_cast<float>(sexa.Py()), static_cast<float>(sexa.Pz()),
                                static_cast<float>(sexa.E()));
    // -- event properties
    fOutput_ChannelA.RunNumber = fInput_Event.RunNumber;
    fOutput_ChannelA.DirNumber = fInput_Event.DirNumber;
    if (!IsMC()) fOutput_ChannelA.DirNumberB = fInput_Event.DirNumberB;
    fOutput_ChannelA.EventNumber = fInput_Event.EventNumber;
    fOutput_ChannelA.MagneticField = fInput_Event.MagneticField;
    // -- primary vertex
    fOutput_ChannelA.PV = fInput_Event.PV;
    // -- fit info
    fOutput_ChannelA.Chi2NDF = static_cast<float>(sexa.Chi2NDF());
    // -- extra info
    fOutput_ChannelA.E_MinusNucleon = static_cast<float>(sexa.E_MinusNucleon());
    fOutput_ChannelA.AntiChannel = anti_channel;

    // `DF::ChannelA` //
    // V0A (`DF::Flat::State`)
    fOutput_ChannelA.V0A.Fill_State(                                                                              //
        static_cast<float>(sexa.V0A.X()), static_cast<float>(sexa.V0A.Y()), static_cast<float>(sexa.V0A.Z()),     //
        static_cast<float>(sexa.V0A.Px()), static_cast<float>(sexa.V0A.Py()), static_cast<float>(sexa.V0A.Pz()),  //
        static_cast<float>(sexa.V0A.E()));
    // -- V0A neg (`DF::Flat::Track`)
    fOutput_ChannelA.V0A.Neg_atV0.Fill_State(static_cast<float>(sexa.V0A.Neg.X()), static_cast<float>(sexa.V0A.Neg.Y()),
                                             static_cast<float>(sexa.V0A.Neg.Z()), static_cast<float>(sexa.V0A.Neg.Px()),
                                             static_cast<float>(sexa.V0A.Neg.Py()), static_cast<float>(sexa.V0A.Neg.Pz()),
                                             static_cast<float>(sexa.V0A.Neg.E()));
    fOutput_ChannelA.V0A.Neg_atV0.Index = sexa.V0A.Neg.Index;
    // -- V0A pos (`DF::Flat::Track`)
    fOutput_ChannelA.V0A.Pos_atV0.Fill_State(static_cast<float>(sexa.V0A.Pos.X()), static_cast<float>(sexa.V0A.Pos.Y()),
                                             static_cast<float>(sexa.V0A.Pos.Z()), static_cast<float>(sexa.V0A.Pos.Px()),
                                             static_cast<float>(sexa.V0A.Pos.Py()), static_cast<float>(sexa.V0A.Pos.Pz()),
                                             static_cast<float>(sexa.V0A.Pos.E()));
    fOutput_ChannelA.V0A.Pos_atV0.Index = sexa.V0A.Pos.Index;
    // -- entry
    fOutput_ChannelA.V0A.Entry = sexa.V0A.Entry;
    // V0A @ PCA (`DF::Flat::Coordinates`)
    fOutput_ChannelA.V0A_atPCA.X = static_cast<float>(sexa.V0A.Pos.X());
    fOutput_ChannelA.V0A_atPCA.Y = static_cast<float>(sexa.V0A.Pos.Y());
    fOutput_ChannelA.V0A_atPCA.Z = static_cast<float>(sexa.V0A.Pos.Z());

    // V0B (`DF::Flat::State`)
    fOutput_ChannelA.V0B.Fill_State(                                                                              //
        static_cast<float>(sexa.V0B.X()), static_cast<float>(sexa.V0B.Y()), static_cast<float>(sexa.V0B.Z()),     //
        static_cast<float>(sexa.V0B.Px()), static_cast<float>(sexa.V0B.Py()), static_cast<float>(sexa.V0B.Pz()),  //
        static_cast<float>(sexa.V0B.E()));
    // -- V0B neg (`DF::Flat::Track`)
    fOutput_ChannelA.V0B.Neg_atV0.Fill_State(static_cast<float>(sexa.V0B.Neg.X()), static_cast<float>(sexa.V0B.Neg.Y()),
                                             static_cast<float>(sexa.V0B.Neg.Z()), static_cast<float>(sexa.V0B.Neg.Px()),
                                             static_cast<float>(sexa.V0B.Neg.Py()), static_cast<float>(sexa.V0B.Neg.Pz()),
                                             static_cast<float>(sexa.V0B.Neg.E()));
    fOutput_ChannelA.V0B.Neg_atV0.Index = sexa.V0B.Neg.Index;
    // -- V0B pos (`DF::Flat::Track`)
    fOutput_ChannelA.V0B.Pos_atV0.Fill_State(static_cast<float>(sexa.V0B.Pos.X()), static_cast<float>(sexa.V0B.Pos.Y()),
                                             static_cast<float>(sexa.V0B.Pos.Z()), static_cast<float>(sexa.V0B.Pos.Px()),
                                             static_cast<float>(sexa.V0B.Pos.Py()), static_cast<float>(sexa.V0B.Pos.Pz()),
                                             static_cast<float>(sexa.V0B.Pos.E()));
    fOutput_ChannelA.V0B.Pos_atV0.Index = sexa.V0B.Pos.Index;
    // -- entry
    fOutput_ChannelA.V0B.Entry = sexa.V0B.Entry;
    // V0B @ PCA (`DF::Flat::Coordinates`)
    fOutput_ChannelA.V0B_atPCA.X = static_cast<float>(sexa.V0B.Pos.X());
    fOutput_ChannelA.V0B_atPCA.Y = static_cast<float>(sexa.V0B.Pos.Y());
    fOutput_ChannelA.V0B_atPCA.Z = static_cast<float>(sexa.V0B.Pos.Z());
}

void Finder::StoreMC(const Ref::ChannelA& sexa) {
    // `Found::MC_Sexaquark` //
    fOutput_MC_ChannelA.Before.Fill_LorentzVector(sexa.BeforePx(), sexa.BeforePy(), sexa.BeforePz(), static_cast<float>(sexa.BeforeE()));
    fOutput_MC_ChannelA.After.Fill_LorentzVector(static_cast<float>(sexa.AsChannelA_AfterPx()), static_cast<float>(sexa.AsChannelA_AfterPy()),  //
                                                 static_cast<float>(sexa.AsChannelA_AfterPz()), static_cast<float>(sexa.AsChannelA_AfterE()));
    fOutput_MC_ChannelA.Nucleon.Fill_LorentzVector(sexa.NucleonPx(), sexa.NucleonPy(), sexa.NucleonPz(), static_cast<float>(sexa.NucleonE()));
    // -- secondary vertex
    fOutput_MC_ChannelA.SV.Fill_Coordinates(sexa.SV_X(), sexa.SV_Y(), sexa.SV_Z());
    // -- event properties
    fOutput_MC_ChannelA.PV = fInput_Event.MC_PV;
    // -- reaction id + flags
    fOutput_MC_ChannelA.ReactionID = sexa.AsChannelA_ReactionID();
    fOutput_MC_ChannelA.IsSignal = sexa.AsChannelA_IsSignal();
    fOutput_MC_ChannelA.IsHybrid = sexa.AsChannelA_IsHybrid();

    // `Found::MC_ChannelA` //
    // V0A (`Flat::MCInfo_V0`)
    // -- `Flat::MCInfo_LV_Mother`
    //    -- `Flat::MCInfo`
    fOutput_MC_ChannelA.V0A.Entry = sexa.v0a.Entry();
    fOutput_MC_ChannelA.V0A.PdgCode = sexa.v0a.PdgCode();
    fOutput_MC_ChannelA.V0A.ReactionID = sexa.v0a.ReactionID();
    fOutput_MC_ChannelA.V0A.IsTrue = sexa.v0a.IsTrue();
    fOutput_MC_ChannelA.V0A.IsSignal = sexa.v0a.IsSignal();
    fOutput_MC_ChannelA.V0A.IsSecondary = sexa.v0a.IsSecondary();
    //    -- `Flat::LorentzVector`
    fOutput_MC_ChannelA.V0A.Fill_LorentzVector(sexa.v0a.Px(), sexa.v0a.Py(), sexa.v0a.Pz(), sexa.v0a.Energy());
    //    -- mother info
    fOutput_MC_ChannelA.V0A.Mother_Entry = sexa.v0a.Mother_Entry();
    fOutput_MC_ChannelA.V0A.Mother_PdgCode = sexa.v0a.Mother_PdgCode();
    // -- V0A neg (`Flat::MCInfo_PxPyPz`)
    //    -- `Flat::MCInfo`
    fOutput_MC_ChannelA.V0A.Neg.Entry = sexa.v0a.Neg_Entry();
    fOutput_MC_ChannelA.V0A.Neg.PdgCode = sexa.v0a.Neg_PdgCode();
    fOutput_MC_ChannelA.V0A.Neg.ReactionID = sexa.v0a.Neg_ReactionID();
    fOutput_MC_ChannelA.V0A.Neg.IsTrue = sexa.v0a.Neg_IsTrue();
    fOutput_MC_ChannelA.V0A.Neg.IsSignal = sexa.v0a.Neg_IsSignal();
    fOutput_MC_ChannelA.V0A.Neg.IsSecondary = sexa.v0a.Neg_IsSecondary();
    //    -- `Flat::PxPyPz`
    fOutput_MC_ChannelA.V0A.Neg.Fill_PxPyPz(sexa.v0a.Neg_Px(), sexa.v0a.Neg_Py(), sexa.v0a.Neg_Pz());
    // -- V0A pos (`Flat::MCInfo_PxPyPz`)
    //    -- `Flat::MCInfo`
    fOutput_MC_ChannelA.V0A.Pos.Entry = sexa.v0a.Pos_Entry();
    fOutput_MC_ChannelA.V0A.Pos.PdgCode = sexa.v0a.Pos_PdgCode();
    fOutput_MC_ChannelA.V0A.Pos.ReactionID = sexa.v0a.Pos_ReactionID();
    fOutput_MC_ChannelA.V0A.Pos.IsTrue = sexa.v0a.Pos_IsTrue();
    fOutput_MC_ChannelA.V0A.Pos.IsSignal = sexa.v0a.Pos_IsSignal();
    fOutput_MC_ChannelA.V0A.Pos.IsSecondary = sexa.v0a.Pos_IsSecondary();
    //    -- `Flat::PxPyPz`
    fOutput_MC_ChannelA.V0A.Pos.Fill_PxPyPz(sexa.v0a.Pos_Px(), sexa.v0a.Pos_Py(), sexa.v0a.Pos_Pz());
    // -- V0A @ decay (`Flat::Coordinates`)
    fOutput_MC_ChannelA.V0A.AtDecay.Fill_Coordinates(sexa.v0a.DecayX(), sexa.v0a.DecayY(), sexa.v0a.DecayZ());
    // -- hybrid flag
    fOutput_MC_ChannelA.V0A.IsHybrid = sexa.v0a.IsHybrid();

    // V0B (`Flat::MCInfo_V0`)
    // -- `Flat::MCInfo_LV_Mother`
    //    -- `Flat::MCInfo`
    fOutput_MC_ChannelA.V0B.Entry = sexa.v0b.Entry();
    fOutput_MC_ChannelA.V0B.PdgCode = sexa.v0b.PdgCode();
    fOutput_MC_ChannelA.V0B.ReactionID = sexa.v0b.ReactionID();
    fOutput_MC_ChannelA.V0B.IsTrue = sexa.v0b.IsTrue();
    fOutput_MC_ChannelA.V0B.IsSignal = sexa.v0b.IsSignal();
    fOutput_MC_ChannelA.V0B.IsSecondary = sexa.v0b.IsSecondary();
    //    -- `Flat::LorentzVector`
    fOutput_MC_ChannelA.V0B.Fill_LorentzVector(sexa.v0b.Px(), sexa.v0b.Py(), sexa.v0b.Pz(), sexa.v0b.Energy());
    //    -- mother info
    fOutput_MC_ChannelA.V0B.Mother_Entry = sexa.v0b.Mother_Entry();
    fOutput_MC_ChannelA.V0B.Mother_PdgCode = sexa.v0b.Mother_PdgCode();
    // -- V0B neg (`Flat::MCInfo_PxPyPz`)
    //    -- `Flat::MCInfo`
    fOutput_MC_ChannelA.V0B.Neg.Entry = sexa.v0b.Neg_Entry();
    fOutput_MC_ChannelA.V0B.Neg.PdgCode = sexa.v0b.Neg_PdgCode();
    fOutput_MC_ChannelA.V0B.Neg.ReactionID = sexa.v0b.Neg_ReactionID();
    fOutput_MC_ChannelA.V0B.Neg.IsTrue = sexa.v0b.Neg_IsTrue();
    fOutput_MC_ChannelA.V0B.Neg.IsSignal = sexa.v0b.Neg_IsSignal();
    fOutput_MC_ChannelA.V0B.Neg.IsSecondary = sexa.v0b.Neg_IsSecondary();
    //    -- `Flat::PxPyPz`
    fOutput_MC_ChannelA.V0B.Neg.Fill_PxPyPz(sexa.v0b.Neg_Px(), sexa.v0b.Neg_Py(), sexa.v0b.Neg_Pz());
    // -- V0B pos (`Flat::MCInfo_PxPyPz`)
    //    -- `Flat::MCInfo`
    fOutput_MC_ChannelA.V0B.Pos.Entry = sexa.v0b.Pos_Entry();
    fOutput_MC_ChannelA.V0B.Pos.PdgCode = sexa.v0b.Pos_PdgCode();
    fOutput_MC_ChannelA.V0B.Pos.ReactionID = sexa.v0b.Pos_ReactionID();
    fOutput_MC_ChannelA.V0B.Pos.IsTrue = sexa.v0b.Pos_IsTrue();
    fOutput_MC_ChannelA.V0B.Pos.IsSignal = sexa.v0b.Pos_IsSignal();
    fOutput_MC_ChannelA.V0B.Pos.IsSecondary = sexa.v0b.Pos_IsSecondary();
    //    -- `Flat::PxPyPz`
    fOutput_MC_ChannelA.V0B.Pos.Fill_PxPyPz(sexa.v0b.Pos_Px(), sexa.v0b.Pos_Py(), sexa.v0b.Pos_Pz());
    // -- V0B @ decay (`Flat::Coordinates`)
    fOutput_MC_ChannelA.V0B.AtDecay.Fill_Coordinates(sexa.v0b.DecayX(), sexa.v0b.DecayY(), sexa.v0b.DecayZ());
    // -- hybrid flag
    fOutput_MC_ChannelA.V0B.IsHybrid = sexa.v0b.IsHybrid();
}

// ## Channel D ZONE ## //

void Finder::FindSexaquarks_ChannelD(bool anti_channel) {

    // determine properties based on anti-channel or not
    // -- v0
    const DF::Packed::V0s* Packed_V0s{&fInput_AntiLambdas};
    const DF::Packed::LinkedV0s* MC_V0s{&fInput_Linked_AntiLambdas};
    EParticle pid_lambda_neg{EParticle::AntiProton};
    EParticle pid_lambda_pos{EParticle::PiPlus};
    if (anti_channel) {
        Packed_V0s = &fInput_Lambdas;
        MC_V0s = &fInput_Linked_Lambdas;
        pid_lambda_neg = EParticle::PiMinus;
        pid_lambda_pos = EParticle::Proton;
    }
    int charge_lambda_neg{Const::Particle_Charge[pid_lambda_neg]};
    int charge_lambda_pos{Const::Particle_Charge[pid_lambda_pos]};
    double mass_lambda_neg{Const::Particle_Mass[pid_lambda_neg]};
    double mass_lambda_pos{Const::Particle_Mass[pid_lambda_neg]};
    // -- kaon
    const DF::Packed::Tracks* Packed_Kaons{&fInput_PosKaons};
    const DF::Packed::LinkedTracks* MC_Kaons{&fInput_Linked_PosKaons};
    EParticle pid_kaon{EParticle::PosKaon};
    if (anti_channel) {
        Packed_Kaons = &fInput_NegKaons;
        MC_Kaons = &fInput_Linked_NegKaons;
        pid_kaon = EParticle::NegKaon;
    }
    double mass_kaon{Const::Particle_Mass[pid_kaon]};
    int charge_kaon{Const::Particle_Charge[pid_kaon]};
    // -- cut flow hist
    TH1D* hist{anti_channel ? fCutFlowHist_Anti.get() : fCutFlowHist.get()};

    // loop over all possible pairs of (anti)lambda + (pos/neg)kaon //
    auto n_v0 = static_cast<int>(Packed_V0s->Entry->size());
    auto n_kaons = static_cast<int>(Packed_Kaons->Index->size());
    for (int idx_v0{0}; idx_v0 < n_v0; ++idx_v0) {

        // unpack (anti)lambda //
        Fit::Track v0_neg{Fit::UnpackTrack(Packed_V0s->Neg, idx_v0, charge_lambda_neg, mass_lambda_neg)};
        Fit::Track v0_pos{Fit::UnpackTrack(Packed_V0s->Pos, idx_v0, charge_lambda_pos, mass_lambda_pos)};
        Fit::V0 v0{Fit::UnpackV0(*Packed_V0s, idx_v0, v0_neg, v0_pos)};

        for (int idx_kaon{0}; idx_kaon < n_kaons; ++idx_kaon) {

            // unpack (pos/neg)kaon //
            Fit::Track kaon{Fit::UnpackTrack(*Packed_Kaons, idx_kaon, charge_kaon, mass_kaon)};

            // sanity check //
            if (v0.Neg.Index == kaon.Index || v0.Pos.Index == kaon.Index) continue;

            // fit //
            Fit::ChannelD sexa{v0, kaon};
            sexa.DoFit(fInput_Event.MagneticField);

            // apply cuts //
            if (!PassesCuts(sexa, hist)) continue;
#ifdef T2S_DEBUG
            Logger::Debug(__FUNCTION__, "idx(v0,neg,pos,kaon)={},{},{},{}", sexa.V0.Entry, sexa.V0.Neg.Index, sexa.V0.Pos.Index, sexa.Kaon.Index);
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
            Logger::Debug(__FUNCTION__, ";cpa_pv={}", sexa.CPA_Point(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z));
#endif

            // store //
            Store(sexa, anti_channel);
            if (IsMC()) {
                Ref::PackedV0 mc_v0{MC_V0s, v0.Entry};
                Ref::PackedBachelor mc_kaon{MC_Kaons, kaon.Index};
                Ref::Injected mc_injected{
                    .source = &fInput_Injected, .mass = fSettings.SexaquarkMass, .nucleon_pid = Const::ReactionNucleonPID[fSettings.ReactionChannel]};
                Ref::ChannelD mc_sexa{mc_injected, mc_v0, mc_kaon};
                StoreMC(mc_sexa);
            }
            fOutputTree->Fill();
        }
    }
}

bool Finder::PassesCuts(const Fit::ChannelD& sexa, TH1D* cut_flow_hist) const {

    cut_flow_hist->Fill(0.);
    if (sexa.Radius2D() < Cuts::ChannelD::Min_Radius2D || sexa.Radius2D() > Cuts::ChannelD::Max_Radius2D) return false;
    cut_flow_hist->Fill(1.);
    if (sexa.AbsRapidity_MinusNucleon() > Cuts::ChannelD::AbsMax_Rapidity) return false;  // PENDING: kinematics, affected by Fermi motion
    cut_flow_hist->Fill(2.);
    if (sexa.CPA_Point(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) < Cuts::ChannelD::Min_CPAwrtPV ||
        sexa.CPA_Point(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) > Cuts::ChannelD::Max_CPAwrtPV) {
        return false;  // PENDING: kinematics, affected by Fermi motion
    }
    cut_flow_hist->Fill(3.);
    if (sexa.DCA_V0_wrt_SV() > Cuts::ChannelD::Max_DCALaSV) return false;
    cut_flow_hist->Fill(4.);
    if (sexa.DCA_Kaon_wrt_SV() > Cuts::ChannelD::Max_DCAKaSV) return false;
    cut_flow_hist->Fill(5.);
    if (sexa.DCA_V0Neg_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelD::Max_DCALaNegSV) return false;
    cut_flow_hist->Fill(6.);
    if (sexa.DCA_V0Pos_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelD::Max_DCALaPosSV) return false;
    cut_flow_hist->Fill(7.);
    if (sexa.DCA_btw_V0_Kaon() > Cuts::ChannelD::Max_DCAKaLa) return false;
    cut_flow_hist->Fill(8.);
    if (sexa.DCA_btw_V0Neg_Kaon(fInput_Event.MagneticField) > Cuts::ChannelD::Max_DCALaNegKa) return false;
    cut_flow_hist->Fill(9.);
    if (sexa.DCA_btw_V0Pos_Kaon(fInput_Event.MagneticField) > Cuts::ChannelD::Max_DCALaPosKa) return false;
    cut_flow_hist->Fill(10.);

    return true;
}

void Finder::Store(const Fit::ChannelD& sexa, bool anti_channel) {
    // `Sexaquark` //
    // -- `Flat::State`
    fOutput_ChannelD.Fill_State(static_cast<float>(sexa.X()), static_cast<float>(sexa.Y()), static_cast<float>(sexa.Z()),  //
                                static_cast<float>(sexa.Px()), static_cast<float>(sexa.Py()), static_cast<float>(sexa.Pz()),
                                static_cast<float>(sexa.E()));
    // -- event properties
    fOutput_ChannelD.PV = fInput_Event.PV;
    fOutput_ChannelD.RunNumber = fInput_Event.RunNumber;
    fOutput_ChannelD.DirNumber = fInput_Event.DirNumber;
    if (!IsMC()) fOutput_ChannelD.DirNumberB = fInput_Event.DirNumberB;
    fOutput_ChannelD.EventNumber = fInput_Event.EventNumber;
    fOutput_ChannelD.MagneticField = fInput_Event.MagneticField;
    // -- fit info
    fOutput_ChannelD.Chi2NDF = static_cast<float>(sexa.Chi2NDF());
    // -- extra info
    fOutput_ChannelD.E_MinusNucleon = static_cast<float>(sexa.E_MinusNucleon());
    fOutput_ChannelD.AntiChannel = anti_channel;

    // `ChannelD` //
    // V0 (`Flat::V0`)
    // -- `Flat::State`
    fOutput_ChannelD.V0.Fill_State(                                                                            //
        static_cast<float>(sexa.V0.X()), static_cast<float>(sexa.V0.Y()), static_cast<float>(sexa.V0.Z()),     //
        static_cast<float>(sexa.V0.Px()), static_cast<float>(sexa.V0.Py()), static_cast<float>(sexa.V0.Pz()),  //
        static_cast<float>(sexa.V0.E()));
    // -- Neg @ V0 (`Flat::Track`)
    fOutput_ChannelD.V0.Neg_atV0.Fill_Coordinates(  //
        static_cast<float>(sexa.V0.Neg.X()), static_cast<float>(sexa.V0.Neg.Y()), static_cast<float>(sexa.V0.Neg.Z()));
    fOutput_ChannelD.V0.Neg_atV0.Fill_PxPyPz(static_cast<float>(sexa.V0.Neg.Px()), static_cast<float>(sexa.V0.Neg.Py()),
                                             static_cast<float>(sexa.V0.Neg.Pz()));
    fOutput_ChannelD.V0.Neg_atV0.Index = sexa.V0.Neg.Index;
    // -- Pos @ V0 (`Flat::Track`)
    fOutput_ChannelD.V0.Pos_atV0.Fill_Coordinates(  //
        static_cast<float>(sexa.V0.Pos.X()), static_cast<float>(sexa.V0.Pos.Y()), static_cast<float>(sexa.V0.Pos.Z()));
    fOutput_ChannelD.V0.Pos_atV0.Fill_PxPyPz(static_cast<float>(sexa.V0.Pos.Px()), static_cast<float>(sexa.V0.Pos.Py()),
                                             static_cast<float>(sexa.V0.Pos.Pz()));
    fOutput_ChannelD.V0.Pos_atV0.Index = sexa.V0.Pos.Index;
    // -- entry
    fOutput_ChannelD.V0.Entry = sexa.V0.Entry;

    // Kaon (`Flat::Bachelor`)
    fOutput_ChannelD.Kaon.Fill_State(                                                                                //
        static_cast<float>(sexa.Kaon.X()), static_cast<float>(sexa.Kaon.Y()), static_cast<float>(sexa.Kaon.Z()),     //
        static_cast<float>(sexa.Kaon.Px()), static_cast<float>(sexa.Kaon.Py()), static_cast<float>(sexa.Kaon.Pz()),  //
        static_cast<float>(sexa.Kaon.E()));
    fOutput_ChannelD.Kaon.Index = sexa.Kaon.Index;

    // V0 @ PCA w.r.t. SV (`Flat::Coordinates`)
    fOutput_ChannelD.V0_atPCA.Fill_Coordinates(static_cast<float>(sexa.V0_PCA_XYZ()[0]), static_cast<float>(sexa.V0_PCA_XYZ()[1]),
                                               static_cast<float>(sexa.V0_PCA_XYZ()[2]));
    // Kaon @ PCA w.r.t. SV (`Flat::Coordinates`)
    fOutput_ChannelD.Kaon_atPCA.Fill_Coordinates(static_cast<float>(sexa.Kaon_PCA_XYZ()[0]), static_cast<float>(sexa.Kaon_PCA_XYZ()[1]),
                                                 static_cast<float>(sexa.Kaon_PCA_XYZ()[2]));
}

void Finder::StoreMC(const Ref::ChannelD& sexa) {
    // `DF::MC_Sexaquark` //
    fOutput_MC_ChannelD.Before.Fill_LorentzVector(sexa.BeforePx(), sexa.BeforePy(), sexa.BeforePz(), static_cast<float>(sexa.BeforeE()));
    fOutput_MC_ChannelD.After.Fill_LorentzVector(static_cast<float>(sexa.AsChannelD_AfterPx()), static_cast<float>(sexa.AsChannelD_AfterPy()),
                                                 static_cast<float>(sexa.AsChannelD_AfterPz()), static_cast<float>(sexa.AsChannelD_AfterE()));
    fOutput_MC_ChannelD.Nucleon.Fill_LorentzVector(sexa.NucleonPx(), sexa.NucleonPy(), sexa.NucleonPz(), static_cast<float>(sexa.NucleonE()));
    // -- secondary vertex
    fOutput_MC_ChannelD.SV.Fill_Coordinates(sexa.SV_X(), sexa.SV_Y(), sexa.SV_Z());
    // -- event properties
    fOutput_MC_ChannelD.PV = fInput_Event.MC_PV;
    // -- reaction id + flags
    fOutput_MC_ChannelD.ReactionID = sexa.AsChannelD_ReactionID();
    fOutput_MC_ChannelD.IsSignal = sexa.AsChannelD_IsSignal();
    fOutput_MC_ChannelD.IsHybrid = sexa.AsChannelD_IsHybrid();

    // `MC_ChannelD` //
    // V0 (`Flat::MCInfo_V0`)
    // -- `Flat::MCInfo_LV_Mother`
    //    -- `Flat::MCInfo`
    fOutput_MC_ChannelD.V0.Entry = sexa.v0.Entry();
    fOutput_MC_ChannelD.V0.PdgCode = sexa.v0.PdgCode();
    fOutput_MC_ChannelD.V0.ReactionID = sexa.v0.ReactionID();
    fOutput_MC_ChannelD.V0.IsTrue = sexa.v0.IsTrue();
    fOutput_MC_ChannelD.V0.IsSignal = sexa.v0.IsSignal();
    fOutput_MC_ChannelD.V0.IsSecondary = sexa.v0.IsSecondary();
    //    -- `Flat::LorentzVector`
    fOutput_MC_ChannelD.V0.Fill_LorentzVector(sexa.v0.Px(), sexa.v0.Py(), sexa.v0.Pz(), sexa.v0.Energy());
    //    -- mother info
    fOutput_MC_ChannelD.V0.Mother_Entry = sexa.v0.Mother_Entry();
    fOutput_MC_ChannelD.V0.Mother_PdgCode = sexa.v0.Mother_PdgCode();
    // -- Neg (`Flat::MCInfo_PxPyPz`)
    //    -- `Flat::MCInfo`
    fOutput_MC_ChannelD.V0.Neg.Entry = sexa.v0.Neg_Entry();
    fOutput_MC_ChannelD.V0.Neg.PdgCode = sexa.v0.Neg_PdgCode();
    fOutput_MC_ChannelD.V0.Neg.ReactionID = sexa.v0.Neg_ReactionID();
    fOutput_MC_ChannelD.V0.Neg.IsTrue = sexa.v0.Neg_IsTrue();
    fOutput_MC_ChannelD.V0.Neg.IsSignal = sexa.v0.Neg_IsSignal();
    fOutput_MC_ChannelD.V0.Neg.IsSecondary = sexa.v0.Neg_IsSecondary();
    //    -- `Flat::PxPyPz`
    fOutput_MC_ChannelD.V0.Neg.Fill_PxPyPz(sexa.v0.Neg_Px(), sexa.v0.Neg_Py(), sexa.v0.Neg_Pz());
    // -- Pos (`Flat::MCInfo_PxPyPz`)
    //    -- `Flat::MCInfo`
    fOutput_MC_ChannelD.V0.Pos.Entry = sexa.v0.Pos_Entry();
    fOutput_MC_ChannelD.V0.Pos.PdgCode = sexa.v0.Pos_PdgCode();
    fOutput_MC_ChannelD.V0.Pos.ReactionID = sexa.v0.Pos_ReactionID();
    fOutput_MC_ChannelD.V0.Pos.IsTrue = sexa.v0.Pos_IsTrue();
    fOutput_MC_ChannelD.V0.Pos.IsSignal = sexa.v0.Pos_IsSignal();
    fOutput_MC_ChannelD.V0.Pos.IsSecondary = sexa.v0.Pos_IsSecondary();
    //    -- `Flat::PxPyPz`
    fOutput_MC_ChannelD.V0.Pos.Fill_PxPyPz(sexa.v0.Pos_Px(), sexa.v0.Pos_Py(), sexa.v0.Pos_Pz());
    // -- V0 @ decay (`Flat::Coordinates`)
    fOutput_MC_ChannelD.V0.AtDecay.Fill_Coordinates(sexa.v0.DecayX(), sexa.v0.DecayY(), sexa.v0.DecayZ());
    // -- hybrid flag
    fOutput_MC_ChannelD.V0.IsHybrid = sexa.v0.IsHybrid();

    // Kaon (`Flat::MCInfo_Bachelor`)
    // -- `Flat::MCInfo_LV_Mother`
    //    -- `Flat::MCInfo`
    fOutput_MC_ChannelD.Kaon.Entry = sexa.kaon.Entry();
    fOutput_MC_ChannelD.Kaon.PdgCode = sexa.kaon.PdgCode();
    fOutput_MC_ChannelD.Kaon.ReactionID = sexa.kaon.ReactionID();
    fOutput_MC_ChannelD.Kaon.IsTrue = sexa.kaon.IsTrue();
    fOutput_MC_ChannelD.Kaon.Mother_Entry = sexa.kaon.Mother_Entry();
    fOutput_MC_ChannelD.Kaon.Mother_PdgCode = sexa.kaon.Mother_PdgCode();
    //    -- `Flat::LorentzVector`
    fOutput_MC_ChannelD.Kaon.Fill_LorentzVector(sexa.kaon.Px(), sexa.kaon.Py(), sexa.kaon.Pz(), sexa.kaon.Energy());
    //    -- mother info
    fOutput_MC_ChannelD.Kaon.IsSignal = sexa.kaon.IsSignal();
    fOutput_MC_ChannelD.Kaon.IsSecondary = sexa.kaon.IsSecondary();
    // -- grandmother info
    fOutput_MC_ChannelD.Kaon.GrandMother_Entry = sexa.kaon.GrandMother_Entry();
    fOutput_MC_ChannelD.Kaon.GrandMother_PdgCode = sexa.kaon.GrandMother_PdgCode();
}

// ## END OF CYCLES ## //

void Finder::EndOfAnalysis() {

    if (IsMC()) {
        fOutputTree_Injected->Write();
        Logger::Info(__FUNCTION__, "TTree \"{}\" has been written onto TFile {}", fOutputTree_Injected->GetName(), fSettings.PathOutputFile);
    }

    fOutputTree->Write();
    Logger::Info(__FUNCTION__, "TTree \"{}\" has been written onto TFile {}", fOutputTree->GetName(), fSettings.PathOutputFile);

    fCutFlowHist->Write();
    Logger::Info(__FUNCTION__, "TH1D \"{}\" has been written into TFile {}", fCutFlowHist->GetName(), fSettings.PathOutputFile);

    fCutFlowHist_Anti->Write();
    Logger::Info(__FUNCTION__, "TH1D \"{}\" has been written into TFile {}", fCutFlowHist_Anti->GetName(), fSettings.PathOutputFile);

    fInputChain_PackedEvents->ResetBranchAddresses();
    fOutputTree->ResetBranchAddresses();

    Logger::Info(__FUNCTION__, "All done.");
}

}  // namespace Tree2Secondaries
