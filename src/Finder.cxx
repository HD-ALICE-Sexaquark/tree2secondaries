#include <filesystem>
#include <memory>

#include "App/Logger.hxx"
#include "DataFormats/PackedEvents.hxx"
#include "Finder/Cuts.hxx"
#include "Finder/Finder.hxx"
#include "KF/Utilities.hxx"
#include "MC/Particles.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries {

bool Finder::Initialize() {

    fInputChain_PackedEvents = std::make_unique<TChain>(Name::PackedEvents.c_str());
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
    df.ReadBranches_PackedV0s(fInputChain_PackedEvents.get(), Particle::Acronym[pid]);
}

void Finder::ReadBranches_Tracks(EParticle pid, DF::Packed::Tracks& df) {
    df.ReadBranches_PackedTracks(fInputChain_PackedEvents.get(), Particle::Acronym[pid]);
}

void Finder::ReadBranches_LinkedV0s(EParticle pid, DF::Packed::LinkedV0s& df) {
    df.ReadBranches_LinkedV0s(fInputChain_PackedEvents.get(), Particle::Acronym[pid]);
}

void Finder::ReadBranches_LinkedTracks(EParticle pid, DF::Packed::LinkedTracks& df) {
    df.ReadBranches_LinkedTracks(fInputChain_PackedEvents.get(), Particle::Acronym[pid]);
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

    std::string tree_name{std::format("FoundSexa_Channel{}", static_cast<char>(fSettings.ReactionChannel))};

    fOutputTree = std::make_unique<TTree>(tree_name.c_str(), "");
    if (!fOutputTree) {
        Logger::Error(__FUNCTION__, "Couldn't create TTree \"{}\"", tree_name);
        return false;
    }

    return true;
}

void Finder::CreateOutputBranches(DF::ChannelA& df) { df.CreateBranches_ChannelA(fOutputTree.get(), IsMC()); }

void Finder::CreateOutputBranches(DF::ChannelD& df) {
    // PENDING
}

void Finder::CreateOutputBranches(DF::MC_ChannelA& df) { df.CreateBranches_MC_ChannelA(fOutputTree.get()); }

void Finder::CreateOutputBranches(DF::MC_ChannelD& df) {
    // PENDING
}

void Finder::CreateCutFlowHistogram() {
    const int x_nbins{20};
    const float x_min{0.};
    const float x_max{20.};
    std::string hist_title{";Cut N;N Passed Cut"};
    fCutFlowHist = std::make_unique<TH1D>("CutFlow", hist_title.c_str(), x_nbins, x_min, x_max);
}

// ## OUTPUT / Injected ZONE ## //

bool Finder::Injected_PrepareOutputTree() {

    std::string tree_name{Name::Injected};

    fOutputTree_Injected = std::make_unique<TTree>(tree_name.c_str(), "");
    if (!fOutputTree_Injected) {
        Logger::Error(__FUNCTION__, "Couldn't create TTree \"{}\"", tree_name);
        return false;
    }

    return true;
}

void Finder::Injected_CreateOutputBranches() { fOutput_Injected.CreateBranches_Flat_Injected(fOutputTree_Injected.get()); };

void Finder::Injected_FlattenAndStore() {

    double nucleon_mass{Particle::Mass[ReactionNucleonPID[fSettings.ReactionChannel]]};

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
        float nucleon_px{fInput_Injected.Nucleon_Px->at(idx_inj)};
        float nucleon_py{fInput_Injected.Nucleon_Py->at(idx_inj)};
        float nucleon_pz{fInput_Injected.Nucleon_Pz->at(idx_inj)};
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

    const DF::Packed::V0s* Packed_Lambdas{anti_channel ? &fInput_Lambdas : &fInput_AntiLambdas};
    const DF::Packed::LinkedV0s* MC_Lambdas{anti_channel ? &fInput_Linked_Lambdas : &fInput_Linked_AntiLambdas};
    EParticle PID_Lambda{anti_channel ? EParticle::Lambda : EParticle::AntiLambda};
    EParticle PID_Lambda_Neg{anti_channel ? EParticle::PiMinus : EParticle::AntiProton};
    EParticle PID_Lambda_Pos{anti_channel ? EParticle::Proton : EParticle::PiPlus};

    // loop over all possible pairs of (anti)lambda + K0S //
    auto n_lambdas = static_cast<int>(Packed_Lambdas->Entry->size());
    auto n_k0s = static_cast<int>(fInput_KaonsZeroShort.Entry->size());
    for (int idx_v0a{0}; idx_v0a < n_lambdas; ++idx_v0a) {

        // unpack (anti)lambda //
        auto v0a = KF::UnpackV0(*Packed_Lambdas, idx_v0a, PID_Lambda, PID_Lambda_Neg, PID_Lambda_Pos);

        for (int idx_v0b{0}; idx_v0b < n_k0s; ++idx_v0b) {

            // unpack K0S //
            auto v0b = KF::UnpackV0(fInput_KaonsZeroShort, idx_v0b, EParticle::KaonZeroShort, EParticle::PiMinus, EParticle::PiPlus);

            // sanity check //
            if (v0a.Neg.idx == v0b.Neg.idx || v0a.Neg.idx == v0b.Pos.idx || v0a.Pos.idx == v0b.Neg.idx || v0a.Pos.idx == v0b.Pos.idx) continue;

            // fit //
            KF::ChannelA sexa{v0a, v0b};
            sexa.DoFit(fInput_Event.MagneticField);

#ifdef T2S_DEBUG
            Logger::Debug(__FUNCTION__, "idx(v0a,neg,pos)={},{},{}", v0a.idx, v0a.Neg.idx, v0a.Pos.idx);
            Logger::Debug(__FUNCTION__, ";x,y,z={},{},{}", v0a.X(), v0a.Y(), v0a.Z());
            Logger::Debug(__FUNCTION__, ";px,py,pz={},{},{}", v0a.Px(), v0a.Py(), v0a.Pz());
            Logger::Debug(__FUNCTION__, ";mass={}", v0a.Mass());

            Logger::Debug(__FUNCTION__, "idx(v0b,neg,pos)={},{},{}", v0b.idx, v0b.Neg.idx, v0b.Pos.idx);
            Logger::Debug(__FUNCTION__, ";x,y,z={},{},{}", v0b.X(), v0b.Y(), v0b.Z());
            Logger::Debug(__FUNCTION__, ";px,py,pz={},{},{}", v0b.Px(), v0b.Py(), v0b.Pz());
            Logger::Debug(__FUNCTION__, ";mass={}", v0b.Mass());

            Logger::Debug(__FUNCTION__, "x,y,z={},{},{}", sexa.X(), sexa.Y(), sexa.Z());
            Logger::Debug(__FUNCTION__, ";x,y,z(v0a)={},{},{}", sexa.V0A_PCA_XYZ()[0], sexa.V0A_PCA_XYZ()[1], sexa.V0A_PCA_XYZ()[2]);
            Logger::Debug(__FUNCTION__, ";x,y,z(v0b)={},{},{}", sexa.V0B_PCA_XYZ()[0], sexa.V0B_PCA_XYZ()[1], sexa.V0B_PCA_XYZ()[2]);
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

            MC::V0 mc_v0a_debug{*MC_Lambdas, v0a.idx};
            MC::V0 mc_v0b_debug{fInput_Linked_KaonsZeroShort, v0b.idx};
            MC::ChannelA mc_sexa_debug{fInput_Injected, fSettings.SexaquarkMass, mc_v0a_debug, mc_v0b_debug};
            if (IsMC()) {
                Logger::Debug(__FUNCTION__, "sexa(is_signal,reaction_id,is_hybrid)={},{},{}", mc_sexa_debug.IsSignal, mc_sexa_debug.ReactionID,
                              mc_sexa_debug.IsHybrid);
                Logger::Debug(__FUNCTION__, ";v0a(entry,pdg_code)={},{}", mc_sexa_debug.V0A.Entry, mc_sexa_debug.V0A.PdgCode);
                Logger::Debug(__FUNCTION__, ";v0a(is_signal,reaction_id,is_hybrid)={},{},{}", mc_sexa_debug.V0A.IsSignal,
                              mc_sexa_debug.V0A.ReactionID, mc_sexa_debug.V0A.IsHybrid);
                Logger::Debug(__FUNCTION__, ";v0b(entry,pdg_code)={},{}", mc_sexa_debug.V0B.Entry, mc_sexa_debug.V0B.PdgCode);
                Logger::Debug(__FUNCTION__, ";v0b(is_signal,reaction_id,is_hybrid)={},{},{}", mc_sexa_debug.V0B.IsSignal,
                              mc_sexa_debug.V0B.ReactionID, mc_sexa_debug.V0B.IsHybrid);
            }
#endif

            // apply cuts //
            if (!PassesCuts(sexa)) continue;

            // store //
            Store(sexa, anti_channel);
            if (IsMC()) {
                MC::V0 mc_v0a{*MC_Lambdas, v0a.idx};
                MC::V0 mc_v0b{fInput_Linked_KaonsZeroShort, v0b.idx};
                MC::ChannelA mc_sexa{fInput_Injected, fSettings.SexaquarkMass, mc_v0a, mc_v0b};
                StoreMC(mc_sexa);
            }
            fOutputTree->Fill();
        }
    }
}

bool Finder::PassesCuts(const KF::ChannelA& sexa) const {

    fCutFlowHist->Fill(0.);
    if (sexa.Radius2D() < Cuts::ChannelA::Min_Radius2D || sexa.Radius2D() > Cuts::ChannelA::Max_Radius2D) return false;
    fCutFlowHist->Fill(1.);
    if (sexa.DecayLength_V0A() > Cuts::ChannelA::Max_DecayLengthLa) return false;
    fCutFlowHist->Fill(2.);
    if (sexa.DecayLength_V0B() > Cuts::ChannelA::Max_DecayLengthK0) return false;
    fCutFlowHist->Fill(3.);
    if (sexa.AbsRapidity_MinusNucleon() > Cuts::ChannelA::AbsMax_Rapidity) return false;  // PENDING: kinematic cut, affected by Fermi motion
    fCutFlowHist->Fill(4.);
    if (sexa.Mass_MinusNucleon() < Cuts::ChannelA::Min_MassMinusNucleon || sexa.Mass_MinusNucleon() > Cuts::ChannelA::Max_MassMinusNucleon) {
        return false;  // PENDING: kinematic cut, affected by Fermi motion
    }
    fCutFlowHist->Fill(5.);
    if (sexa.CPA_Point(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) < Cuts::ChannelA::Min_CPAwrtPV ||
        sexa.CPA_Point(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) > Cuts::ChannelA::Max_CPAwrtPV) {
        return false;  // PENDING: kinematic cut, affected by Fermi motion
    }
    fCutFlowHist->Fill(6.);
    if (sexa.DCA_V0ANeg_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelA::Max_DCALaNegSV) return false;
    fCutFlowHist->Fill(7.);
    if (sexa.DCA_V0APos_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelA::Max_DCALaPosSV) return false;
    fCutFlowHist->Fill(8.);
    if (sexa.DCA_V0BNeg_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelA::Max_DCAK0NegSV) return false;
    fCutFlowHist->Fill(9.);
    if (sexa.DCA_V0BPos_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelA::Max_DCAK0PosSV) return false;
    fCutFlowHist->Fill(10.);
    if (sexa.DCA_V0A_wrt_SV() > Cuts::ChannelA::Max_DCALaSV) return false;
    fCutFlowHist->Fill(11.);
    if (sexa.DCA_V0B_wrt_SV() > Cuts::ChannelA::Max_DCAK0SV) return false;
    fCutFlowHist->Fill(12.);
    if (sexa.DCA_btw_V0s() > Cuts::ChannelA::Max_DCAbtwV0s) return false;
    fCutFlowHist->Fill(13.);

    return true;
}

void Finder::Store(const KF::ChannelA& sexa, bool anti_channel) {
    // `DF::Sexaquark`
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
    fOutput_ChannelA.Chi2NDF = Const::DummyFloat;  // PENDING!
    // -- extra info
    fOutput_ChannelA.E_MinusNucleon = static_cast<float>(sexa.E_MinusNucleon());
    fOutput_ChannelA.AntiChannel = anti_channel;
    // `DF::ChannelA`
    // -- V0A
    fOutput_ChannelA.V0A_atPCA.Fill_State(Const::DummyFloat, Const::DummyFloat, Const::DummyFloat,                        //
                                          Const::DummyFloat, Const::DummyFloat, Const::DummyFloat, Const::DummyFloat);    // PENDING!
    fOutput_ChannelA.V0A_atDecay.Fill_State(Const::DummyFloat, Const::DummyFloat, Const::DummyFloat,                      //
                                            Const::DummyFloat, Const::DummyFloat, Const::DummyFloat, Const::DummyFloat);  // PENDING!
    // -- V0B
    fOutput_ChannelA.V0B_atPCA.Fill_State(Const::DummyFloat, Const::DummyFloat, Const::DummyFloat,                        //
                                          Const::DummyFloat, Const::DummyFloat, Const::DummyFloat, Const::DummyFloat);    // PENDING!
    fOutput_ChannelA.V0B_atDecay.Fill_State(Const::DummyFloat, Const::DummyFloat, Const::DummyFloat,                      //
                                            Const::DummyFloat, Const::DummyFloat, Const::DummyFloat, Const::DummyFloat);  // PENDING!
    // -- V0A daughters
    fOutput_ChannelA.V0A_Neg_atPCA.Fill_State(Const::DummyFloat, Const::DummyFloat, Const::DummyFloat,                      //
                                              Const::DummyFloat, Const::DummyFloat, Const::DummyFloat, Const::DummyFloat);  // PENDING!
    fOutput_ChannelA.V0A_Neg_atV0.Fill_State(Const::DummyFloat, Const::DummyFloat, Const::DummyFloat,                       //
                                             Const::DummyFloat, Const::DummyFloat, Const::DummyFloat, Const::DummyFloat);   // PENDING!
    fOutput_ChannelA.V0A_Pos_atPCA.Fill_State(Const::DummyFloat, Const::DummyFloat, Const::DummyFloat,                      //
                                              Const::DummyFloat, Const::DummyFloat, Const::DummyFloat, Const::DummyFloat);  // PENDING!
    fOutput_ChannelA.V0A_Pos_atV0.Fill_State(Const::DummyFloat, Const::DummyFloat, Const::DummyFloat,                       //
                                             Const::DummyFloat, Const::DummyFloat, Const::DummyFloat, Const::DummyFloat);   // PENDING!
    // -- V0B daughters
    fOutput_ChannelA.V0B_Neg_atPCA.Fill_State(Const::DummyFloat, Const::DummyFloat, Const::DummyFloat,                      //
                                              Const::DummyFloat, Const::DummyFloat, Const::DummyFloat, Const::DummyFloat);  // PENDING!
    fOutput_ChannelA.V0B_Neg_atV0.Fill_State(Const::DummyFloat, Const::DummyFloat, Const::DummyFloat,                       //
                                             Const::DummyFloat, Const::DummyFloat, Const::DummyFloat, Const::DummyFloat);   // PENDING!
    fOutput_ChannelA.V0B_Pos_atPCA.Fill_State(Const::DummyFloat, Const::DummyFloat, Const::DummyFloat,                      //
                                              Const::DummyFloat, Const::DummyFloat, Const::DummyFloat, Const::DummyFloat);  // PENDING!
    fOutput_ChannelA.V0B_Pos_atV0.Fill_State(Const::DummyFloat, Const::DummyFloat, Const::DummyFloat,                       //
                                             Const::DummyFloat, Const::DummyFloat, Const::DummyFloat, Const::DummyFloat);   // PENDING!
    // -- indices
    fOutput_ChannelA.V0A_Entry = sexa.V0A.idx;
    fOutput_ChannelA.V0B_Entry = sexa.V0B.idx;
    fOutput_ChannelA.V0A_Neg_Entry = sexa.V0A.Neg.idx;
    fOutput_ChannelA.V0A_Pos_Entry = sexa.V0A.Pos.idx;
    fOutput_ChannelA.V0B_Neg_Entry = sexa.V0B.Neg.idx;
    fOutput_ChannelA.V0B_Pos_Entry = sexa.V0B.Pos.idx;
}

void Finder::StoreMC(const MC::ChannelA& sexa) {
    // `Found::MC_Sexaquark`
    fOutput_MC_ChannelA.Before.Fill_LorentzVector(static_cast<float>(sexa.BeforePx), static_cast<float>(sexa.BeforePy),
                                                  static_cast<float>(sexa.BeforePz), static_cast<float>(sexa.BeforeE));
    fOutput_MC_ChannelA.After.Fill_LorentzVector(static_cast<float>(sexa.AfterPx), static_cast<float>(sexa.AfterPy),  //
                                                 static_cast<float>(sexa.AfterPz), static_cast<float>(sexa.AfterE));
    fOutput_MC_ChannelA.Nucleon.Fill_LorentzVector(static_cast<float>(sexa.NucleonPx), static_cast<float>(sexa.NucleonPy),
                                                   static_cast<float>(sexa.NucleonPz), static_cast<float>(sexa.NucleonE));
    // -- secondary vertex
    fOutput_MC_ChannelA.SV.X = static_cast<float>(sexa.X);
    fOutput_MC_ChannelA.SV.Y = static_cast<float>(sexa.Y);
    fOutput_MC_ChannelA.SV.Z = static_cast<float>(sexa.Z);
    // -- event properties
    fOutput_MC_ChannelA.PV = fInput_Event.MC_PV;
    // -- reaction id + flags
    fOutput_MC_ChannelA.ReactionID = sexa.ReactionID;
    fOutput_MC_ChannelA.IsSignal = sexa.IsSignal;
    fOutput_MC_ChannelA.IsHybrid = sexa.IsHybrid;
    // `Found::MC_ChannelA`
    // -- V0A
    fOutput_MC_ChannelA.V0A.Fill_LorentzVector(sexa.V0A.Px, sexa.V0A.Py, sexa.V0A.Pz, sexa.V0A.Energy);
    fOutput_MC_ChannelA.V0A.Entry = sexa.V0A.Entry;
    fOutput_MC_ChannelA.V0A.PdgCode = sexa.V0A.PdgCode;
    fOutput_MC_ChannelA.V0A.Mother_Entry = sexa.V0A.Mother_Entry;
    fOutput_MC_ChannelA.V0A.Mother_PdgCode = sexa.V0A.Mother_PdgCode;
    fOutput_MC_ChannelA.V0A.ReactionID = sexa.V0A.ReactionID;
    fOutput_MC_ChannelA.V0A.IsTrue = sexa.V0A.IsTrue;
    fOutput_MC_ChannelA.V0A.IsSignal = sexa.V0A.IsSignal;
    fOutput_MC_ChannelA.V0A.IsSecondary = sexa.V0A.IsSecondary;
    // -- V0A neg
    fOutput_MC_ChannelA.V0A_Neg.Fill_LorentzVector(Const::DummyInt, Const::DummyInt, Const::DummyInt, Const::DummyInt);  // PENDING!
    fOutput_MC_ChannelA.V0A_Neg.Entry = sexa.V0A.neg.Entry;
    fOutput_MC_ChannelA.V0A_Neg.PdgCode = sexa.V0A.neg.PdgCode;
    fOutput_MC_ChannelA.V0A_Neg.Mother_Entry = Const::DummyInt;    // PENDING!
    fOutput_MC_ChannelA.V0A_Neg.Mother_PdgCode = Const::DummyInt;  // PENDING!
    fOutput_MC_ChannelA.V0A_Neg.ReactionID = Const::DummyInt;      // PENDING!
    fOutput_MC_ChannelA.V0A_Neg.IsTrue = Const::DummyInt;          // PENDING!
    fOutput_MC_ChannelA.V0A_Neg.IsSignal = Const::DummyInt;        // PENDING!
    fOutput_MC_ChannelA.V0A_Neg.IsSecondary = Const::DummyInt;     // PENDING!
    // -- V0A pos
    fOutput_MC_ChannelA.V0A_Pos.Fill_LorentzVector(Const::DummyInt, Const::DummyInt, Const::DummyInt, Const::DummyInt);  // PENDING!
    fOutput_MC_ChannelA.V0A_Pos.Entry = sexa.V0A.pos.Entry;
    fOutput_MC_ChannelA.V0A_Pos.PdgCode = sexa.V0A.pos.PdgCode;
    fOutput_MC_ChannelA.V0A_Pos.Mother_Entry = Const::DummyInt;    // PENDING!
    fOutput_MC_ChannelA.V0A_Pos.Mother_PdgCode = Const::DummyInt;  // PENDING!
    fOutput_MC_ChannelA.V0A_Pos.ReactionID = Const::DummyInt;      // PENDING!
    fOutput_MC_ChannelA.V0A_Pos.IsTrue = Const::DummyInt;          // PENDING!
    fOutput_MC_ChannelA.V0A_Pos.IsSignal = Const::DummyInt;        // PENDING!
    fOutput_MC_ChannelA.V0A_Pos.IsSecondary = Const::DummyInt;     // PENDING!
    // -- V0B
    fOutput_MC_ChannelA.V0B.Fill_LorentzVector(sexa.V0B.Px, sexa.V0B.Py, sexa.V0B.Pz, sexa.V0B.Energy);
    fOutput_MC_ChannelA.V0B.Entry = sexa.V0B.Entry;
    fOutput_MC_ChannelA.V0B.PdgCode = sexa.V0B.PdgCode;
    fOutput_MC_ChannelA.V0B.Mother_Entry = sexa.V0B.Mother_Entry;
    fOutput_MC_ChannelA.V0B.Mother_PdgCode = sexa.V0B.Mother_PdgCode;
    fOutput_MC_ChannelA.V0B.ReactionID = sexa.V0B.ReactionID;
    fOutput_MC_ChannelA.V0B.IsTrue = sexa.V0B.IsTrue;
    fOutput_MC_ChannelA.V0B.IsSignal = sexa.V0B.IsSignal;
    fOutput_MC_ChannelA.V0B.IsSecondary = sexa.V0B.IsSecondary;
    // -- V0A neg
    fOutput_MC_ChannelA.V0B_Neg.Fill_LorentzVector(Const::DummyInt, Const::DummyInt, Const::DummyInt, Const::DummyInt);  // PENDING!
    fOutput_MC_ChannelA.V0B_Neg.Entry = sexa.V0B.neg.Entry;
    fOutput_MC_ChannelA.V0B_Neg.PdgCode = sexa.V0B.neg.PdgCode;
    fOutput_MC_ChannelA.V0B_Neg.Mother_Entry = Const::DummyInt;    // PENDING!
    fOutput_MC_ChannelA.V0B_Neg.Mother_PdgCode = Const::DummyInt;  // PENDING!
    fOutput_MC_ChannelA.V0B_Neg.ReactionID = Const::DummyInt;      // PENDING!
    fOutput_MC_ChannelA.V0B_Neg.IsTrue = Const::DummyInt;          // PENDING!
    fOutput_MC_ChannelA.V0B_Neg.IsSignal = Const::DummyInt;        // PENDING!
    fOutput_MC_ChannelA.V0B_Neg.IsSecondary = Const::DummyInt;     // PENDING!
    // -- V0B pos
    fOutput_MC_ChannelA.V0B_Pos.Fill_LorentzVector(Const::DummyInt, Const::DummyInt, Const::DummyInt, Const::DummyInt);  // PENDING!
    fOutput_MC_ChannelA.V0B_Pos.Entry = sexa.V0B.pos.Entry;
    fOutput_MC_ChannelA.V0B_Pos.PdgCode = sexa.V0B.pos.PdgCode;
    fOutput_MC_ChannelA.V0B_Pos.Mother_Entry = Const::DummyInt;    // PENDING!
    fOutput_MC_ChannelA.V0B_Pos.Mother_PdgCode = Const::DummyInt;  // PENDING!
    fOutput_MC_ChannelA.V0B_Pos.ReactionID = Const::DummyInt;      // PENDING!
    fOutput_MC_ChannelA.V0B_Pos.IsTrue = Const::DummyInt;          // PENDING!
    fOutput_MC_ChannelA.V0B_Pos.IsSignal = Const::DummyInt;        // PENDING!
    fOutput_MC_ChannelA.V0B_Pos.IsSecondary = Const::DummyInt;     // PENDING!
    // -- V0A @ decay
    fOutput_MC_ChannelA.V0A_atDecay.Fill_Coordinates(Const::DummyInt, Const::DummyInt, Const::DummyInt);
    // -- V0B @ decay
    fOutput_MC_ChannelA.V0B_atDecay.Fill_Coordinates(Const::DummyInt, Const::DummyInt, Const::DummyInt);
    // hybrid flags
    fOutput_MC_ChannelA.V0A_IsHybrid = sexa.V0A.IsHybrid;
    fOutput_MC_ChannelA.V0B_IsHybrid = sexa.V0B.IsHybrid;
}

// ## Channel D ZONE ## //

void Finder::FindSexaquarks_ChannelD(bool anti_channel) {

    const DF::Packed::V0s* Packed_Lambdas{anti_channel ? &fInput_Lambdas : &fInput_AntiLambdas};
    const DF::Packed::Tracks* Packed_Kaons{anti_channel ? &fInput_NegKaons : &fInput_PosKaons};
    const DF::Packed::LinkedV0s* Linked_Lambdas{anti_channel ? &fInput_Linked_Lambdas : &fInput_Linked_AntiLambdas};
    const DF::Packed::LinkedTracks* Linked_Kaons{anti_channel ? &fInput_Linked_NegKaons : &fInput_Linked_PosKaons};
    EParticle PID_Lambda{anti_channel ? EParticle::Lambda : EParticle::AntiLambda};
    EParticle PID_Lambda_Neg{anti_channel ? EParticle::PiMinus : EParticle::AntiProton};
    EParticle PID_Lambda_Pos{anti_channel ? EParticle::Proton : EParticle::PiPlus};
    EParticle PID_Kaon{anti_channel ? EParticle::NegKaon : EParticle::PosKaon};
    double Mass_Kaon{Particle::Mass[PID_Kaon]};
    int Charge_Kaon{Particle::Charge[PID_Kaon]};

    // loop over all possible pairs of (anti)lambda + (pos/neg)kaon //
    auto n_lambdas = static_cast<int>(Packed_Lambdas->Entry->size());
    auto n_kaons = static_cast<int>(Packed_Kaons->Entry->size());
    for (int idx_v0{0}; idx_v0 < n_lambdas; ++idx_v0) {

        // unpack (anti)lambda //
        auto v0 = KF::UnpackV0(*Packed_Lambdas, idx_v0, PID_Lambda, PID_Lambda_Neg, PID_Lambda_Pos);

        for (int idx_kaon{0}; idx_kaon < n_kaons; ++idx_kaon) {

            // unpack kaon //
            auto kaon = KF::UnpackTrack(*Packed_Kaons, Charge_Kaon, Mass_Kaon, idx_kaon);

            // sanity check //
            if (v0.Neg.idx == kaon.idx || v0.Pos.idx == kaon.idx) continue;

            // fit //
            KF::ChannelD sexa{v0, kaon};
            sexa.DoFit(fInput_Event.MagneticField);

            // apply cuts //
            if (!PassesCuts(sexa)) continue;
#ifdef T2S_DEBUG
            Logger::Debug(__FUNCTION__, "idx(sexa,v0,neg,pos,kaon)={},{},{},{}", sexa.V0.idx, sexa.V0.Neg.idx, sexa.V0.Pos.idx, sexa.Kaon.idx);
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
                MC::V0 mc_v0{*Linked_Lambdas, v0.idx};
                MC::Track mc_kaon{*Linked_Kaons, kaon.idx};
                MC::ChannelD mc_sexa{fInput_Injected, fSettings.SexaquarkMass, mc_v0, mc_kaon};
                StoreMC(mc_sexa);
            }
            fOutputTree->Fill();
        }
    }
}

bool Finder::PassesCuts(const KF::ChannelD& sexa) const {

    fCutFlowHist->Fill(0.);
    if (sexa.Radius2D() < Cuts::ChannelD::Min_Radius2D || sexa.Radius2D() > Cuts::ChannelD::Max_Radius2D) return false;
    fCutFlowHist->Fill(1.);
    if (sexa.AbsRapidity_MinusNucleon() > Cuts::ChannelD::AbsMax_Rapidity) return false;  // PENDING: kinematics, affected by Fermi motion
    fCutFlowHist->Fill(2.);
    if (sexa.CPA_Point(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) < Cuts::ChannelD::Min_CPAwrtPV ||
        sexa.CPA_Point(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) > Cuts::ChannelD::Max_CPAwrtPV) {
        return false;  // PENDING: kinematics, affected by Fermi motion
    }
    fCutFlowHist->Fill(3.);
    if (sexa.DCA_V0_wrt_SV() > Cuts::ChannelD::Max_DCALaSV) return false;
    fCutFlowHist->Fill(4.);
    if (sexa.DCA_Kaon_wrt_SV() > Cuts::ChannelD::Max_DCAKaSV) return false;
    fCutFlowHist->Fill(5.);
    if (sexa.DCA_V0Neg_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelD::Max_DCALaNegSV) return false;
    fCutFlowHist->Fill(6.);
    if (sexa.DCA_V0Pos_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelD::Max_DCALaPosSV) return false;
    fCutFlowHist->Fill(7.);
    if (sexa.DCA_btw_V0_Kaon() > Cuts::ChannelD::Max_DCAKaLa) return false;
    fCutFlowHist->Fill(8.);
    if (sexa.DCA_btw_V0Neg_Kaon(fInput_Event.MagneticField) > Cuts::ChannelD::Max_DCALaNegKa) return false;
    fCutFlowHist->Fill(9.);
    if (sexa.DCA_btw_V0Pos_Kaon(fInput_Event.MagneticField) > Cuts::ChannelD::Max_DCALaPosKa) return false;
    fCutFlowHist->Fill(10.);

    return true;
}

void Finder::Store(const KF::ChannelD& sexa, bool anti_channel) {}

void Finder::StoreMC(const MC::ChannelD& sexa) {}

// ## END OF CYCLES ## //

void Finder::EndOfAnalysis() {

    if (IsMC()) {
        fOutputTree_Injected->Write();
        Logger::Info(__FUNCTION__, "TTree \"{}\" has been written onto TFile {}", fOutputTree_Injected->GetName(), fSettings.PathOutputFile);
    }

    fOutputTree->Write();
    Logger::Info(__FUNCTION__, "TTree \"{}\" has been written onto TFile {}", fOutputTree->GetName(), fSettings.PathOutputFile);

    fCutFlowHist->Write();

    fInputChain_PackedEvents->ResetBranchAddresses();
    fOutputTree->ResetBranchAddresses();

    Logger::Info(__FUNCTION__, "All done.");
}

}  // namespace Tree2Secondaries
