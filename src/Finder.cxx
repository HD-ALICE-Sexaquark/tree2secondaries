#include <filesystem>
#include <memory>

#include "App/Logger.hxx"
#include "Finder/Finder.hxx"
#include "Finder/FinderCuts.hxx"
#include "Math/Constants.hxx"
#include "Storage/Vector/VectorTracks.hxx"
#include "Storage/Vector/VectorV0s.hxx"
#include "Truth/TruthSexaquark.hxx"
#include "View/MC/ViewMcInjected.hxx"
#include "View/MC/ViewMcPackedTrack.hxx"
#include "View/MC/ViewMcPackedV0.hxx"
#ifdef T2S_LEGACY_KF
#include "Fit/Legacy/FitChannelA_Legacy.hxx"
#include "Fit/Legacy/FitChannelD_Legacy.hxx"
#include "Fit/Legacy/FitTrack_Legacy.hxx"
#include "Math/Legacy/BaseMath_Legacy.hxx"
#include "View/Reconstructed/Legacy/ViewTrack_Legacy.hxx"
#else
#include "Fit/FitChannelA.hxx"
#include "Fit/FitChannelD.hxx"
#include "Fit/FitTrack.hxx"
#include "Math/BaseMath.hxx"
#include "View/Reconstructed/ViewTrack.hxx"
#endif

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

    PrepareOutputHistograms();

    if (!PrepareOutputTree()) return false;
    CreateOutputBranches();

    if (IsMC()) {
        if (!Injected_PrepareOutputTree()) return false;
        fOutput_Injected.CreateBranches_FlatInjected(fOutputTree_Injected.get());
    }

    Logger::Info(__FUNCTION__, "Finder initialized successfully.");

    return true;
}

// ## INPUT ZONE ## //

void Finder::ReadInputBranches() {
    // by default, turn off all branches
    fInputChain_PackedEvents->SetBranchStatus("*", false);
    // connect input branches to memory
    fInput_Event.ReadBranches_FlatEvent(fInputChain_PackedEvents.get(), IsMC());
    if (IsMC()) {
        fInput_MC_PV.ReadBranches_FlatCoordinates(fInputChain_PackedEvents.get(), "MC_PV", "v");
        fInput_Injected.ReadBranches_VectorInjected(fInputChain_PackedEvents.get(), true);
    }
    // -- depending on reaction channels
    switch (GetReactionChannel()) {
        case EReactionChannel::A:
            fInput_AntiLambdas.ReadBranches_VectorV0s(fInputChain_PackedEvents.get(), Const::Particle_Acronym[EParticle::AntiLambda]);
            fInput_Lambdas.ReadBranches_VectorV0s(fInputChain_PackedEvents.get(), Const::Particle_Acronym[EParticle::Lambda]);
            fInput_KaonsZeroShort.ReadBranches_VectorV0s(fInputChain_PackedEvents.get(), Const::Particle_Acronym[EParticle::KaonZeroShort]);
            if (IsMC()) {
                fInput_MC_AntiLambdas.ReadBranches_VectorMC_V0s(fInputChain_PackedEvents.get(), Const::Particle_Acronym[EParticle::AntiLambda]);
                fInput_MC_Lambdas.ReadBranches_VectorMC_V0s(fInputChain_PackedEvents.get(), Const::Particle_Acronym[EParticle::Lambda]);
                fInput_MC_KaonsZeroShort.ReadBranches_VectorMC_V0s(fInputChain_PackedEvents.get(), Const::Particle_Acronym[EParticle::KaonZeroShort]);
            }
            break;
        case EReactionChannel::D:
            fInput_AntiLambdas.ReadBranches_VectorV0s(fInputChain_PackedEvents.get(), Const::Particle_Acronym[EParticle::AntiLambda]);
            fInput_Lambdas.ReadBranches_VectorV0s(fInputChain_PackedEvents.get(), Const::Particle_Acronym[EParticle::Lambda]);
            fInput_NegKaons.ReadBranches_VectorTracks(fInputChain_PackedEvents.get(), IsMC(), true, Const::Particle_Acronym[EParticle::NegKaon]);
            fInput_PosKaons.ReadBranches_VectorTracks(fInputChain_PackedEvents.get(), IsMC(), true, Const::Particle_Acronym[EParticle::PosKaon]);
            if (IsMC()) {
                fInput_MC_AntiLambdas.ReadBranches_VectorMC_V0s(fInputChain_PackedEvents.get(), Const::Particle_Acronym[EParticle::AntiLambda]);
                fInput_MC_Lambdas.ReadBranches_VectorMC_V0s(fInputChain_PackedEvents.get(), Const::Particle_Acronym[EParticle::Lambda]);
                fInput_MC_NegKaons.ReadBranches_VectorMC_Tracks(fInputChain_PackedEvents.get(), Const::Particle_Acronym[EParticle::NegKaon]);
                fInput_MC_PosKaons.ReadBranches_VectorMC_Tracks(fInputChain_PackedEvents.get(), Const::Particle_Acronym[EParticle::PosKaon]);
            }
            break;
        case EReactionChannel::E:
            fInput_AntiLambdas.ReadBranches_VectorV0s(fInputChain_PackedEvents.get(), Const::Particle_Acronym[EParticle::AntiLambda]);
            fInput_Lambdas.ReadBranches_VectorV0s(fInputChain_PackedEvents.get(), Const::Particle_Acronym[EParticle::Lambda]);
            fInput_NegKaons.ReadBranches_VectorTracks(fInputChain_PackedEvents.get(), IsMC(), true, Const::Particle_Acronym[EParticle::NegKaon]);
            fInput_PosKaons.ReadBranches_VectorTracks(fInputChain_PackedEvents.get(), IsMC(), true, Const::Particle_Acronym[EParticle::PosKaon]);
            fInput_PiMinus.ReadBranches_VectorTracks(fInputChain_PackedEvents.get(), IsMC(), true, Const::Particle_Acronym[EParticle::PiMinus]);
            fInput_PiPlus.ReadBranches_VectorTracks(fInputChain_PackedEvents.get(), IsMC(), true, Const::Particle_Acronym[EParticle::PiPlus]);
            if (IsMC()) {
                fInput_MC_AntiLambdas.ReadBranches_VectorMC_V0s(fInputChain_PackedEvents.get(), Const::Particle_Acronym[EParticle::AntiLambda]);
                fInput_MC_Lambdas.ReadBranches_VectorMC_V0s(fInputChain_PackedEvents.get(), Const::Particle_Acronym[EParticle::Lambda]);
                fInput_MC_NegKaons.ReadBranches_VectorMC_Tracks(fInputChain_PackedEvents.get(), Const::Particle_Acronym[EParticle::NegKaon]);
                fInput_MC_PosKaons.ReadBranches_VectorMC_Tracks(fInputChain_PackedEvents.get(), Const::Particle_Acronym[EParticle::PosKaon]);
                fInput_MC_PiMinus.ReadBranches_VectorMC_Tracks(fInputChain_PackedEvents.get(), Const::Particle_Acronym[EParticle::PiMinus]);
                fInput_MC_PiPlus.ReadBranches_VectorMC_Tracks(fInputChain_PackedEvents.get(), Const::Particle_Acronym[EParticle::PiPlus]);
            }
            break;
        case EReactionChannel::H:
            fInput_NegKaons.ReadBranches_VectorTracks(fInputChain_PackedEvents.get(), IsMC(), true, Const::Particle_Acronym[EParticle::NegKaon]);
            fInput_PosKaons.ReadBranches_VectorTracks(fInputChain_PackedEvents.get(), IsMC(), true, Const::Particle_Acronym[EParticle::PosKaon]);
            if (IsMC()) {
                fInput_MC_NegKaons.ReadBranches_VectorMC_Tracks(fInputChain_PackedEvents.get(), Const::Particle_Acronym[EParticle::NegKaon]);
                fInput_MC_PosKaons.ReadBranches_VectorMC_Tracks(fInputChain_PackedEvents.get(), Const::Particle_Acronym[EParticle::PosKaon]);
            }
            break;
    }  // end of switch statement
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

void Finder::PrepareOutputHistograms() {
    // event counter //
    fHist_EventCounter = std::make_unique<TH1D>("N_Events", ";;N_Events", 1, 0, 1);
    // cut flows //
    const int x_nbins{20};
    const float x_min{0.};
    const float x_max{20.};
    std::string hist_title{";Cut N;N Passed Cut"};
    fHist_CutFlow = std::make_unique<TH1D>("CutFlow", hist_title.c_str(), x_nbins, x_min, x_max);
    fHist_CutFlow_AntiChannel = std::make_unique<TH1D>("CutFlow_Anti", hist_title.c_str(), x_nbins, x_min, x_max);
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

void Finder::CreateOutputBranches() {
    switch (GetReactionChannel()) {
        case EReactionChannel::A:
            fOutput_ChannelA.CreateBranches_FlatChannelA(fOutputTree.get(), IsMC());
            if (IsMC()) fOutput_MC_ChannelA.CreateBranches_FlatMC_ChannelA(fOutputTree.get());
            break;
        case EReactionChannel::D:
            fOutput_ChannelD.CreateBranches_FlatChannelD(fOutputTree.get(), IsMC());
            if (IsMC()) fOutput_MC_ChannelD.CreateBranches_FlatMC_ChannelD(fOutputTree.get());
            break;
        default:
            break;
    }
}

// ## Event ZONE ## //

void Finder::ProcessEvent() {  //
    fHist_EventCounter->Fill(0.);
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

void Finder::ProcessInjected() {

    auto n_injected = static_cast<int>(fInput_Injected.ReactionID->size());
    for (int entry{0}; entry < n_injected; ++entry) {

        View::MC::Injected inj_view{&fInput_Injected, entry};

        // `Flat::Injected`
        // -- `Flat::LV`
        fOutput_Injected.Px = inj_view.Px();
        fOutput_Injected.Py = inj_view.Py();
        fOutput_Injected.Pz = inj_view.Pz();
        fOutput_Injected.Energy = Truth::Sexaquark::AsInjected_Energy(inj_view, fSettings.SexaquarkMass);
        // -- SV (`Flat::Coordinates`)
        fOutput_Injected.SV.X = inj_view.SV_X();
        fOutput_Injected.SV.Y = inj_view.SV_Y();
        fOutput_Injected.SV.Z = inj_view.SV_Z();
        // -- Nucleon (`Flat::LV`)
        fOutput_Injected.Nucleon.Px = inj_view.Nucleon_Px();
        fOutput_Injected.Nucleon.Py = inj_view.Nucleon_Py();
        fOutput_Injected.Nucleon.Pz = inj_view.Nucleon_Pz();
        fOutput_Injected.Nucleon.Energy = Truth::Sexaquark::AsInjected_NucleonEnergy(inj_view, fSettings.ReactionChannel);
        // -- event properties
        fOutput_Injected.RunNumber = fInput_Event.RunNumber;
        fOutput_Injected.DirNumber = fInput_Event.DirNumber;
        fOutput_Injected.EventNumber = fInput_Event.EventNumber;
        // -- reaction id
        fOutput_Injected.ReactionID = inj_view.ReactionID();
        fOutputTree_Injected->Fill();
    }
}

// ## Channel A ZONE ## //

void Finder::FindSexaquarks_ChannelA(bool anti_channel) {

    // determine properties based on anti-channel or not
    // -- v0a
    const Storage::Vector::V0s* Packed_V0A{&fInput_AntiLambdas};
    const Storage::Vector::MC_V0s* MC_V0A{&fInput_MC_AntiLambdas};
    EParticle pid_v0a_neg{EParticle::AntiProton};
    EParticle pid_v0a_pos{EParticle::PiPlus};
    if (anti_channel) {
        Packed_V0A = &fInput_Lambdas;
        MC_V0A = &fInput_MC_Lambdas;
        pid_v0a_neg = EParticle::PiMinus;
        pid_v0a_pos = EParticle::Proton;
    }
    double mass_v0a_neg{Const::Particle_Mass[pid_v0a_neg]};
    double mass_v0a_pos{Const::Particle_Mass[pid_v0a_pos]};
    // -- v0b
    const Storage::Vector::V0s* Packed_V0B{&fInput_KaonsZeroShort};
    const Storage::Vector::MC_V0s* MC_V0B{&fInput_MC_KaonsZeroShort};
    EParticle pid_v0b_neg{EParticle::PiMinus};
    EParticle pid_v0b_pos{EParticle::PiPlus};
    double mass_v0b_neg{Const::Particle_Mass[pid_v0b_neg]};
    double mass_v0b_pos{Const::Particle_Mass[pid_v0b_pos]};
    // -- cut flow hist
    TH1D* hist{anti_channel ? fHist_CutFlow_AntiChannel.get() : fHist_CutFlow.get()};

    // loop over all possible pairs of (anti)lambda + K0S //
    size_t n_v0a{Packed_V0A->X->size()};
    size_t n_v0b{Packed_V0B->X->size()};
    for (size_t v0a_entry{0}; v0a_entry < n_v0a; ++v0a_entry) {

        // unpack (anti)lambda //
        View::Rec::V0 v0a_view{Packed_V0A, static_cast<int>(v0a_entry)};
        Fit::V0 v0a{v0a_view, mass_v0a_neg, mass_v0a_pos};

        for (size_t v0b_entry{0}; v0b_entry < n_v0b; ++v0b_entry) {

            // unpack K0S //
            View::Rec::V0 v0b_view{Packed_V0B, static_cast<int>(v0b_entry)};
            Fit::V0 v0b{v0b_view, mass_v0b_neg, mass_v0b_pos};

            // sanity check //
            if (v0a_view.Neg.Index() == v0b_view.Neg.Index() || v0a_view.Neg.Index() == v0b_view.Pos.Index() ||
                v0a_view.Pos.Index() == v0b_view.Neg.Index() || v0a_view.Pos.Index() == v0b_view.Pos.Index()) {
                continue;
            }

            // fit //
            Fit::ChannelA sexa{v0a, v0b};
            sexa.DoFit(fInput_Event.MagneticField);

#ifdef T2S_DEBUG
            Logger::Debug(__FUNCTION__, "idx(v0a,neg,pos)={},{},{}", v0a_entry, v0a.Neg.View.Index(), v0a.Pos.View.Index());
            Logger::Debug(__FUNCTION__, ";x,y,z={},{},{}", v0a.X(), v0a.Y(), v0a.Z());
            Logger::Debug(__FUNCTION__, ";px,py,pz={},{},{}", v0a.Px(), v0a.Py(), v0a.Pz());
#ifdef T2S_LEGACY_KF
            Logger::Debug(__FUNCTION__, ";mass={}", v0a.GetMass());
#else
            Logger::Debug(__FUNCTION__, ";mass={}", v0a.Mass());
#endif
            Logger::Debug(__FUNCTION__, "idx(v0b,neg,pos)={},{},{}", v0b_entry, v0b.Neg.View.Index(), v0b.Pos.View.Index());
            Logger::Debug(__FUNCTION__, ";x,y,z={},{},{}", v0b.X(), v0b.Y(), v0b.Z());
            Logger::Debug(__FUNCTION__, ";px,py,pz={},{},{}", v0b.Px(), v0b.Py(), v0b.Pz());
#ifdef T2S_LEGACY_KF
            Logger::Debug(__FUNCTION__, ";mass={}", v0b.GetMass());
#else
            Logger::Debug(__FUNCTION__, ";mass={}", v0b.Mass());
#endif
            Logger::Debug(__FUNCTION__, "x,y,z={},{},{}", sexa.X(), sexa.Y(), sexa.Z());
// Logger::Debug(__FUNCTION__, ";x,y,z(v0a)={},{},{}", sexa.V0A_PCA_XYZ()[0], sexa.V0A_PCA_XYZ()[1], sexa.V0A_PCA_XYZ()[2]);
// Logger::Debug(__FUNCTION__, ";x,y,z(v0b)={},{},{}", sexa.V0B_PCA_XYZ()[0], sexa.V0B_PCA_XYZ()[1], sexa.V0B_PCA_XYZ()[2]);
#ifdef T2S_LEGACY_KF
            Logger::Debug(__FUNCTION__, ";mass={}", sexa.GetMass());
#else
            Logger::Debug(__FUNCTION__, ";mass={}", sexa.Mass());
#endif
            Logger::Debug(__FUNCTION__, ";mass_minus_n={}", sexa.Mass_MinusNucleon());
            Logger::Debug(__FUNCTION__, ";dca_btw_v0s={}", sexa.DCA_btw_V0s());
            // Logger::Debug(__FUNCTION__, ";radius={}", sexa.Radius2D()); // PENDING
            Logger::Debug(__FUNCTION__, ";dca_v0a={}", sexa.DCA_V0A_wrt_SV());
            Logger::Debug(__FUNCTION__, ";dca_v0b={}", sexa.DCA_V0B_wrt_SV());
#ifdef T2S_LEGACY_KF
            Logger::Debug(__FUNCTION__, ";pt={}", sexa.GetPt());
            Logger::Debug(__FUNCTION__, ";eta={}", sexa.GetEta());
#else
            Logger::Debug(__FUNCTION__, ";pt={}", sexa.Pt());
            Logger::Debug(__FUNCTION__, ";eta={}", sexa.Eta());
#endif
            Logger::Debug(__FUNCTION__, ";decay_length(v0a,v0b)={},{}", sexa.DecayLength_V0A(), sexa.DecayLength_V0B());
            // Logger::Debug(__FUNCTION__, ";cpa_pv={}", sexa.CPA_Point(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z)); // PENDING

            // PENDING: to rewrite... chore
            // if (IsMC()) {
            // View::PackedV0 mc_v0a_debug{*MC_V0A, v0a.Entry};
            // View::PackedV0 mc_v0b_debug{*MC_V0B, v0b.Entry};
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
                View::MC::ChannelA mc_sexa{&fInput_Injected, MC_V0A, MC_V0B, static_cast<int>(v0a_entry), static_cast<int>(v0b_entry)};
                StoreMC(mc_sexa);
            }
            fOutputTree->Fill();
        }
    }
}

bool Finder::PassesCuts(const Fit::ChannelA& sexa, TH1D* cut_flow_hist) const {

    cut_flow_hist->Fill(0.);

    // if (sexa.Radius2D() < Cuts::ChannelA::Min_Radius2D) return false; // PENDING
    // cut_flow_hist->Fill(1.); // PENDING

    if (sexa.DCA_V0A_wrt_SV() > Cuts::ChannelA::Max_DCALaSV) return false;
    cut_flow_hist->Fill(2.);

    if (sexa.DCA_V0B_wrt_SV() > Cuts::ChannelA::Max_DCAK0SV) return false;
    cut_flow_hist->Fill(3.);

    // if (sexa.CPA_wrt(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) < Cuts::ChannelA::Min_CPAwrtPV) return false;
    // cut_flow_hist->Fill(4.);

    // if (sexa.DCA_wrt(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) > Cuts::ChannelA::Max_DCAwrtPV) return false;
    // cut_flow_hist->Fill(5.);

    if (sexa.DCA_btw_V0s() > Cuts::ChannelA::Max_DCAbtwV0s) return false;
    cut_flow_hist->Fill(4.);

    if (sexa.CPA_V0A_wrt_SV() < Cuts::ChannelA::Min_La_CPAwrtSV) return false;
    cut_flow_hist->Fill(5.);

    if (sexa.CPA_V0B_wrt_SV() < Cuts::ChannelA::Min_K0S_CPAwrtSV) return false;
    cut_flow_hist->Fill(6.);

    return true;
}

void Finder::Store(const Fit::ChannelA& sexa, bool anti_channel) {
    // `Flat::ChannelA` //

    // `Flat::Sexaquark`
    // -- `Flat::State`
    fOutput_ChannelA.X = static_cast<float>(sexa.X());
    fOutput_ChannelA.Y = static_cast<float>(sexa.Y());
    fOutput_ChannelA.Z = static_cast<float>(sexa.Z());
    fOutput_ChannelA.Px = static_cast<float>(sexa.Px());
    fOutput_ChannelA.Py = static_cast<float>(sexa.Py());
    fOutput_ChannelA.Pz = static_cast<float>(sexa.Pz());
    fOutput_ChannelA.Energy = static_cast<float>(sexa.E());
    // -- event properties
    fOutput_ChannelA.Event = fInput_Event;
    // -- fit info
#ifdef T2S_LEGACY_KF
    fOutput_ChannelA.Chi2NDF = static_cast<float>(double(sexa.Chi2()) / double(sexa.NDF()));
#else
    fOutput_ChannelA.Chi2NDF = static_cast<float>(sexa.Chi2NDF());
#endif
    // -- extra info
    fOutput_ChannelA.E_MinusNucleon = static_cast<float>(sexa.E_MinusNucleon());
    fOutput_ChannelA.AntiChannel = anti_channel;

    // V0A (`Flat::V0`)
    // -- `Flat::State`
    fOutput_ChannelA.V0A.X = sexa.V0A.View.X();
    fOutput_ChannelA.V0A.Y = sexa.V0A.View.Y();
    fOutput_ChannelA.V0A.Z = sexa.V0A.View.Z();
    fOutput_ChannelA.V0A.Px = sexa.V0A.View.Px();
    fOutput_ChannelA.V0A.Py = sexa.V0A.View.Py();
    fOutput_ChannelA.V0A.Pz = sexa.V0A.View.Pz();
    fOutput_ChannelA.V0A.Energy = sexa.V0A.View.Energy();
    // -- Fit info
    fOutput_ChannelA.V0A.Chi2NDF = sexa.V0A.View.Chi2NDF();
    // -- @ PCA (`Flat::State_NoE`)
    fOutput_ChannelA.V0A_atPCA.X = static_cast<float>(sexa.V0A_PCA_XYZ()[0]);
    fOutput_ChannelA.V0A_atPCA.Y = static_cast<float>(sexa.V0A_PCA_XYZ()[1]);
    fOutput_ChannelA.V0A_atPCA.Z = static_cast<float>(sexa.V0A_PCA_XYZ()[2]);
    fOutput_ChannelA.V0A_atPCA.Px = static_cast<float>(sexa.V0A_PCA_PxPyPz()[0]);
    fOutput_ChannelA.V0A_atPCA.Py = static_cast<float>(sexa.V0A_PCA_PxPyPz()[1]);
    fOutput_ChannelA.V0A_atPCA.Pz = static_cast<float>(sexa.V0A_PCA_PxPyPz()[2]);

    // V0A Neg
    // -- `Flat::Track`
    fOutput_ChannelA.V0A.Neg.X = sexa.V0A.Neg.View.X();
    fOutput_ChannelA.V0A.Neg.Y = sexa.V0A.Neg.View.Y();
    fOutput_ChannelA.V0A.Neg.Z = sexa.V0A.Neg.View.Z();
    fOutput_ChannelA.V0A.Neg.Px = sexa.V0A.Neg.View.Px();
    fOutput_ChannelA.V0A.Neg.Py = sexa.V0A.Neg.View.Py();
    fOutput_ChannelA.V0A.Neg.Pz = sexa.V0A.Neg.View.Pz();
    fOutput_ChannelA.V0A.Neg.Index = sexa.V0A.Neg.View.Index();
    fOutput_ChannelA.V0A.Neg.DCAxy = sexa.V0A.Neg.View.DCAxy();
    fOutput_ChannelA.V0A.Neg.DCAz = sexa.V0A.Neg.View.DCAz();
    fOutput_ChannelA.V0A.Neg.TPCSignal = sexa.V0A.Neg.View.TPCSignal();
    fOutput_ChannelA.V0A.Neg.NSigmaPion = sexa.V0A.Neg.View.NSigmaPion();
    fOutput_ChannelA.V0A.Neg.NSigmaKaon = sexa.V0A.Neg.View.NSigmaKaon();
    fOutput_ChannelA.V0A.Neg.NSigmaProton = sexa.V0A.Neg.View.NSigmaProton();
    // -- @ V0 `Flat::State_NoE`
    fOutput_ChannelA.V0A.Neg_atV0.X = sexa.V0A.View.Neg_atV0_X();
    fOutput_ChannelA.V0A.Neg_atV0.Y = sexa.V0A.View.Neg_atV0_Y();
    fOutput_ChannelA.V0A.Neg_atV0.Z = sexa.V0A.View.Neg_atV0_Z();
    fOutput_ChannelA.V0A.Neg_atV0.Px = sexa.V0A.View.Neg_atV0_Px();
    fOutput_ChannelA.V0A.Neg_atV0.Py = sexa.V0A.View.Neg_atV0_Py();
    fOutput_ChannelA.V0A.Neg_atV0.Pz = sexa.V0A.View.Neg_atV0_Pz();

    // V0A Pos
    // -- `Flat::Track`
    fOutput_ChannelA.V0A.Pos.X = sexa.V0A.Pos.View.X();
    fOutput_ChannelA.V0A.Pos.Y = sexa.V0A.Pos.View.Y();
    fOutput_ChannelA.V0A.Pos.Z = sexa.V0A.Pos.View.Z();
    fOutput_ChannelA.V0A.Pos.Px = sexa.V0A.Pos.View.Px();
    fOutput_ChannelA.V0A.Pos.Py = sexa.V0A.Pos.View.Py();
    fOutput_ChannelA.V0A.Pos.Pz = sexa.V0A.Pos.View.Pz();
    fOutput_ChannelA.V0A.Pos.Index = sexa.V0A.Pos.View.Index();
    fOutput_ChannelA.V0A.Pos.DCAxy = sexa.V0A.Pos.View.DCAxy();
    fOutput_ChannelA.V0A.Pos.DCAz = sexa.V0A.Pos.View.DCAz();
    fOutput_ChannelA.V0A.Pos.TPCSignal = sexa.V0A.Pos.View.TPCSignal();
    fOutput_ChannelA.V0A.Pos.NSigmaPion = sexa.V0A.Pos.View.NSigmaPion();
    fOutput_ChannelA.V0A.Pos.NSigmaKaon = sexa.V0A.Pos.View.NSigmaKaon();
    fOutput_ChannelA.V0A.Pos.NSigmaProton = sexa.V0A.Pos.View.NSigmaProton();
    // -- @ V0 (`Flat::State_NoE`)
    fOutput_ChannelA.V0A.Pos_atV0.X = sexa.V0A.View.Pos_atV0_X();
    fOutput_ChannelA.V0A.Pos_atV0.Y = sexa.V0A.View.Pos_atV0_Y();
    fOutput_ChannelA.V0A.Pos_atV0.Z = sexa.V0A.View.Pos_atV0_Z();
    fOutput_ChannelA.V0A.Pos_atV0.Px = sexa.V0A.View.Pos_atV0_Px();
    fOutput_ChannelA.V0A.Pos_atV0.Py = sexa.V0A.View.Pos_atV0_Py();
    fOutput_ChannelA.V0A.Pos_atV0.Pz = sexa.V0A.View.Pos_atV0_Pz();

    // V0B
    // -- `Flat::State`
    fOutput_ChannelA.V0B.X = sexa.V0B.View.X();
    fOutput_ChannelA.V0B.Y = sexa.V0B.View.Y();
    fOutput_ChannelA.V0B.Z = sexa.V0B.View.Z();
    fOutput_ChannelA.V0B.Px = sexa.V0B.View.Px();
    fOutput_ChannelA.V0B.Py = sexa.V0B.View.Py();
    fOutput_ChannelA.V0B.Pz = sexa.V0B.View.Pz();
    fOutput_ChannelA.V0B.Energy = sexa.V0B.View.Energy();
    // -- Fit info
    fOutput_ChannelA.V0B.Chi2NDF = sexa.V0B.View.Chi2NDF();
    // -- @ PCA (`Flat::State_NoE`)
    fOutput_ChannelA.V0B_atPCA.X = static_cast<float>(sexa.V0B_PCA_XYZ()[0]);
    fOutput_ChannelA.V0B_atPCA.Y = static_cast<float>(sexa.V0B_PCA_XYZ()[1]);
    fOutput_ChannelA.V0B_atPCA.Z = static_cast<float>(sexa.V0B_PCA_XYZ()[2]);
    fOutput_ChannelA.V0B_atPCA.Px = static_cast<float>(sexa.V0B_PCA_PxPyPz()[0]);
    fOutput_ChannelA.V0B_atPCA.Py = static_cast<float>(sexa.V0B_PCA_PxPyPz()[1]);
    fOutput_ChannelA.V0B_atPCA.Pz = static_cast<float>(sexa.V0B_PCA_PxPyPz()[2]);

    // V0B Neg
    // -- `Flat::Track`
    fOutput_ChannelA.V0B.Neg.X = sexa.V0B.Neg.View.X();
    fOutput_ChannelA.V0B.Neg.Y = sexa.V0B.Neg.View.Y();
    fOutput_ChannelA.V0B.Neg.Z = sexa.V0B.Neg.View.Z();
    fOutput_ChannelA.V0B.Neg.Px = sexa.V0B.Neg.View.Px();
    fOutput_ChannelA.V0B.Neg.Py = sexa.V0B.Neg.View.Py();
    fOutput_ChannelA.V0B.Neg.Pz = sexa.V0B.Neg.View.Pz();
    fOutput_ChannelA.V0B.Neg.Index = sexa.V0B.Neg.View.Index();
    fOutput_ChannelA.V0B.Neg.DCAxy = sexa.V0B.Neg.View.DCAxy();
    fOutput_ChannelA.V0B.Neg.DCAz = sexa.V0B.Neg.View.DCAz();
    fOutput_ChannelA.V0B.Neg.TPCSignal = sexa.V0B.Neg.View.TPCSignal();
    fOutput_ChannelA.V0B.Neg.NSigmaPion = sexa.V0B.Neg.View.NSigmaPion();
    fOutput_ChannelA.V0B.Neg.NSigmaKaon = sexa.V0B.Neg.View.NSigmaKaon();
    fOutput_ChannelA.V0B.Neg.NSigmaProton = sexa.V0B.Neg.View.NSigmaProton();
    // -- @ V0 (`Flat::State_NoE`)
    fOutput_ChannelA.V0B.Neg_atV0.X = sexa.V0B.View.Neg_atV0_X();
    fOutput_ChannelA.V0B.Neg_atV0.Y = sexa.V0B.View.Neg_atV0_Y();
    fOutput_ChannelA.V0B.Neg_atV0.Z = sexa.V0B.View.Neg_atV0_Z();
    fOutput_ChannelA.V0B.Neg_atV0.Px = sexa.V0B.View.Neg_atV0_Px();
    fOutput_ChannelA.V0B.Neg_atV0.Py = sexa.V0B.View.Neg_atV0_Py();
    fOutput_ChannelA.V0B.Neg_atV0.Pz = sexa.V0B.View.Neg_atV0_Pz();

    // V0B Pos
    // -- `Flat::Track`
    fOutput_ChannelA.V0B.Pos.X = sexa.V0B.Pos.View.X();
    fOutput_ChannelA.V0B.Pos.Y = sexa.V0B.Pos.View.Y();
    fOutput_ChannelA.V0B.Pos.Z = sexa.V0B.Pos.View.Z();
    fOutput_ChannelA.V0B.Pos.Px = sexa.V0B.Pos.View.Px();
    fOutput_ChannelA.V0B.Pos.Py = sexa.V0B.Pos.View.Py();
    fOutput_ChannelA.V0B.Pos.Pz = sexa.V0B.Pos.View.Pz();
    fOutput_ChannelA.V0B.Pos.Index = sexa.V0B.Pos.View.Index();
    fOutput_ChannelA.V0B.Pos.DCAxy = sexa.V0B.Pos.View.DCAxy();
    fOutput_ChannelA.V0B.Pos.DCAz = sexa.V0B.Pos.View.DCAz();
    fOutput_ChannelA.V0B.Pos.TPCSignal = sexa.V0B.Pos.View.TPCSignal();
    fOutput_ChannelA.V0B.Pos.NSigmaPion = sexa.V0B.Pos.View.NSigmaPion();
    fOutput_ChannelA.V0B.Pos.NSigmaKaon = sexa.V0B.Pos.View.NSigmaKaon();
    fOutput_ChannelA.V0B.Pos.NSigmaProton = sexa.V0B.Pos.View.NSigmaProton();
    // -- @ V0 (`Flat::State_NoE`)
    fOutput_ChannelA.V0B.Pos_atV0.X = sexa.V0B.View.Pos_atV0_X();
    fOutput_ChannelA.V0B.Pos_atV0.Y = sexa.V0B.View.Pos_atV0_Y();
    fOutput_ChannelA.V0B.Pos_atV0.Z = sexa.V0B.View.Pos_atV0_Z();
    fOutput_ChannelA.V0B.Pos_atV0.Px = sexa.V0B.View.Pos_atV0_Px();
    fOutput_ChannelA.V0B.Pos_atV0.Py = sexa.V0B.View.Pos_atV0_Py();
    fOutput_ChannelA.V0B.Pos_atV0.Pz = sexa.V0B.View.Pos_atV0_Pz();
}

void Finder::StoreMC(const View::MC::ChannelA& view) {
    // `Flat::MC_ChannelA` //

    // `Flat::MC_Sexaquark`
    // -- Before (`Flat::LV`)
    fOutput_MC_ChannelA.Before.Px = view.Px();
    fOutput_MC_ChannelA.Before.Py = view.Py();
    fOutput_MC_ChannelA.Before.Pz = view.Pz();
    fOutput_MC_ChannelA.Before.Energy = Truth::Sexaquark::AsInjected_Energy(view, fSettings.SexaquarkMass);
    // -- After (`Flat::LV`)
    fOutput_MC_ChannelA.After.Px = static_cast<float>(Truth::Sexaquark::AsChannelA_AfterPx(view));
    fOutput_MC_ChannelA.After.Py = static_cast<float>(Truth::Sexaquark::AsChannelA_AfterPy(view));
    fOutput_MC_ChannelA.After.Pz = static_cast<float>(Truth::Sexaquark::AsChannelA_AfterPz(view));
    fOutput_MC_ChannelA.After.Energy = static_cast<float>(Truth::Sexaquark::AsChannelA_AfterE(view));
    // -- Nucleon (`Flat::LV`)
    fOutput_MC_ChannelA.Nucleon.Px = view.Nucleon_Px();
    fOutput_MC_ChannelA.Nucleon.Py = view.Nucleon_Py();
    fOutput_MC_ChannelA.Nucleon.Pz = view.Nucleon_Pz();
    fOutput_MC_ChannelA.Nucleon.Energy = Truth::Sexaquark::AsInjected_NucleonEnergy(view, fSettings.ReactionChannel);
    // -- PV (`Flat::Coordinates`)
    fOutput_MC_ChannelA.PV = fInput_MC_PV;
    // -- SV (`Flat::Coordinates`)
    fOutput_MC_ChannelA.SV.X = view.SV_X();
    fOutput_MC_ChannelA.SV.Y = view.SV_Y();
    fOutput_MC_ChannelA.SV.Z = view.SV_Z();
    // -- reaction id + flags
    fOutput_MC_ChannelA.ReactionID = view.ReactionID();
    fOutput_MC_ChannelA.IsSignal = Truth::Sexaquark::AsChannelA_IsSignal(view);
    fOutput_MC_ChannelA.IsHybrid = Truth::Sexaquark::AsChannelA_IsHybrid(view);

    // V0A (`Flat::MC_V0`)
    // -- `Flat::MC`
    fOutput_MC_ChannelA.V0A.McEntry = view.V0A.Entry;
    fOutput_MC_ChannelA.V0A.PdgCode = view.V0A.PdgCode();
    fOutput_MC_ChannelA.V0A.ReactionID = view.V0A.ReactionID();
    fOutput_MC_ChannelA.V0A.IsTrue = view.V0A.IsTrue();
    fOutput_MC_ChannelA.V0A.IsSignal = view.V0A.IsSignal();
    fOutput_MC_ChannelA.V0A.IsSecondary = view.V0A.IsSecondary();
    // -- `Flat::LV`
    fOutput_MC_ChannelA.V0A.Px = view.V0A.Px();
    fOutput_MC_ChannelA.V0A.Py = view.V0A.Py();
    fOutput_MC_ChannelA.V0A.Pz = view.V0A.Pz();
    fOutput_MC_ChannelA.V0A.Energy = view.V0A.Energy();
    // -- mother info (`Flat::MC_Id`)
    fOutput_MC_ChannelA.V0A.Mother.McEntry = view.V0A.Mother_Entry();
    fOutput_MC_ChannelA.V0A.Mother.PdgCode = view.V0A.Mother_PdgCode();
    // -- @ decay (`Flat::Coordinates`)
    fOutput_MC_ChannelA.V0A.AtDecay.X = view.V0A.Decay_X();
    fOutput_MC_ChannelA.V0A.AtDecay.Y = view.V0A.Decay_Y();
    fOutput_MC_ChannelA.V0A.AtDecay.Z = view.V0A.Decay_Z();
    // -- hybrid flag
    fOutput_MC_ChannelA.V0A.IsHybrid = view.V0A.IsHybrid();

    // V0A neg
    // -- `Flat::MC`
    fOutput_MC_ChannelA.V0A.Neg.McEntry = view.V0A.Neg_Entry();
    fOutput_MC_ChannelA.V0A.Neg.PdgCode = view.V0A.Neg_PdgCode();
    fOutput_MC_ChannelA.V0A.Neg.ReactionID = view.V0A.Neg_ReactionID();
    fOutput_MC_ChannelA.V0A.Neg.IsTrue = view.V0A.Neg_IsTrue();
    fOutput_MC_ChannelA.V0A.Neg.IsSignal = view.V0A.Neg_IsSignal();
    fOutput_MC_ChannelA.V0A.Neg.IsSecondary = view.V0A.Neg_IsSecondary();
    // -- `Flat::PxPyPz`
    fOutput_MC_ChannelA.V0A.Neg_Momentum.Px = view.V0A.Neg_Px();
    fOutput_MC_ChannelA.V0A.Neg_Momentum.Py = view.V0A.Neg_Py();
    fOutput_MC_ChannelA.V0A.Neg_Momentum.Pz = view.V0A.Neg_Pz();

    // V0A pos
    // -- `Flat::MC`
    fOutput_MC_ChannelA.V0A.Pos.McEntry = view.V0A.Pos_Entry();
    fOutput_MC_ChannelA.V0A.Pos.PdgCode = view.V0A.Pos_PdgCode();
    fOutput_MC_ChannelA.V0A.Pos.ReactionID = view.V0A.Pos_ReactionID();
    fOutput_MC_ChannelA.V0A.Pos.IsTrue = view.V0A.Pos_IsTrue();
    fOutput_MC_ChannelA.V0A.Pos.IsSignal = view.V0A.Pos_IsSignal();
    fOutput_MC_ChannelA.V0A.Pos.IsSecondary = view.V0A.Pos_IsSecondary();
    // -- `Flat::PxPyPz`
    fOutput_MC_ChannelA.V0A.Pos_Momentum.Px = view.V0A.Pos_Px();
    fOutput_MC_ChannelA.V0A.Pos_Momentum.Py = view.V0A.Pos_Py();
    fOutput_MC_ChannelA.V0A.Pos_Momentum.Pz = view.V0A.Pos_Pz();

    // V0B (`Flat::MC_V0`)
    // -- `Flat::MC`
    fOutput_MC_ChannelA.V0B.McEntry = view.V0B.Entry;
    fOutput_MC_ChannelA.V0B.PdgCode = view.V0B.PdgCode();
    fOutput_MC_ChannelA.V0B.ReactionID = view.V0B.ReactionID();
    fOutput_MC_ChannelA.V0B.IsTrue = view.V0B.IsTrue();
    fOutput_MC_ChannelA.V0B.IsSignal = view.V0B.IsSignal();
    fOutput_MC_ChannelA.V0B.IsSecondary = view.V0B.IsSecondary();
    // -- `Flat::LV`
    fOutput_MC_ChannelA.V0B.Px = view.V0B.Px();
    fOutput_MC_ChannelA.V0B.Py = view.V0B.Py();
    fOutput_MC_ChannelA.V0B.Pz = view.V0B.Pz();
    fOutput_MC_ChannelA.V0B.Energy = view.V0B.Energy();
    // -- mother info (`Flat::MC_Id`)
    fOutput_MC_ChannelA.V0B.Mother.McEntry = view.V0B.Mother_Entry();
    fOutput_MC_ChannelA.V0B.Mother.PdgCode = view.V0B.Mother_PdgCode();
    // -- @ decay (`Flat::Coordinates`)
    fOutput_MC_ChannelA.V0B.AtDecay.X = view.V0B.Decay_X();
    fOutput_MC_ChannelA.V0B.AtDecay.Y = view.V0B.Decay_Y();
    fOutput_MC_ChannelA.V0B.AtDecay.Z = view.V0B.Decay_Z();
    // -- hybrid flag
    fOutput_MC_ChannelA.V0B.IsHybrid = view.V0B.IsHybrid();

    // V0B neg
    // -- `Flat::MC`
    fOutput_MC_ChannelA.V0B.Neg.McEntry = view.V0B.Neg_Entry();
    fOutput_MC_ChannelA.V0B.Neg.PdgCode = view.V0B.Neg_PdgCode();
    fOutput_MC_ChannelA.V0B.Neg.ReactionID = view.V0B.Neg_ReactionID();
    fOutput_MC_ChannelA.V0B.Neg.IsTrue = view.V0B.Neg_IsTrue();
    fOutput_MC_ChannelA.V0B.Neg.IsSignal = view.V0B.Neg_IsSignal();
    fOutput_MC_ChannelA.V0B.Neg.IsSecondary = view.V0B.Neg_IsSecondary();
    // -- `Flat::PxPyPz`
    fOutput_MC_ChannelA.V0B.Neg_Momentum.Px = view.V0B.Neg_Px();
    fOutput_MC_ChannelA.V0B.Neg_Momentum.Py = view.V0B.Neg_Py();
    fOutput_MC_ChannelA.V0B.Neg_Momentum.Pz = view.V0B.Neg_Pz();

    // V0B pos
    // -- `Flat::MC`
    fOutput_MC_ChannelA.V0B.Pos.McEntry = view.V0B.Pos_Entry();
    fOutput_MC_ChannelA.V0B.Pos.PdgCode = view.V0B.Pos_PdgCode();
    fOutput_MC_ChannelA.V0B.Pos.ReactionID = view.V0B.Pos_ReactionID();
    fOutput_MC_ChannelA.V0B.Pos.IsTrue = view.V0B.Pos_IsTrue();
    fOutput_MC_ChannelA.V0B.Pos.IsSignal = view.V0B.Pos_IsSignal();
    fOutput_MC_ChannelA.V0B.Pos.IsSecondary = view.V0B.Pos_IsSecondary();
    // -- `Flat::PxPyPz`
    fOutput_MC_ChannelA.V0B.Pos_Momentum.Px = view.V0B.Pos_Px();
    fOutput_MC_ChannelA.V0B.Pos_Momentum.Py = view.V0B.Pos_Py();
    fOutput_MC_ChannelA.V0B.Pos_Momentum.Pz = view.V0B.Pos_Pz();
}

// ## Channel D ZONE ## //

void Finder::FindSexaquarks_ChannelD(bool anti_channel) {

    // determine properties based on anti-channel or not
    // -- v0
    const Storage::Vector::V0s* Packed_V0s{&fInput_AntiLambdas};
    const Storage::Vector::MC_V0s* MC_V0s{&fInput_MC_AntiLambdas};
    EParticle pid_v0_neg{EParticle::AntiProton};
    EParticle pid_v0_pos{EParticle::PiPlus};
    if (anti_channel) {
        Packed_V0s = &fInput_Lambdas;
        MC_V0s = &fInput_MC_Lambdas;
        pid_v0_neg = EParticle::PiMinus;
        pid_v0_pos = EParticle::Proton;
    }
    double mass_v0_neg{Const::Particle_Mass[pid_v0_neg]};
    double mass_v0_pos{Const::Particle_Mass[pid_v0_pos]};
    // -- kaon
    const Storage::Vector::Tracks* Packed_Kaons{&fInput_PosKaons};
    const Storage::Vector::MC_Tracks* MC_Kaons{&fInput_MC_PosKaons};
    EParticle pid_kaon{EParticle::PosKaon};
    if (anti_channel) {
        Packed_Kaons = &fInput_NegKaons;
        MC_Kaons = &fInput_MC_NegKaons;
        pid_kaon = EParticle::NegKaon;
    }
    double mass_kaon{Const::Particle_Mass[pid_kaon]};
    // -- cut flow hist
    TH1D* hist{anti_channel ? fHist_CutFlow_AntiChannel.get() : fHist_CutFlow.get()};

    // loop over all possible pairs of (anti)lambda + (pos/neg)kaon //
    size_t n_v0{Packed_V0s->X->size()};
    size_t n_kaons{Packed_Kaons->X->size()};
    for (size_t v0_entry{0}; v0_entry < n_v0; ++v0_entry) {

        // unpack (anti)lambda //
        View::Rec::V0 v0_view{Packed_V0s, static_cast<int>(v0_entry)};
        Fit::V0 v0{v0_view, mass_v0_neg, mass_v0_pos};

        for (size_t kaon_entry{0}; kaon_entry < n_kaons; ++kaon_entry) {

            // unpack (pos/neg)kaon //
            View::Rec::Track kaon_view{Packed_Kaons, static_cast<int>(kaon_entry)};
            Fit::Track kaon{kaon_view, mass_kaon};

            // sanity check //
            if (v0.Neg.View.Index() == kaon.View.Index() || v0.Pos.View.Index() == kaon.View.Index()) continue;

            // fit //
            Fit::ChannelD sexa{v0, kaon};
            sexa.DoFit(fInput_Event.MagneticField);

            // apply cuts //
            if (!PassesCuts(sexa, hist)) continue;
#ifdef T2S_DEBUG
            Logger::Debug(__FUNCTION__, "idx(v0,neg,pos,kaon)={},{},{},{}", v0_entry, sexa.V0.Neg.View.Index(), sexa.V0.Pos.View.Index(),
                          sexa.Kaon.View.Index());
            Logger::Debug(__FUNCTION__, ";x,y,z={},{},{}", sexa.X(), sexa.Y(), sexa.Z());
            Logger::Debug(__FUNCTION__, ";x,y,z(v0)={},{},{}", sexa.V0_PCA_XYZ()[0], sexa.V0_PCA_XYZ()[1], sexa.V0_PCA_XYZ()[2]);
            Logger::Debug(__FUNCTION__, ";x,y,z(kaon)={},{},{}", sexa.Kaon_PCA_XYZ()[0], sexa.Kaon_PCA_XYZ()[1], sexa.Kaon_PCA_XYZ()[2]);
#ifdef T2S_LEGACY_KF
            Logger::Debug(__FUNCTION__, ";mass={}", sexa.GetMass());
#else
            Logger::Debug(__FUNCTION__, ";mass={}", sexa.Mass());
#endif
            Logger::Debug(__FUNCTION__, ";dca_v0_kaon={}", sexa.DCA_btw_V0_Kaon());
            // Logger::Debug(__FUNCTION__, ";radius={}", sexa.Radius2D()); // PENDING
            Logger::Debug(__FUNCTION__, ";dca_v0={}", sexa.DCA_V0_wrt_SV());
            Logger::Debug(__FUNCTION__, ";dca_kaon={}", sexa.DCA_Kaon_wrt_SV());
#ifdef T2S_LEGACY_KF
            Logger::Debug(__FUNCTION__, ";pt={}", sexa.GetPt());
            Logger::Debug(__FUNCTION__, ";eta={}", sexa.GetEta());
#else
            Logger::Debug(__FUNCTION__, ";pt={}", sexa.Pt());
            Logger::Debug(__FUNCTION__, ";eta={}", sexa.Eta());
#endif
            // Logger::Debug(__FUNCTION__, ";cpa_pv={}", sexa.CPA_Point(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z)); // PENDING
#endif

            // store //
            Store(sexa, anti_channel);

            if (IsMC()) {
                View::MC::ChannelD mc_sexa{&fInput_Injected, MC_V0s, MC_Kaons, static_cast<int>(v0_entry), static_cast<int>(kaon_entry)};
                StoreMC(mc_sexa);
            }
            fOutputTree->Fill();
        }
    }
}

bool Finder::PassesCuts(const Fit::ChannelD& sexa, TH1D* cut_flow_hist) const {

    cut_flow_hist->Fill(0.);

    // if (sexa.Radius2D() < Cuts::ChannelD::Min_Radius2D || sexa.Radius2D() > Cuts::ChannelD::Max_Radius2D) return false; // PENDING
    // cut_flow_hist->Fill(1.); // PENDING

    // if (sexa.AbsRapidity_MinusNucleon() > Cuts::ChannelD::AbsMax_Rapidity) return false;  // PENDING: kinematics, affected by Fermi motion
    cut_flow_hist->Fill(2.);

    // if (sexa.CPA_Point(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) < Cuts::ChannelD::Min_CPAwrtPV ||
    // sexa.CPA_Point(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) > Cuts::ChannelD::Max_CPAwrtPV) {
    // return false;  // PENDING: kinematics, affected by Fermi motion
    // }
    cut_flow_hist->Fill(3.);

    if (sexa.DCA_V0_wrt_SV() > Cuts::ChannelD::Max_DCALaSV) return false;
    cut_flow_hist->Fill(4.);

    if (sexa.DCA_Kaon_wrt_SV() > Cuts::ChannelD::Max_DCAKaSV) return false;
    cut_flow_hist->Fill(5.);

    // if (sexa.DCA_V0Neg_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelD::Max_DCALaNegSV) return false;
    cut_flow_hist->Fill(6.);

    // if (sexa.DCA_V0Pos_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelD::Max_DCALaPosSV) return false;
    cut_flow_hist->Fill(7.);

    if (sexa.DCA_btw_V0_Kaon() > Cuts::ChannelD::Max_DCAKaLa) return false;
    cut_flow_hist->Fill(8.);

    return true;
}

void Finder::Store(const Fit::ChannelD& sexa, bool anti_channel) {
    // `Flat::ChannelD` //

    // `Flat::Sexaquark`
    // -- `Flat::State`
    fOutput_ChannelD.X = static_cast<float>(sexa.X());
    fOutput_ChannelD.Y = static_cast<float>(sexa.Y());
    fOutput_ChannelD.Z = static_cast<float>(sexa.Z());
    fOutput_ChannelD.Px = static_cast<float>(sexa.Px());
    fOutput_ChannelD.Py = static_cast<float>(sexa.Py());
    fOutput_ChannelD.Pz = static_cast<float>(sexa.Pz());
    fOutput_ChannelD.Energy = static_cast<float>(sexa.E());
    // -- `Flat::Event`
    fOutput_ChannelD.Event = fInput_Event;
// -- fit info
#ifdef T2S_LEGACY_KF
    fOutput_ChannelD.Chi2NDF = static_cast<float>(double(sexa.Chi2()) / double(sexa.NDF()));
#else
    fOutput_ChannelD.Chi2NDF = static_cast<float>(sexa.Chi2NDF());
#endif
    // -- extra info
    fOutput_ChannelD.E_MinusNucleon = static_cast<float>(sexa.E_MinusNucleon());
    fOutput_ChannelD.AntiChannel = anti_channel;

    // V0 (`Flat::V0`)
    // -- `Flat::State`
    fOutput_ChannelD.V0.X = sexa.V0.View.X();
    fOutput_ChannelD.V0.Y = sexa.V0.View.Y();
    fOutput_ChannelD.V0.Z = sexa.V0.View.Z();
    fOutput_ChannelD.V0.Px = sexa.V0.View.Px();
    fOutput_ChannelD.V0.Py = sexa.V0.View.Py();
    fOutput_ChannelD.V0.Pz = sexa.V0.View.Pz();
    fOutput_ChannelD.V0.Energy = sexa.V0.View.Energy();
    // -- fit info
    fOutput_ChannelD.V0.Chi2NDF = sexa.V0.View.Chi2NDF();
    // -- @ PCA (`Flat::State_NoE`)
    fOutput_ChannelD.V0_atPCA.X = static_cast<float>(sexa.V0_PCA_XYZ()[0]);
    fOutput_ChannelD.V0_atPCA.Y = static_cast<float>(sexa.V0_PCA_XYZ()[1]);
    fOutput_ChannelD.V0_atPCA.Z = static_cast<float>(sexa.V0_PCA_XYZ()[2]);
    fOutput_ChannelD.V0_atPCA.Px = static_cast<float>(sexa.V0_PCA_PxPyPz()[0]);
    fOutput_ChannelD.V0_atPCA.Py = static_cast<float>(sexa.V0_PCA_PxPyPz()[1]);
    fOutput_ChannelD.V0_atPCA.Pz = static_cast<float>(sexa.V0_PCA_PxPyPz()[2]);

    // V0 Neg (`Flat::Track`)
    // -- `Flat::State_NoE`
    fOutput_ChannelD.V0.Neg.X = sexa.V0.Neg.View.X();
    fOutput_ChannelD.V0.Neg.Y = sexa.V0.Neg.View.Y();
    fOutput_ChannelD.V0.Neg.Z = sexa.V0.Neg.View.Z();
    fOutput_ChannelD.V0.Neg.Px = sexa.V0.Neg.View.Px();
    fOutput_ChannelD.V0.Neg.Py = sexa.V0.Neg.View.Py();
    fOutput_ChannelD.V0.Neg.Pz = sexa.V0.Neg.View.Pz();
    // -- Rest of info
    fOutput_ChannelD.V0.Neg.Index = sexa.V0.Neg.View.Index();
    fOutput_ChannelD.V0.Neg.DCAxy = sexa.V0.Neg.View.DCAxy();
    fOutput_ChannelD.V0.Neg.DCAz = sexa.V0.Neg.View.DCAz();
    fOutput_ChannelD.V0.Neg.TPCSignal = sexa.V0.Neg.View.TPCSignal();
    fOutput_ChannelD.V0.Neg.NSigmaPion = sexa.V0.Neg.View.NSigmaPion();
    fOutput_ChannelD.V0.Neg.NSigmaKaon = sexa.V0.Neg.View.NSigmaKaon();
    fOutput_ChannelD.V0.Neg.NSigmaProton = sexa.V0.Neg.View.NSigmaProton();

    // V0 Pos (`Flat::Track`)
    // -- `Flat::State_NoE`
    fOutput_ChannelD.V0.Pos.X = sexa.V0.Pos.View.X();
    fOutput_ChannelD.V0.Pos.Y = sexa.V0.Pos.View.Y();
    fOutput_ChannelD.V0.Pos.Z = sexa.V0.Pos.View.Z();
    fOutput_ChannelD.V0.Pos.Px = sexa.V0.Pos.View.Px();
    fOutput_ChannelD.V0.Pos.Py = sexa.V0.Pos.View.Py();
    fOutput_ChannelD.V0.Pos.Pz = sexa.V0.Pos.View.Pz();
    // -- Rest of info
    fOutput_ChannelD.V0.Pos.Index = sexa.V0.Pos.View.Index();
    fOutput_ChannelD.V0.Pos.DCAxy = sexa.V0.Pos.View.DCAxy();
    fOutput_ChannelD.V0.Pos.DCAz = sexa.V0.Pos.View.DCAz();
    fOutput_ChannelD.V0.Pos.TPCSignal = sexa.V0.Pos.View.TPCSignal();
    fOutput_ChannelD.V0.Pos.NSigmaPion = sexa.V0.Pos.View.NSigmaPion();
    fOutput_ChannelD.V0.Pos.NSigmaKaon = sexa.V0.Pos.View.NSigmaKaon();
    fOutput_ChannelD.V0.Pos.NSigmaProton = sexa.V0.Pos.View.NSigmaProton();

    // Kaon (`Flat::Track`)
    // -- `Flat::State_NoE`
    fOutput_ChannelD.Kaon.X = sexa.Kaon.View.X();
    fOutput_ChannelD.Kaon.Y = sexa.Kaon.View.Y();
    fOutput_ChannelD.Kaon.Z = sexa.Kaon.View.Z();
    fOutput_ChannelD.Kaon.Px = sexa.Kaon.View.Px();
    fOutput_ChannelD.Kaon.Py = sexa.Kaon.View.Py();
    fOutput_ChannelD.Kaon.Pz = sexa.Kaon.View.Pz();
    // -- Rest of info
    fOutput_ChannelD.Kaon.Index = sexa.Kaon.View.Index();
    fOutput_ChannelD.Kaon.DCAxy = sexa.Kaon.View.DCAxy();
    fOutput_ChannelD.Kaon.DCAz = sexa.Kaon.View.DCAz();
    fOutput_ChannelD.Kaon.TPCSignal = sexa.Kaon.View.TPCSignal();
    fOutput_ChannelD.Kaon.NSigmaPion = sexa.Kaon.View.NSigmaPion();
    fOutput_ChannelD.Kaon.NSigmaKaon = sexa.Kaon.View.NSigmaKaon();
    fOutput_ChannelD.Kaon.NSigmaProton = sexa.Kaon.View.NSigmaProton();
    // -- @ PCA (`Flat::State_NoE`)
    fOutput_ChannelD.Kaon_atPCA.X = static_cast<float>(sexa.Kaon_PCA_XYZ()[0]);
    fOutput_ChannelD.Kaon_atPCA.Y = static_cast<float>(sexa.Kaon_PCA_XYZ()[1]);
    fOutput_ChannelD.Kaon_atPCA.Z = static_cast<float>(sexa.Kaon_PCA_XYZ()[2]);
    fOutput_ChannelD.Kaon_atPCA.Px = static_cast<float>(sexa.Kaon_PCA_PxPyPz()[0]);
    fOutput_ChannelD.Kaon_atPCA.Py = static_cast<float>(sexa.Kaon_PCA_PxPyPz()[1]);
    fOutput_ChannelD.Kaon_atPCA.Pz = static_cast<float>(sexa.Kaon_PCA_PxPyPz()[2]);
}

void Finder::StoreMC(const View::MC::ChannelD& view) {
    // `Flat::MC_ChannelD` //

    // `Flat::MC_Sexaquark`
    // -- Before (`Flat::LV`)
    fOutput_MC_ChannelD.Before.Px = view.Px();
    fOutput_MC_ChannelD.Before.Py = view.Py();
    fOutput_MC_ChannelD.Before.Pz = view.Pz();
    fOutput_MC_ChannelD.Before.Energy = Truth::Sexaquark::AsInjected_Energy(view, fSettings.SexaquarkMass);
    // -- After (`Flat::LV`)
    fOutput_MC_ChannelD.After.Px = static_cast<float>(Truth::Sexaquark::AsChannelD_AfterPx(view));
    fOutput_MC_ChannelD.After.Py = static_cast<float>(Truth::Sexaquark::AsChannelD_AfterPy(view));
    fOutput_MC_ChannelD.After.Pz = static_cast<float>(Truth::Sexaquark::AsChannelD_AfterPz(view));
    fOutput_MC_ChannelD.After.Energy = static_cast<float>(Truth::Sexaquark::AsChannelD_AfterE(view));
    // -- Nucleon (`Flat::LV`)
    fOutput_MC_ChannelD.Nucleon.Px = view.Nucleon_Px();
    fOutput_MC_ChannelD.Nucleon.Py = view.Nucleon_Py();
    fOutput_MC_ChannelD.Nucleon.Pz = view.Nucleon_Pz();
    fOutput_MC_ChannelD.Nucleon.Energy = Truth::Sexaquark::AsInjected_NucleonEnergy(view, fSettings.ReactionChannel);
    // -- PV (`Flat::Coordinates`)
    fOutput_MC_ChannelD.PV = fInput_MC_PV;
    // -- SV (`Flat::Coordinates`)
    fOutput_MC_ChannelD.SV.X = view.SV_X();
    fOutput_MC_ChannelD.SV.Y = view.SV_Y();
    fOutput_MC_ChannelD.SV.Z = view.SV_Z();
    // -- Reaction ID + Flags
    fOutput_MC_ChannelD.ReactionID = view.ReactionID();
    fOutput_MC_ChannelD.IsSignal = Truth::Sexaquark::AsChannelD_IsSignal(view);
    fOutput_MC_ChannelD.IsHybrid = Truth::Sexaquark::AsChannelD_IsHybrid(view);

    // V0 (`Flat::MC_V0`)
    // -- `Flat::MC`
    fOutput_MC_ChannelD.V0.McEntry = view.V0.Entry;
    fOutput_MC_ChannelD.V0.PdgCode = view.V0.PdgCode();
    fOutput_MC_ChannelD.V0.ReactionID = view.V0.ReactionID();
    fOutput_MC_ChannelD.V0.IsTrue = view.V0.IsTrue();
    fOutput_MC_ChannelD.V0.IsSignal = view.V0.IsSignal();
    fOutput_MC_ChannelD.V0.IsSecondary = view.V0.IsSecondary();
    // -- `Flat::LV`
    fOutput_MC_ChannelD.V0.Px = view.V0.Px();
    fOutput_MC_ChannelD.V0.Py = view.V0.Py();
    fOutput_MC_ChannelD.V0.Pz = view.V0.Pz();
    fOutput_MC_ChannelD.V0.Energy = view.V0.Energy();
    // -- Mother (`Flat::MC_Id`)
    fOutput_MC_ChannelD.V0.Mother.McEntry = view.V0.Mother_Entry();
    fOutput_MC_ChannelD.V0.Mother.PdgCode = view.V0.Mother_PdgCode();
    // -- @ decay (`Flat::Coordinates`)
    fOutput_MC_ChannelD.V0.AtDecay.X = view.V0.Decay_X();
    fOutput_MC_ChannelD.V0.AtDecay.Y = view.V0.Decay_Y();
    fOutput_MC_ChannelD.V0.AtDecay.Z = view.V0.Decay_Z();
    // -- hybrid flag
    fOutput_MC_ChannelD.V0.IsHybrid = view.V0.IsHybrid();

    // V0 Neg
    // -- `Flat::MC`
    fOutput_MC_ChannelD.V0.Neg.McEntry = view.V0.Neg_Entry();
    fOutput_MC_ChannelD.V0.Neg.PdgCode = view.V0.Neg_PdgCode();
    fOutput_MC_ChannelD.V0.Neg.ReactionID = view.V0.Neg_ReactionID();
    fOutput_MC_ChannelD.V0.Neg.IsTrue = view.V0.Neg_IsTrue();
    fOutput_MC_ChannelD.V0.Neg.IsSignal = view.V0.Neg_IsSignal();
    fOutput_MC_ChannelD.V0.Neg.IsSecondary = view.V0.Neg_IsSecondary();
    // -- `Flat::PxPyPz`
    fOutput_MC_ChannelD.V0.Neg_Momentum.Px = view.V0.Neg_Px();
    fOutput_MC_ChannelD.V0.Neg_Momentum.Py = view.V0.Neg_Py();
    fOutput_MC_ChannelD.V0.Neg_Momentum.Pz = view.V0.Neg_Pz();

    // V0 Pos
    // -- `Flat::MC`
    fOutput_MC_ChannelD.V0.Pos.McEntry = view.V0.Pos_Entry();
    fOutput_MC_ChannelD.V0.Pos.PdgCode = view.V0.Pos_PdgCode();
    fOutput_MC_ChannelD.V0.Pos.ReactionID = view.V0.Pos_ReactionID();
    fOutput_MC_ChannelD.V0.Pos.IsTrue = view.V0.Pos_IsTrue();
    fOutput_MC_ChannelD.V0.Pos.IsSignal = view.V0.Pos_IsSignal();
    fOutput_MC_ChannelD.V0.Pos.IsSecondary = view.V0.Pos_IsSecondary();
    // -- `Flat::PxPyPz`
    fOutput_MC_ChannelD.V0.Pos_Momentum.Px = view.V0.Pos_Px();
    fOutput_MC_ChannelD.V0.Pos_Momentum.Py = view.V0.Pos_Py();
    fOutput_MC_ChannelD.V0.Pos_Momentum.Pz = view.V0.Pos_Pz();

    // Kaon (`Flat::MC_Track`)
    // -- `Flat::LV`
    fOutput_MC_ChannelD.Kaon.Px = view.Kaon.Px();
    fOutput_MC_ChannelD.Kaon.Py = view.Kaon.Py();
    fOutput_MC_ChannelD.Kaon.Pz = view.Kaon.Pz();
    fOutput_MC_ChannelD.Kaon.Energy = view.Kaon.Energy();
    // -- `Flat::MC`
    fOutput_MC_ChannelD.Kaon.McEntry = view.Kaon.Entry;
    fOutput_MC_ChannelD.Kaon.PdgCode = view.Kaon.PdgCode();
    fOutput_MC_ChannelD.Kaon.ReactionID = view.Kaon.ReactionID();
    fOutput_MC_ChannelD.Kaon.IsTrue = view.Kaon.IsTrue();
    fOutput_MC_ChannelD.Kaon.IsSignal = view.Kaon.IsSignal();
    fOutput_MC_ChannelD.Kaon.IsSecondary = view.Kaon.IsSecondary();
    // -- Mother (`Flat::MC_Id`)
    fOutput_MC_ChannelD.Kaon.Mother.McEntry = view.Kaon.Mother_Entry();
    fOutput_MC_ChannelD.Kaon.Mother.PdgCode = view.Kaon.Mother_PdgCode();
    // -- GrandMother (`Flat::MC_Id`)
    fOutput_MC_ChannelD.Kaon.GrandMother.McEntry = view.Kaon.GrandMother_Entry();
    fOutput_MC_ChannelD.Kaon.GrandMother.PdgCode = view.Kaon.GrandMother_PdgCode();
}

// ## END OF CYCLES ## //

void Finder::EndOfAnalysis() {

    Logger::Info(__FUNCTION__, "The following objects have been written into TFile {}:", fSettings.PathOutputFile);

    // write trees
    if (IsMC()) {
        fOutputTree_Injected->Write();
        Logger::Info(__FUNCTION__, "- TTree \"{}\"", fOutputTree_Injected->GetName());
    }
    fOutputTree->Write();
    Logger::Info(__FUNCTION__, "- TTree \"{}\"", fOutputTree->GetName());

    // write histograms
    // -- event counter
    fHist_EventCounter->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_EventCounter->GetName());
    // -- cut flows
    fHist_CutFlow->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow->GetName());
    fHist_CutFlow_AntiChannel->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow_AntiChannel->GetName());

    fInputChain_PackedEvents->ResetBranchAddresses();
    fOutputTree->ResetBranchAddresses();

    Logger::Info(__FUNCTION__, "All done.");
}

}  // namespace Tree2Secondaries
