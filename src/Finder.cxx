#include <filesystem>
#include <memory>

#include "App/Logger.hxx"
#include "Finder/Finder.hxx"
#include "Finder/FinderCuts.hxx"
#include "KalmanFitter/BaseKalmanFitter.hxx"
#include "KalmanFitter/KalmanFitterChannelA.hxx"
#include "KalmanFitter/KalmanFitterChannelD.hxx"
#include "Math/BaseMath.hxx"
#include "Math/Constants.hxx"
#include "Seeder/SeederHelixLine.hxx"
#include "Seeder/SeederLineLine.hxx"
#include "Storage/Vector/VectorTracks.hxx"
#include "Storage/Vector/VectorV0s.hxx"
#include "Truth/TruthSexaquark.hxx"
#include "View/MC/ViewMcInjected.hxx"
#include "View/MC/ViewMcPackedTrack.hxx"
#include "View/MC/ViewMcPackedV0.hxx"
#include "View/Reconstructed/ViewRecTrack.hxx"

namespace Tree2Secondaries {

namespace KF = KalmanFitter;

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
            fInput_AntiLambdas.ReadBranches_VectorV0s(fInputChain_PackedEvents.get(), Const::V0_Acronym[PID_V0::AntiLambda]);
            fInput_Lambdas.ReadBranches_VectorV0s(fInputChain_PackedEvents.get(), Const::V0_Acronym[PID_V0::Lambda]);
            fInput_KaonsZeroShort.ReadBranches_VectorV0s(fInputChain_PackedEvents.get(), Const::V0_Acronym[PID_V0::KaonZeroShort]);
            if (IsMC()) {
                fInput_MC_AntiLambdas.ReadBranches_VectorMC_V0s(fInputChain_PackedEvents.get(), Const::V0_Acronym[PID_V0::AntiLambda]);
                fInput_MC_Lambdas.ReadBranches_VectorMC_V0s(fInputChain_PackedEvents.get(), Const::V0_Acronym[PID_V0::Lambda]);
                fInput_MC_KaonsZeroShort.ReadBranches_VectorMC_V0s(fInputChain_PackedEvents.get(), Const::V0_Acronym[PID_V0::KaonZeroShort]);
            }
            break;
        case EReactionChannel::D:
            fInput_AntiLambdas.ReadBranches_VectorV0s(fInputChain_PackedEvents.get(), Const::V0_Acronym[PID_V0::AntiLambda]);
            fInput_Lambdas.ReadBranches_VectorV0s(fInputChain_PackedEvents.get(), Const::V0_Acronym[PID_V0::Lambda]);
            fInput_NegKaons.ReadBranches_VectorTracks(fInputChain_PackedEvents.get(), IsMC(), true,
                                                      Const::Particle_Acronym[PID_StableParticle::NegKaon]);
            fInput_PosKaons.ReadBranches_VectorTracks(fInputChain_PackedEvents.get(), IsMC(), true,
                                                      Const::Particle_Acronym[PID_StableParticle::PosKaon]);
            if (IsMC()) {
                fInput_MC_AntiLambdas.ReadBranches_VectorMC_V0s(fInputChain_PackedEvents.get(), Const::V0_Acronym[PID_V0::AntiLambda]);
                fInput_MC_Lambdas.ReadBranches_VectorMC_V0s(fInputChain_PackedEvents.get(), Const::V0_Acronym[PID_V0::Lambda]);
                fInput_MC_NegKaons.ReadBranches_VectorMC_Tracks(fInputChain_PackedEvents.get(), Const::Particle_Acronym[PID_StableParticle::NegKaon]);
                fInput_MC_PosKaons.ReadBranches_VectorMC_Tracks(fInputChain_PackedEvents.get(), Const::Particle_Acronym[PID_StableParticle::PosKaon]);
            }
            break;
        case EReactionChannel::E:
            fInput_AntiLambdas.ReadBranches_VectorV0s(fInputChain_PackedEvents.get(), Const::V0_Acronym[PID_V0::AntiLambda]);
            fInput_Lambdas.ReadBranches_VectorV0s(fInputChain_PackedEvents.get(), Const::V0_Acronym[PID_V0::Lambda]);
            fInput_NegKaons.ReadBranches_VectorTracks(fInputChain_PackedEvents.get(), IsMC(), true,
                                                      Const::Particle_Acronym[PID_StableParticle::NegKaon]);
            fInput_PosKaons.ReadBranches_VectorTracks(fInputChain_PackedEvents.get(), IsMC(), true,
                                                      Const::Particle_Acronym[PID_StableParticle::PosKaon]);
            fInput_PiMinus.ReadBranches_VectorTracks(fInputChain_PackedEvents.get(), IsMC(), true,
                                                     Const::Particle_Acronym[PID_StableParticle::PiMinus]);
            fInput_PiPlus.ReadBranches_VectorTracks(fInputChain_PackedEvents.get(), IsMC(), true,
                                                    Const::Particle_Acronym[PID_StableParticle::PiPlus]);
            if (IsMC()) {
                fInput_MC_AntiLambdas.ReadBranches_VectorMC_V0s(fInputChain_PackedEvents.get(), Const::V0_Acronym[PID_V0::AntiLambda]);
                fInput_MC_Lambdas.ReadBranches_VectorMC_V0s(fInputChain_PackedEvents.get(), Const::V0_Acronym[PID_V0::Lambda]);
                fInput_MC_NegKaons.ReadBranches_VectorMC_Tracks(fInputChain_PackedEvents.get(), Const::Particle_Acronym[PID_StableParticle::NegKaon]);
                fInput_MC_PosKaons.ReadBranches_VectorMC_Tracks(fInputChain_PackedEvents.get(), Const::Particle_Acronym[PID_StableParticle::PosKaon]);
                fInput_MC_PiMinus.ReadBranches_VectorMC_Tracks(fInputChain_PackedEvents.get(), Const::Particle_Acronym[PID_StableParticle::PiMinus]);
                fInput_MC_PiPlus.ReadBranches_VectorMC_Tracks(fInputChain_PackedEvents.get(), Const::Particle_Acronym[PID_StableParticle::PiPlus]);
            }
            break;
        case EReactionChannel::H:
            fInput_NegKaons.ReadBranches_VectorTracks(fInputChain_PackedEvents.get(), IsMC(), true,
                                                      Const::Particle_Acronym[PID_StableParticle::NegKaon]);
            fInput_PosKaons.ReadBranches_VectorTracks(fInputChain_PackedEvents.get(), IsMC(), true,
                                                      Const::Particle_Acronym[PID_StableParticle::PosKaon]);
            if (IsMC()) {
                fInput_MC_NegKaons.ReadBranches_VectorMC_Tracks(fInputChain_PackedEvents.get(), Const::Particle_Acronym[PID_StableParticle::NegKaon]);
                fInput_MC_PosKaons.ReadBranches_VectorMC_Tracks(fInputChain_PackedEvents.get(), Const::Particle_Acronym[PID_StableParticle::PosKaon]);
            }
            break;
        default:
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
    std::string hist_title = ";Cut N;N Passed Cut";
    fHist_CutFlow = std::make_unique<TH1D>("CutFlow", hist_title.c_str(), x_nbins, x_min, x_max);
    fHist_CutFlow_AntiChannel = std::make_unique<TH1D>("CutFlow_Anti", hist_title.c_str(), x_nbins, x_min, x_max);
}

bool Finder::PrepareOutputTree() {

    std::string tree_name = std::format("FoundChannel{}", Const::ReactionChannel_Char[fSettings.ReactionChannel]);

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

    fOutputTree_Injected = std::make_unique<TTree>(Const::TreeName_Injected.c_str(), "");
    if (fOutputTree_Injected == nullptr) {
        Logger::Error(__FUNCTION__, "Couldn't create TTree \"{}\"", Const::TreeName_Injected);
        return false;
    }

    return true;
}

void Finder::ProcessInjected() {

    double nucleon_mass = Const::Particle_Mass[Const::ReactionChannel_NucleonPID[fSettings.ReactionChannel]];

    auto n_injected = static_cast<int>(NumberInjected());
    for (int entry = 0; entry < n_injected; ++entry) {

        View::MC::Injected sexa{&fInput_Injected, entry};  // NOTE: unguarded because contiguous access

        // `Flat::Injected` //

        // -- `Flat::LV`
        fOutput_Injected.Px = sexa.Px();
        fOutput_Injected.Py = sexa.Py();
        fOutput_Injected.Pz = sexa.Pz();
        fOutput_Injected.Energy = static_cast<float>(Truth::Sexaquark::Energy(sexa, fSettings.SexaquarkMass));
        // -- SV (`Flat::Coordinates`)
        fOutput_Injected.SV.X = sexa.SV_X();
        fOutput_Injected.SV.Y = sexa.SV_Y();
        fOutput_Injected.SV.Z = sexa.SV_Z();
        // -- Nucleon (`Flat::LV`)
        fOutput_Injected.Nucleon.Px = sexa.Nucleon_Px();
        fOutput_Injected.Nucleon.Py = sexa.Nucleon_Py();
        fOutput_Injected.Nucleon.Pz = sexa.Nucleon_Pz();
        fOutput_Injected.Nucleon.Energy = static_cast<float>(Truth::Sexaquark::NucleonEnergy(sexa, nucleon_mass));
        // -- event properties
        fOutput_Injected.RunNumber = fInput_Event.RunNumber;
        fOutput_Injected.DirNumber = fInput_Event.DirNumber;
        fOutput_Injected.EventNumber = fInput_Event.EventNumber;
        // -- reaction id
        fOutput_Injected.ReactionID = Truth::Sexaquark::ReactionID(sexa);

        // fill tree //

        fOutputTree_Injected->Fill();
    }
}

// ## Channel A ZONE ## //

void Finder::FindSexaquarks_ChannelA(bool anti_channel) {

    // determine properties based on anti-channel or not
    // -- v0a
    const Storage::Vector::V0s* Packed_V0A{&fInput_AntiLambdas};
    const Storage::Vector::MC_V0s* MC_V0A{&fInput_MC_AntiLambdas};
    if (anti_channel) {
        Packed_V0A = &fInput_Lambdas;
        MC_V0A = &fInput_MC_Lambdas;
    }
    // -- v0b
    const Storage::Vector::V0s* Packed_V0B = &fInput_KaonsZeroShort;
    const Storage::Vector::MC_V0s* MC_V0B = &fInput_MC_KaonsZeroShort;
    // -- cut flow hist
    TH1D* hist = anti_channel ? fHist_CutFlow_AntiChannel.get() : fHist_CutFlow.get();

    // loop over all possible pairs of (anti)lambda + K0S //
    size_t n_v0a = Packed_V0A->X->size();
    size_t n_v0b = Packed_V0B->X->size();
    for (size_t v0a_entry = 0; v0a_entry < n_v0a; ++v0a_entry) {
        for (size_t v0b_entry = 0; v0b_entry < n_v0b; ++v0b_entry) {

            // get views //
            View::Rec::V0 v0a{Packed_V0A, v0a_entry};
            View::Rec::V0 v0b{Packed_V0B, v0b_entry};

            // sanity check //
            if (v0a.Neg.Index() == v0b.Neg.Index() || v0a.Neg.Index() == v0b.Pos.Index() || v0a.Pos.Index() == v0b.Neg.Index() ||
                v0a.Pos.Index() == v0b.Pos.Index()) {
                continue;
            }

            // PCAs //
            Seeder::LineLine::Cache pca_cache;
            auto [seed_v0a, seed_v0b] = Seeder::LineLine::FastPCAs(v0a, v0b, &pca_cache);

            // apply cuts (1) //
            if (!FastCuts_ChannelA(seed_v0a, seed_v0b, hist)) continue;

            // PCAs derivatives //
            auto [deriv_v0a, deriv_v0b] = Seeder::LineLine::ComputeDerivatives(pca_cache);

            // fit vertex //
            auto fit = KF::FitVertex(v0a, v0b, {seed_v0a, deriv_v0a}, {seed_v0b, deriv_v0b});
            KF::ChannelA sexa(fit, seed_v0a.pca, seed_v0b.pca, v0a, v0b);

            // apply cuts (2) //
            if (!SlowCuts(sexa, hist)) continue;

            // store //
            Store(sexa, anti_channel);

            if (IsMC()) {
                View::MC::PackedV0 mc_v0a(MC_V0A, v0a_entry);
                View::MC::PackedV0 mc_v0b(MC_V0B, v0b_entry);
                View::MC::ChannelA mc_sexa(&fInput_Injected, mc_v0a, mc_v0b);
                if (View::IsValid(mc_sexa)) {
                    StoreMC(mc_sexa);
                } else {
                    StoreDummyMC_ChannelA();
                }
            }

            // fill flat tree //
            fOutputTree->Fill();
        }
    }
}

bool Finder::FastCuts_ChannelA(const Seeder::Seed& seed_v0a, const Seeder::Seed& seed_v0b, TH1D* cut_flow_hist) const {

    cut_flow_hist->Fill(0.);

    if (Math::SquaredDistance(seed_v0a.pca.xyz, seed_v0b.pca.xyz) > Cuts::ChannelA::Max_DCAbtwV0s * Cuts::ChannelA::Max_DCAbtwV0s) return false;
    cut_flow_hist->Fill(1.);

    return true;
}

bool Finder::SlowCuts(const KF::ChannelA& sexa, TH1D* cut_flow_hist) const {

    if (sexa.SquaredRadius2D() < Cuts::ChannelA::Min_Radius2D * Cuts::ChannelA::Min_Radius2D) return false;
    cut_flow_hist->Fill(2.);

    if (sexa.SquaredDCA_V0A_SV() > Cuts::ChannelA::Max_DCALaSV * Cuts::ChannelA::Max_DCALaSV) return false;
    cut_flow_hist->Fill(3.);

    if (sexa.SquaredDCA_V0B_SV() > Cuts::ChannelA::Max_DCAK0SV * Cuts::ChannelA::Max_DCAK0SV) return false;
    cut_flow_hist->Fill(4.);

    // if (sexa.CPA_wrt(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) < Cuts::ChannelA::Min_CPAwrtPV) return false;
    // cut_flow_hist->Fill(5.);

    // if (sexa.DCA_wrt(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) > Cuts::ChannelA::Max_DCAwrtPV) return false;
    // cut_flow_hist->Fill(6.);

    if (sexa.CPA_V0A_SV() < Cuts::ChannelA::Min_La_CPAwrtSV) return false;
    cut_flow_hist->Fill(7.);

    if (sexa.CPA_V0B_SV() < Cuts::ChannelA::Min_K0S_CPAwrtSV) return false;
    cut_flow_hist->Fill(8.);

    return true;
}

void Finder::Store(const KF::ChannelA& sexa, bool anti_channel) {
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
    fOutput_ChannelA.Chi2NDF = static_cast<float>(sexa.Chi2NDF());
    // -- extra info
    fOutput_ChannelA.E_MinusNucleon = static_cast<float>(sexa.E_MinusNucleon());
    fOutput_ChannelA.AntiChannel = anti_channel;

    // V0A (`Flat::V0`)
    // -- `Flat::State`
    fOutput_ChannelA.V0A.X = sexa.V0A.X();
    fOutput_ChannelA.V0A.Y = sexa.V0A.Y();
    fOutput_ChannelA.V0A.Z = sexa.V0A.Z();
    fOutput_ChannelA.V0A.Px = sexa.V0A.Px();
    fOutput_ChannelA.V0A.Py = sexa.V0A.Py();
    fOutput_ChannelA.V0A.Pz = sexa.V0A.Pz();
    fOutput_ChannelA.V0A.Energy = sexa.V0A.Energy();
    // -- FitVertex info
    fOutput_ChannelA.V0A.Chi2NDF = sexa.V0A.Chi2NDF();
    // -- @ PCA (`Flat::State_NoE`)
    fOutput_ChannelA.V0A_atPCA.X = static_cast<float>(sexa.V0A_at_PCA.X());
    fOutput_ChannelA.V0A_atPCA.Y = static_cast<float>(sexa.V0A_at_PCA.Y());
    fOutput_ChannelA.V0A_atPCA.Z = static_cast<float>(sexa.V0A_at_PCA.Z());
    fOutput_ChannelA.V0A_atPCA.Px = static_cast<float>(sexa.V0A_at_PCA.Px());
    fOutput_ChannelA.V0A_atPCA.Py = static_cast<float>(sexa.V0A_at_PCA.Py());
    fOutput_ChannelA.V0A_atPCA.Pz = static_cast<float>(sexa.V0A_at_PCA.Pz());

    // V0A Neg
    // -- `Flat::Track`
    fOutput_ChannelA.V0A.Neg.X = sexa.V0A.Neg.X();
    fOutput_ChannelA.V0A.Neg.Y = sexa.V0A.Neg.Y();
    fOutput_ChannelA.V0A.Neg.Z = sexa.V0A.Neg.Z();
    fOutput_ChannelA.V0A.Neg.Px = sexa.V0A.Neg.Px();
    fOutput_ChannelA.V0A.Neg.Py = sexa.V0A.Neg.Py();
    fOutput_ChannelA.V0A.Neg.Pz = sexa.V0A.Neg.Pz();
    fOutput_ChannelA.V0A.Neg.Index = sexa.V0A.Neg.Index();
    fOutput_ChannelA.V0A.Neg.DCAxy = sexa.V0A.Neg.DCAxy();
    fOutput_ChannelA.V0A.Neg.DCAz = sexa.V0A.Neg.DCAz();
    fOutput_ChannelA.V0A.Neg.TPCSignal = sexa.V0A.Neg.TPCSignal();
    fOutput_ChannelA.V0A.Neg.NSigmaPion = sexa.V0A.Neg.NSigmaPion();
    fOutput_ChannelA.V0A.Neg.NSigmaKaon = sexa.V0A.Neg.NSigmaKaon();
    fOutput_ChannelA.V0A.Neg.NSigmaProton = sexa.V0A.Neg.NSigmaProton();
    // -- @ V0 `Flat::State_NoE`
    fOutput_ChannelA.V0A.Neg_atV0.X = sexa.V0A.Neg_atPCA_X();
    fOutput_ChannelA.V0A.Neg_atV0.Y = sexa.V0A.Neg_atPCA_Y();
    fOutput_ChannelA.V0A.Neg_atV0.Z = sexa.V0A.Neg_atPCA_Z();
    fOutput_ChannelA.V0A.Neg_atV0.Px = sexa.V0A.Neg_atPCA_Px();
    fOutput_ChannelA.V0A.Neg_atV0.Py = sexa.V0A.Neg_atPCA_Py();
    fOutput_ChannelA.V0A.Neg_atV0.Pz = sexa.V0A.Neg_atPCA_Pz();

    // V0A Pos
    // -- `Flat::Track`
    fOutput_ChannelA.V0A.Pos.X = sexa.V0A.Pos.X();
    fOutput_ChannelA.V0A.Pos.Y = sexa.V0A.Pos.Y();
    fOutput_ChannelA.V0A.Pos.Z = sexa.V0A.Pos.Z();
    fOutput_ChannelA.V0A.Pos.Px = sexa.V0A.Pos.Px();
    fOutput_ChannelA.V0A.Pos.Py = sexa.V0A.Pos.Py();
    fOutput_ChannelA.V0A.Pos.Pz = sexa.V0A.Pos.Pz();
    fOutput_ChannelA.V0A.Pos.Index = sexa.V0A.Pos.Index();
    fOutput_ChannelA.V0A.Pos.DCAxy = sexa.V0A.Pos.DCAxy();
    fOutput_ChannelA.V0A.Pos.DCAz = sexa.V0A.Pos.DCAz();
    fOutput_ChannelA.V0A.Pos.TPCSignal = sexa.V0A.Pos.TPCSignal();
    fOutput_ChannelA.V0A.Pos.NSigmaPion = sexa.V0A.Pos.NSigmaPion();
    fOutput_ChannelA.V0A.Pos.NSigmaKaon = sexa.V0A.Pos.NSigmaKaon();
    fOutput_ChannelA.V0A.Pos.NSigmaProton = sexa.V0A.Pos.NSigmaProton();
    // -- @ V0 (`Flat::State_NoE`)
    fOutput_ChannelA.V0A.Pos_atV0.X = sexa.V0A.Pos_atPCA_X();
    fOutput_ChannelA.V0A.Pos_atV0.Y = sexa.V0A.Pos_atPCA_Y();
    fOutput_ChannelA.V0A.Pos_atV0.Z = sexa.V0A.Pos_atPCA_Z();
    fOutput_ChannelA.V0A.Pos_atV0.Px = sexa.V0A.Pos_atPCA_Px();
    fOutput_ChannelA.V0A.Pos_atV0.Py = sexa.V0A.Pos_atPCA_Py();
    fOutput_ChannelA.V0A.Pos_atV0.Pz = sexa.V0A.Pos_atPCA_Pz();

    // V0B
    // -- `Flat::State`
    fOutput_ChannelA.V0B.X = sexa.V0B.X();
    fOutput_ChannelA.V0B.Y = sexa.V0B.Y();
    fOutput_ChannelA.V0B.Z = sexa.V0B.Z();
    fOutput_ChannelA.V0B.Px = sexa.V0B.Px();
    fOutput_ChannelA.V0B.Py = sexa.V0B.Py();
    fOutput_ChannelA.V0B.Pz = sexa.V0B.Pz();
    fOutput_ChannelA.V0B.Energy = sexa.V0B.Energy();
    // -- FitVertex info
    fOutput_ChannelA.V0B.Chi2NDF = sexa.V0B.Chi2NDF();
    // -- @ PCA (`Flat::State_NoE`)
    fOutput_ChannelA.V0B_atPCA.X = static_cast<float>(sexa.V0B_at_PCA.X());
    fOutput_ChannelA.V0B_atPCA.Y = static_cast<float>(sexa.V0B_at_PCA.Y());
    fOutput_ChannelA.V0B_atPCA.Z = static_cast<float>(sexa.V0B_at_PCA.Z());
    fOutput_ChannelA.V0B_atPCA.Px = static_cast<float>(sexa.V0B_at_PCA.Px());
    fOutput_ChannelA.V0B_atPCA.Py = static_cast<float>(sexa.V0B_at_PCA.Py());
    fOutput_ChannelA.V0B_atPCA.Pz = static_cast<float>(sexa.V0B_at_PCA.Pz());

    // V0B Neg
    // -- `Flat::Track`
    fOutput_ChannelA.V0B.Neg.X = sexa.V0B.Neg.X();
    fOutput_ChannelA.V0B.Neg.Y = sexa.V0B.Neg.Y();
    fOutput_ChannelA.V0B.Neg.Z = sexa.V0B.Neg.Z();
    fOutput_ChannelA.V0B.Neg.Px = sexa.V0B.Neg.Px();
    fOutput_ChannelA.V0B.Neg.Py = sexa.V0B.Neg.Py();
    fOutput_ChannelA.V0B.Neg.Pz = sexa.V0B.Neg.Pz();
    fOutput_ChannelA.V0B.Neg.Index = sexa.V0B.Neg.Index();
    fOutput_ChannelA.V0B.Neg.DCAxy = sexa.V0B.Neg.DCAxy();
    fOutput_ChannelA.V0B.Neg.DCAz = sexa.V0B.Neg.DCAz();
    fOutput_ChannelA.V0B.Neg.TPCSignal = sexa.V0B.Neg.TPCSignal();
    fOutput_ChannelA.V0B.Neg.NSigmaPion = sexa.V0B.Neg.NSigmaPion();
    fOutput_ChannelA.V0B.Neg.NSigmaKaon = sexa.V0B.Neg.NSigmaKaon();
    fOutput_ChannelA.V0B.Neg.NSigmaProton = sexa.V0B.Neg.NSigmaProton();
    // -- @ V0 (`Flat::State_NoE`)
    fOutput_ChannelA.V0B.Neg_atV0.X = sexa.V0B.Neg_atPCA_X();
    fOutput_ChannelA.V0B.Neg_atV0.Y = sexa.V0B.Neg_atPCA_Y();
    fOutput_ChannelA.V0B.Neg_atV0.Z = sexa.V0B.Neg_atPCA_Z();
    fOutput_ChannelA.V0B.Neg_atV0.Px = sexa.V0B.Neg_atPCA_Px();
    fOutput_ChannelA.V0B.Neg_atV0.Py = sexa.V0B.Neg_atPCA_Py();
    fOutput_ChannelA.V0B.Neg_atV0.Pz = sexa.V0B.Neg_atPCA_Pz();

    // V0B Pos
    // -- `Flat::Track`
    fOutput_ChannelA.V0B.Pos.X = sexa.V0B.Pos.X();
    fOutput_ChannelA.V0B.Pos.Y = sexa.V0B.Pos.Y();
    fOutput_ChannelA.V0B.Pos.Z = sexa.V0B.Pos.Z();
    fOutput_ChannelA.V0B.Pos.Px = sexa.V0B.Pos.Px();
    fOutput_ChannelA.V0B.Pos.Py = sexa.V0B.Pos.Py();
    fOutput_ChannelA.V0B.Pos.Pz = sexa.V0B.Pos.Pz();
    fOutput_ChannelA.V0B.Pos.Index = sexa.V0B.Pos.Index();
    fOutput_ChannelA.V0B.Pos.DCAxy = sexa.V0B.Pos.DCAxy();
    fOutput_ChannelA.V0B.Pos.DCAz = sexa.V0B.Pos.DCAz();
    fOutput_ChannelA.V0B.Pos.TPCSignal = sexa.V0B.Pos.TPCSignal();
    fOutput_ChannelA.V0B.Pos.NSigmaPion = sexa.V0B.Pos.NSigmaPion();
    fOutput_ChannelA.V0B.Pos.NSigmaKaon = sexa.V0B.Pos.NSigmaKaon();
    fOutput_ChannelA.V0B.Pos.NSigmaProton = sexa.V0B.Pos.NSigmaProton();
    // -- @ V0 (`Flat::State_NoE`)
    fOutput_ChannelA.V0B.Pos_atV0.X = sexa.V0B.Pos_atPCA_X();
    fOutput_ChannelA.V0B.Pos_atV0.Y = sexa.V0B.Pos_atPCA_Y();
    fOutput_ChannelA.V0B.Pos_atV0.Z = sexa.V0B.Pos_atPCA_Z();
    fOutput_ChannelA.V0B.Pos_atV0.Px = sexa.V0B.Pos_atPCA_Px();
    fOutput_ChannelA.V0B.Pos_atV0.Py = sexa.V0B.Pos_atPCA_Py();
    fOutput_ChannelA.V0B.Pos_atV0.Pz = sexa.V0B.Pos_atPCA_Pz();
}

void Finder::StoreMC(const View::MC::ChannelA& sexa) {
    // `Flat::MC_ChannelA` //

    // `Flat::MC_Sexaquark`
    // -- Before (`Flat::LV`)
    fOutput_MC_ChannelA.Before.Px = sexa.Px();
    fOutput_MC_ChannelA.Before.Py = sexa.Py();
    fOutput_MC_ChannelA.Before.Pz = sexa.Pz();
    fOutput_MC_ChannelA.Before.Energy = static_cast<float>(Truth::Sexaquark::Energy(sexa, fSettings.SexaquarkMass));
    // -- After (`Flat::LV`)
    fOutput_MC_ChannelA.After.Px = static_cast<float>(Truth::Sexaquark::AfterPx(sexa));
    fOutput_MC_ChannelA.After.Py = static_cast<float>(Truth::Sexaquark::AfterPy(sexa));
    fOutput_MC_ChannelA.After.Pz = static_cast<float>(Truth::Sexaquark::AfterPz(sexa));
    fOutput_MC_ChannelA.After.Energy = static_cast<float>(Truth::Sexaquark::AfterE(sexa));
    // -- Nucleon (`Flat::LV`)
    fOutput_MC_ChannelA.Nucleon.Px = sexa.Nucleon_Px();
    fOutput_MC_ChannelA.Nucleon.Py = sexa.Nucleon_Py();
    fOutput_MC_ChannelA.Nucleon.Pz = sexa.Nucleon_Pz();
    fOutput_MC_ChannelA.Nucleon.Energy =
        static_cast<float>(Truth::Sexaquark::NucleonEnergy(sexa, Const::Particle_Mass[Const::ReactionChannel_NucleonPID[fSettings.ReactionChannel]]));
    // -- PV (`Flat::Coordinates`)
    fOutput_MC_ChannelA.PV = fInput_MC_PV;
    // -- SV (`Flat::Coordinates`)
    fOutput_MC_ChannelA.SV.X = sexa.SV_X();
    fOutput_MC_ChannelA.SV.Y = sexa.SV_Y();
    fOutput_MC_ChannelA.SV.Z = sexa.SV_Z();
    // -- reaction id + flags
    fOutput_MC_ChannelA.ReactionID = Truth::Sexaquark::ReactionID(sexa);
    fOutput_MC_ChannelA.IsSignal = Truth::Sexaquark::IsSignal(sexa);
    fOutput_MC_ChannelA.IsHybrid = Truth::Sexaquark::IsHybrid(sexa);

    // V0A (`Flat::MC_V0`)
    // -- `Flat::MC`
    fOutput_MC_ChannelA.V0A.McEntry = static_cast<int>(sexa.V0A.Entry);
    fOutput_MC_ChannelA.V0A.PdgCode = sexa.V0A.PdgCode();
    fOutput_MC_ChannelA.V0A.ReactionID = sexa.V0A.ReactionID();
    fOutput_MC_ChannelA.V0A.IsTrue = sexa.V0A.IsTrue();
    fOutput_MC_ChannelA.V0A.IsSignal = sexa.V0A.IsSignal();
    fOutput_MC_ChannelA.V0A.IsSecondary = sexa.V0A.IsSecondary();
    // -- `Flat::LV`
    fOutput_MC_ChannelA.V0A.Px = sexa.V0A.Px();
    fOutput_MC_ChannelA.V0A.Py = sexa.V0A.Py();
    fOutput_MC_ChannelA.V0A.Pz = sexa.V0A.Pz();
    fOutput_MC_ChannelA.V0A.Energy = sexa.V0A.Energy();
    // -- mother info (`Flat::MC_Id`)
    fOutput_MC_ChannelA.V0A.Mother.McEntry = sexa.V0A.Mother_Entry();
    fOutput_MC_ChannelA.V0A.Mother.PdgCode = sexa.V0A.Mother_PdgCode();
    // -- @ decay (`Flat::Coordinates`)
    fOutput_MC_ChannelA.V0A.AtDecay.X = sexa.V0A.Decay_X();
    fOutput_MC_ChannelA.V0A.AtDecay.Y = sexa.V0A.Decay_Y();
    fOutput_MC_ChannelA.V0A.AtDecay.Z = sexa.V0A.Decay_Z();
    // -- hybrid flag
    fOutput_MC_ChannelA.V0A.IsHybrid = sexa.V0A.IsHybrid();

    // V0A neg
    // -- `Flat::MC`
    fOutput_MC_ChannelA.V0A.Neg.McEntry = sexa.V0A.Neg_Entry();
    fOutput_MC_ChannelA.V0A.Neg.PdgCode = sexa.V0A.Neg_PdgCode();
    fOutput_MC_ChannelA.V0A.Neg.ReactionID = sexa.V0A.Neg_ReactionID();
    fOutput_MC_ChannelA.V0A.Neg.IsTrue = sexa.V0A.Neg_IsTrue();
    fOutput_MC_ChannelA.V0A.Neg.IsSignal = sexa.V0A.Neg_IsSignal();
    fOutput_MC_ChannelA.V0A.Neg.IsSecondary = sexa.V0A.Neg_IsSecondary();
    // -- `Flat::PxPyPz`
    fOutput_MC_ChannelA.V0A.Neg_Momentum.Px = sexa.V0A.Neg_Px();
    fOutput_MC_ChannelA.V0A.Neg_Momentum.Py = sexa.V0A.Neg_Py();
    fOutput_MC_ChannelA.V0A.Neg_Momentum.Pz = sexa.V0A.Neg_Pz();

    // V0A pos
    // -- `Flat::MC`
    fOutput_MC_ChannelA.V0A.Pos.McEntry = sexa.V0A.Pos_Entry();
    fOutput_MC_ChannelA.V0A.Pos.PdgCode = sexa.V0A.Pos_PdgCode();
    fOutput_MC_ChannelA.V0A.Pos.ReactionID = sexa.V0A.Pos_ReactionID();
    fOutput_MC_ChannelA.V0A.Pos.IsTrue = sexa.V0A.Pos_IsTrue();
    fOutput_MC_ChannelA.V0A.Pos.IsSignal = sexa.V0A.Pos_IsSignal();
    fOutput_MC_ChannelA.V0A.Pos.IsSecondary = sexa.V0A.Pos_IsSecondary();
    // -- `Flat::PxPyPz`
    fOutput_MC_ChannelA.V0A.Pos_Momentum.Px = sexa.V0A.Pos_Px();
    fOutput_MC_ChannelA.V0A.Pos_Momentum.Py = sexa.V0A.Pos_Py();
    fOutput_MC_ChannelA.V0A.Pos_Momentum.Pz = sexa.V0A.Pos_Pz();

    // V0B (`Flat::MC_V0`)
    // -- `Flat::MC`
    fOutput_MC_ChannelA.V0B.McEntry = static_cast<int>(sexa.V0B.Entry);
    fOutput_MC_ChannelA.V0B.PdgCode = sexa.V0B.PdgCode();
    fOutput_MC_ChannelA.V0B.ReactionID = sexa.V0B.ReactionID();
    fOutput_MC_ChannelA.V0B.IsTrue = sexa.V0B.IsTrue();
    fOutput_MC_ChannelA.V0B.IsSignal = sexa.V0B.IsSignal();
    fOutput_MC_ChannelA.V0B.IsSecondary = sexa.V0B.IsSecondary();
    // -- `Flat::LV`
    fOutput_MC_ChannelA.V0B.Px = sexa.V0B.Px();
    fOutput_MC_ChannelA.V0B.Py = sexa.V0B.Py();
    fOutput_MC_ChannelA.V0B.Pz = sexa.V0B.Pz();
    fOutput_MC_ChannelA.V0B.Energy = sexa.V0B.Energy();
    // -- mother info (`Flat::MC_Id`)
    fOutput_MC_ChannelA.V0B.Mother.McEntry = sexa.V0B.Mother_Entry();
    fOutput_MC_ChannelA.V0B.Mother.PdgCode = sexa.V0B.Mother_PdgCode();
    // -- @ decay (`Flat::Coordinates`)
    fOutput_MC_ChannelA.V0B.AtDecay.X = sexa.V0B.Decay_X();
    fOutput_MC_ChannelA.V0B.AtDecay.Y = sexa.V0B.Decay_Y();
    fOutput_MC_ChannelA.V0B.AtDecay.Z = sexa.V0B.Decay_Z();
    // -- hybrid flag
    fOutput_MC_ChannelA.V0B.IsHybrid = sexa.V0B.IsHybrid();

    // V0B neg
    // -- `Flat::MC`
    fOutput_MC_ChannelA.V0B.Neg.McEntry = sexa.V0B.Neg_Entry();
    fOutput_MC_ChannelA.V0B.Neg.PdgCode = sexa.V0B.Neg_PdgCode();
    fOutput_MC_ChannelA.V0B.Neg.ReactionID = sexa.V0B.Neg_ReactionID();
    fOutput_MC_ChannelA.V0B.Neg.IsTrue = sexa.V0B.Neg_IsTrue();
    fOutput_MC_ChannelA.V0B.Neg.IsSignal = sexa.V0B.Neg_IsSignal();
    fOutput_MC_ChannelA.V0B.Neg.IsSecondary = sexa.V0B.Neg_IsSecondary();
    // -- `Flat::PxPyPz`
    fOutput_MC_ChannelA.V0B.Neg_Momentum.Px = sexa.V0B.Neg_Px();
    fOutput_MC_ChannelA.V0B.Neg_Momentum.Py = sexa.V0B.Neg_Py();
    fOutput_MC_ChannelA.V0B.Neg_Momentum.Pz = sexa.V0B.Neg_Pz();

    // V0B pos
    // -- `Flat::MC`
    fOutput_MC_ChannelA.V0B.Pos.McEntry = sexa.V0B.Pos_Entry();
    fOutput_MC_ChannelA.V0B.Pos.PdgCode = sexa.V0B.Pos_PdgCode();
    fOutput_MC_ChannelA.V0B.Pos.ReactionID = sexa.V0B.Pos_ReactionID();
    fOutput_MC_ChannelA.V0B.Pos.IsTrue = sexa.V0B.Pos_IsTrue();
    fOutput_MC_ChannelA.V0B.Pos.IsSignal = sexa.V0B.Pos_IsSignal();
    fOutput_MC_ChannelA.V0B.Pos.IsSecondary = sexa.V0B.Pos_IsSecondary();
    // -- `Flat::PxPyPz`
    fOutput_MC_ChannelA.V0B.Pos_Momentum.Px = sexa.V0B.Pos_Px();
    fOutput_MC_ChannelA.V0B.Pos_Momentum.Py = sexa.V0B.Pos_Py();
    fOutput_MC_ChannelA.V0B.Pos_Momentum.Pz = sexa.V0B.Pos_Pz();
}

void Finder::StoreDummyMC_ChannelA() {
    // `Flat::MC_ChannelA` //

    // `Flat::MC_Sexaquark`
    // -- Before (`Flat::LV`)
    fOutput_MC_ChannelA.Before.Px = Const::DummyFloat;
    fOutput_MC_ChannelA.Before.Py = Const::DummyFloat;
    fOutput_MC_ChannelA.Before.Pz = Const::DummyFloat;
    fOutput_MC_ChannelA.Before.Energy = Const::DummyFloat;
    // -- After (`Flat::LV`)
    fOutput_MC_ChannelA.After.Px = Const::DummyFloat;
    fOutput_MC_ChannelA.After.Py = Const::DummyFloat;
    fOutput_MC_ChannelA.After.Pz = Const::DummyFloat;
    fOutput_MC_ChannelA.After.Energy = Const::DummyFloat;
    // -- Nucleon (`Flat::LV`)
    fOutput_MC_ChannelA.Nucleon.Px = Const::DummyFloat;
    fOutput_MC_ChannelA.Nucleon.Py = Const::DummyFloat;
    fOutput_MC_ChannelA.Nucleon.Pz = Const::DummyFloat;
    fOutput_MC_ChannelA.Nucleon.Energy = Const::DummyFloat;
    // -- PV (`Flat::Coordinates`)
    fOutput_MC_ChannelA.PV.X = Const::DummyFloat;
    fOutput_MC_ChannelA.PV.Y = Const::DummyFloat;
    fOutput_MC_ChannelA.PV.Z = Const::DummyFloat;
    // -- SV (`Flat::Coordinates`)
    fOutput_MC_ChannelA.SV.X = Const::DummyFloat;
    fOutput_MC_ChannelA.SV.Y = Const::DummyFloat;
    fOutput_MC_ChannelA.SV.Z = Const::DummyFloat;
    // -- reaction id + flags
    fOutput_MC_ChannelA.ReactionID = Const::DummyInt;
    fOutput_MC_ChannelA.IsSignal = false;
    fOutput_MC_ChannelA.IsHybrid = false;

    // V0A (`Flat::MC_V0`)
    // -- `Flat::MC`
    fOutput_MC_ChannelA.V0A.McEntry = Const::DummyInt;
    fOutput_MC_ChannelA.V0A.PdgCode = Const::DummyInt;
    fOutput_MC_ChannelA.V0A.ReactionID = Const::DummyInt;
    fOutput_MC_ChannelA.V0A.IsTrue = false;
    fOutput_MC_ChannelA.V0A.IsSignal = false;
    fOutput_MC_ChannelA.V0A.IsSecondary = false;
    // -- `Flat::LV`
    fOutput_MC_ChannelA.V0A.Px = Const::DummyFloat;
    fOutput_MC_ChannelA.V0A.Py = Const::DummyFloat;
    fOutput_MC_ChannelA.V0A.Pz = Const::DummyFloat;
    fOutput_MC_ChannelA.V0A.Energy = Const::DummyFloat;
    // -- mother info (`Flat::MC_Id`)
    fOutput_MC_ChannelA.V0A.Mother.McEntry = Const::DummyInt;
    fOutput_MC_ChannelA.V0A.Mother.PdgCode = Const::DummyInt;
    // -- @ decay (`Flat::Coordinates`)
    fOutput_MC_ChannelA.V0A.AtDecay.X = Const::DummyFloat;
    fOutput_MC_ChannelA.V0A.AtDecay.Y = Const::DummyFloat;
    fOutput_MC_ChannelA.V0A.AtDecay.Z = Const::DummyFloat;
    // -- hybrid flag
    fOutput_MC_ChannelA.V0A.IsHybrid = false;

    // V0A neg
    // -- `Flat::MC`
    fOutput_MC_ChannelA.V0A.Neg.McEntry = Const::DummyInt;
    fOutput_MC_ChannelA.V0A.Neg.PdgCode = Const::DummyInt;
    fOutput_MC_ChannelA.V0A.Neg.ReactionID = Const::DummyInt;
    fOutput_MC_ChannelA.V0A.Neg.IsTrue = false;
    fOutput_MC_ChannelA.V0A.Neg.IsSignal = false;
    fOutput_MC_ChannelA.V0A.Neg.IsSecondary = false;
    // -- `Flat::PxPyPz`
    fOutput_MC_ChannelA.V0A.Neg_Momentum.Px = Const::DummyFloat;
    fOutput_MC_ChannelA.V0A.Neg_Momentum.Py = Const::DummyFloat;
    fOutput_MC_ChannelA.V0A.Neg_Momentum.Pz = Const::DummyFloat;

    // V0A pos
    // -- `Flat::MC`
    fOutput_MC_ChannelA.V0A.Pos.McEntry = Const::DummyInt;
    fOutput_MC_ChannelA.V0A.Pos.PdgCode = Const::DummyInt;
    fOutput_MC_ChannelA.V0A.Pos.ReactionID = Const::DummyInt;
    fOutput_MC_ChannelA.V0A.Pos.IsTrue = false;
    fOutput_MC_ChannelA.V0A.Pos.IsSignal = false;
    fOutput_MC_ChannelA.V0A.Pos.IsSecondary = false;
    // -- `Flat::PxPyPz`
    fOutput_MC_ChannelA.V0A.Pos_Momentum.Px = Const::DummyFloat;
    fOutput_MC_ChannelA.V0A.Pos_Momentum.Py = Const::DummyFloat;
    fOutput_MC_ChannelA.V0A.Pos_Momentum.Pz = Const::DummyFloat;

    // V0B (`Flat::MC_V0`)
    // -- `Flat::MC`
    fOutput_MC_ChannelA.V0B.McEntry = Const::DummyInt;
    fOutput_MC_ChannelA.V0B.PdgCode = Const::DummyInt;
    fOutput_MC_ChannelA.V0B.ReactionID = Const::DummyInt;
    fOutput_MC_ChannelA.V0B.IsTrue = false;
    fOutput_MC_ChannelA.V0B.IsSignal = false;
    fOutput_MC_ChannelA.V0B.IsSecondary = false;
    // -- `Flat::LV`
    fOutput_MC_ChannelA.V0B.Px = Const::DummyFloat;
    fOutput_MC_ChannelA.V0B.Py = Const::DummyFloat;
    fOutput_MC_ChannelA.V0B.Pz = Const::DummyFloat;
    fOutput_MC_ChannelA.V0B.Energy = Const::DummyFloat;
    // -- mother info (`Flat::MC_Id`)
    fOutput_MC_ChannelA.V0B.Mother.McEntry = Const::DummyInt;
    fOutput_MC_ChannelA.V0B.Mother.PdgCode = Const::DummyInt;
    // -- @ decay (`Flat::Coordinates`)
    fOutput_MC_ChannelA.V0B.AtDecay.X = Const::DummyFloat;
    fOutput_MC_ChannelA.V0B.AtDecay.Y = Const::DummyFloat;
    fOutput_MC_ChannelA.V0B.AtDecay.Z = Const::DummyFloat;
    // -- hybrid flag
    fOutput_MC_ChannelA.V0B.IsHybrid = false;

    // V0B neg
    // -- `Flat::MC`
    fOutput_MC_ChannelA.V0B.Neg.McEntry = Const::DummyInt;
    fOutput_MC_ChannelA.V0B.Neg.PdgCode = Const::DummyInt;
    fOutput_MC_ChannelA.V0B.Neg.ReactionID = Const::DummyInt;
    fOutput_MC_ChannelA.V0B.Neg.IsTrue = false;
    fOutput_MC_ChannelA.V0B.Neg.IsSignal = false;
    fOutput_MC_ChannelA.V0B.Neg.IsSecondary = false;
    // -- `Flat::PxPyPz`
    fOutput_MC_ChannelA.V0B.Neg_Momentum.Px = Const::DummyFloat;
    fOutput_MC_ChannelA.V0B.Neg_Momentum.Py = Const::DummyFloat;
    fOutput_MC_ChannelA.V0B.Neg_Momentum.Pz = Const::DummyFloat;

    // V0B pos
    // -- `Flat::MC`
    fOutput_MC_ChannelA.V0B.Pos.McEntry = Const::DummyInt;
    fOutput_MC_ChannelA.V0B.Pos.PdgCode = Const::DummyInt;
    fOutput_MC_ChannelA.V0B.Pos.ReactionID = Const::DummyInt;
    fOutput_MC_ChannelA.V0B.Pos.IsTrue = false;
    fOutput_MC_ChannelA.V0B.Pos.IsSignal = false;
    fOutput_MC_ChannelA.V0B.Pos.IsSecondary = false;
    // -- `Flat::PxPyPz`
    fOutput_MC_ChannelA.V0B.Pos_Momentum.Px = Const::DummyFloat;
    fOutput_MC_ChannelA.V0B.Pos_Momentum.Py = Const::DummyFloat;
    fOutput_MC_ChannelA.V0B.Pos_Momentum.Pz = Const::DummyFloat;
}

// ## Channel D ZONE ## //

void Finder::FindSexaquarks_ChannelD(bool anti_channel) {

    // determine properties based on anti-channel or not
    // -- kaon
    const Storage::Vector::Tracks* Packed_Kaons{&fInput_PosKaons};
    const Storage::Vector::MC_Tracks* MC_Kaons{&fInput_MC_PosKaons};
    const double ka_mass = Const::Particle_Mass[PID_StableParticle::NegKaon];  // NOTE: same as mass of pos. kaon
    if (anti_channel) {
        Packed_Kaons = &fInput_NegKaons;
        MC_Kaons = &fInput_MC_NegKaons;
    }
    // -- v0
    const Storage::Vector::V0s* Packed_V0s{&fInput_AntiLambdas};
    const Storage::Vector::MC_V0s* MC_V0s{&fInput_MC_AntiLambdas};
    if (anti_channel) {
        Packed_V0s = &fInput_Lambdas;
        MC_V0s = &fInput_MC_Lambdas;
    }
    // -- cut flow hist
    TH1D* hist{anti_channel ? fHist_CutFlow_AntiChannel.get() : fHist_CutFlow.get()};

    // loop over all possible pairs of (anti)lambda + (pos/neg)kaon //
    size_t n_kaons{Packed_Kaons->X->size()};
    size_t n_v0{Packed_V0s->X->size()};
    for (size_t ka_entry = 0; ka_entry < n_kaons; ++ka_entry) {
        for (size_t v0_entry = 0; v0_entry < n_v0; ++v0_entry) {

            // get views //
            View::Rec::Track ka{Packed_Kaons, ka_entry};
            View::Rec::V0 v0{Packed_V0s, v0_entry};

            // sanity check //
            if (v0.Neg.Index() == ka.Index() || v0.Pos.Index() == ka.Index()) continue;

            // PCAs (1) //
            auto [seed_ka, seed_v0, pca_cache] = Seeder::HelixLine::FastCorrectPCAs(ka, v0, fInput_Event.MagneticField);

            // apply cuts (1) //
            if (!FastCuts_ChannelD(seed_ka, seed_v0, hist)) continue;

            // PCAs derivatives //
            auto [deriv_v0, deriv_ka] = Seeder::HelixLine::ComputeDerivatives(seed_ka, seed_v0, pca_cache);

            // fit vertex //
            auto fit = KF::FitVertex(ka, v0, ka_mass,                           //
                                     {seed_ka, deriv_ka}, {seed_v0, deriv_v0},  //
                                     fInput_Event.MagneticField);
            KF::ChannelD sexa(fit, seed_ka.pca, seed_v0.pca, ka, v0);

            // apply cuts (2) //
            if (!SlowCuts(sexa, hist)) continue;

            // store //
            Store(sexa, anti_channel);

            if (IsMC()) {
                View::MC::PackedTrack mc_ka(MC_Kaons, ka_entry);
                View::MC::PackedV0 mc_v0(MC_V0s, v0_entry);
                View::MC::ChannelD mc_sexa(&fInput_Injected, mc_ka, mc_v0);
                if (View::IsValid(mc_sexa)) {
                    StoreMC(mc_sexa);
                } else {
                    StoreDummyMC_ChannelD();
                }
            }

            // fill flat tree //
            fOutputTree->Fill();
        }
    }
}

bool Finder::FastCuts_ChannelD(const Seeder::Seed& seed_ka, const Seeder::Seed& seed_v0, TH1D* cut_flow_hist) const {

    cut_flow_hist->Fill(0.);

    if (Math::SquaredDistance(seed_ka.pca.xyz, seed_v0.pca.xyz) > Cuts::ChannelD::Max_DCAKaLa * Cuts::ChannelD::Max_DCAKaLa) return false;
    cut_flow_hist->Fill(1.);

    return true;
}

bool Finder::SlowCuts(const KF::ChannelD& sexa, TH1D* cut_flow_hist) const {

    double sq_radius_2d = sexa.SquaredRadius2D();
    if (sq_radius_2d < Cuts::ChannelD::Min_Radius2D * Cuts::ChannelD::Min_Radius2D ||
        sq_radius_2d > Cuts::ChannelD::Max_Radius2D * Cuts::ChannelD::Max_Radius2D) {
        return false;
    }
    cut_flow_hist->Fill(2.);

    // if (sexa.AbsRapidity_MinusNucleon() > Cuts::ChannelD::AbsMax_Rapidity) return false;  // PENDING: kinematics, affected by Fermi motion
    // cut_flow_hist->Fill(3.);

    // if (sexa.CPA_Vertex(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) < Cuts::ChannelD::Min_CPAwrtPV ||
    // sexa.CPA_Vertex(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) > Cuts::ChannelD::Max_CPAwrtPV) {
    // return false;  // PENDING: kinematics, affected by Fermi motion
    // }
    // cut_flow_hist->Fill(3.);

    if (sexa.SquaredDCA_V0_SV() > Cuts::ChannelD::Max_DCALaSV * Cuts::ChannelD::Max_DCALaSV) return false;
    cut_flow_hist->Fill(4.);

    if (sexa.SquaredDCA_Kaon_SV() > Cuts::ChannelD::Max_DCAKaSV * Cuts::ChannelD::Max_DCAKaSV) return false;
    cut_flow_hist->Fill(5.);

    // if (sexa.DCA_V0Neg_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelD::Max_DCALaNegSV) return false;
    // cut_flow_hist->Fill(6.);

    // if (sexa.DCA_V0Pos_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelD::Max_DCALaPosSV) return false;
    // cut_flow_hist->Fill(7.);

    return true;
}

void Finder::Store(const KF::ChannelD& sexa, bool anti_channel) {
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
    // -- Fit info
    fOutput_ChannelD.Chi2NDF = static_cast<float>(sexa.Chi2NDF());
    // -- Extra info
    fOutput_ChannelD.E_MinusNucleon = static_cast<float>(sexa.E_MinusNucleon());
    fOutput_ChannelD.AntiChannel = anti_channel;

    // Kaon (`Flat::Track`)
    // -- `Flat::State_NoE`
    fOutput_ChannelD.Kaon.X = sexa.Kaon.X();
    fOutput_ChannelD.Kaon.Y = sexa.Kaon.Y();
    fOutput_ChannelD.Kaon.Z = sexa.Kaon.Z();
    fOutput_ChannelD.Kaon.Px = sexa.Kaon.Px();
    fOutput_ChannelD.Kaon.Py = sexa.Kaon.Py();
    fOutput_ChannelD.Kaon.Pz = sexa.Kaon.Pz();
    // -- Rest of info
    fOutput_ChannelD.Kaon.Index = sexa.Kaon.Index();
    fOutput_ChannelD.Kaon.DCAxy = sexa.Kaon.DCAxy();
    fOutput_ChannelD.Kaon.DCAz = sexa.Kaon.DCAz();
    fOutput_ChannelD.Kaon.TPCSignal = sexa.Kaon.TPCSignal();
    fOutput_ChannelD.Kaon.NSigmaPion = sexa.Kaon.NSigmaPion();
    fOutput_ChannelD.Kaon.NSigmaKaon = sexa.Kaon.NSigmaKaon();
    fOutput_ChannelD.Kaon.NSigmaProton = sexa.Kaon.NSigmaProton();
    // -- @ PCA (`Flat::State_NoE`)
    fOutput_ChannelD.Kaon_atPCA.X = static_cast<float>(sexa.Kaon_at_PCA.X());
    fOutput_ChannelD.Kaon_atPCA.Y = static_cast<float>(sexa.Kaon_at_PCA.Y());
    fOutput_ChannelD.Kaon_atPCA.Z = static_cast<float>(sexa.Kaon_at_PCA.Z());
    fOutput_ChannelD.Kaon_atPCA.Px = static_cast<float>(sexa.Kaon_at_PCA.Px());
    fOutput_ChannelD.Kaon_atPCA.Py = static_cast<float>(sexa.Kaon_at_PCA.Py());
    fOutput_ChannelD.Kaon_atPCA.Pz = static_cast<float>(sexa.Kaon_at_PCA.Pz());

    // V0 (`Flat::V0`)
    // -- `Flat::State`
    fOutput_ChannelD.V0.X = sexa.V0.X();
    fOutput_ChannelD.V0.Y = sexa.V0.Y();
    fOutput_ChannelD.V0.Z = sexa.V0.Z();
    fOutput_ChannelD.V0.Px = sexa.V0.Px();
    fOutput_ChannelD.V0.Py = sexa.V0.Py();
    fOutput_ChannelD.V0.Pz = sexa.V0.Pz();
    fOutput_ChannelD.V0.Energy = sexa.V0.Energy();
    // -- Fit info
    fOutput_ChannelD.V0.Chi2NDF = sexa.V0.Chi2NDF();
    // -- @ PCA (`Flat::State_NoE`)
    fOutput_ChannelD.V0_atPCA.X = static_cast<float>(sexa.V0_at_PCA.X());
    fOutput_ChannelD.V0_atPCA.Y = static_cast<float>(sexa.V0_at_PCA.Y());
    fOutput_ChannelD.V0_atPCA.Z = static_cast<float>(sexa.V0_at_PCA.Z());
    fOutput_ChannelD.V0_atPCA.Px = static_cast<float>(sexa.V0_at_PCA.Px());
    fOutput_ChannelD.V0_atPCA.Py = static_cast<float>(sexa.V0_at_PCA.Py());
    fOutput_ChannelD.V0_atPCA.Pz = static_cast<float>(sexa.V0_at_PCA.Pz());

    // V0 Neg (`Flat::Track`)
    // -- `Flat::State_NoE`
    fOutput_ChannelD.V0.Neg.X = sexa.V0.Neg.X();
    fOutput_ChannelD.V0.Neg.Y = sexa.V0.Neg.Y();
    fOutput_ChannelD.V0.Neg.Z = sexa.V0.Neg.Z();
    fOutput_ChannelD.V0.Neg.Px = sexa.V0.Neg.Px();
    fOutput_ChannelD.V0.Neg.Py = sexa.V0.Neg.Py();
    fOutput_ChannelD.V0.Neg.Pz = sexa.V0.Neg.Pz();
    // -- Rest of info
    fOutput_ChannelD.V0.Neg.Index = sexa.V0.Neg.Index();
    fOutput_ChannelD.V0.Neg.DCAxy = sexa.V0.Neg.DCAxy();
    fOutput_ChannelD.V0.Neg.DCAz = sexa.V0.Neg.DCAz();
    fOutput_ChannelD.V0.Neg.TPCSignal = sexa.V0.Neg.TPCSignal();
    fOutput_ChannelD.V0.Neg.NSigmaPion = sexa.V0.Neg.NSigmaPion();
    fOutput_ChannelD.V0.Neg.NSigmaKaon = sexa.V0.Neg.NSigmaKaon();
    fOutput_ChannelD.V0.Neg.NSigmaProton = sexa.V0.Neg.NSigmaProton();
    // -- V0 Neg @ V0 (`Flat::State_NoE`)
    fOutput_ChannelD.V0.Neg_atV0.X = sexa.V0.Neg_atPCA_X();
    fOutput_ChannelD.V0.Neg_atV0.Y = sexa.V0.Neg_atPCA_Y();
    fOutput_ChannelD.V0.Neg_atV0.Z = sexa.V0.Neg_atPCA_Z();
    fOutput_ChannelD.V0.Neg_atV0.Px = sexa.V0.Neg_atPCA_Px();
    fOutput_ChannelD.V0.Neg_atV0.Py = sexa.V0.Neg_atPCA_Py();
    fOutput_ChannelD.V0.Neg_atV0.Pz = sexa.V0.Neg_atPCA_Pz();

    // V0 Pos (`Flat::Track`)
    // -- `Flat::State_NoE`
    fOutput_ChannelD.V0.Pos.X = sexa.V0.Pos.X();
    fOutput_ChannelD.V0.Pos.Y = sexa.V0.Pos.Y();
    fOutput_ChannelD.V0.Pos.Z = sexa.V0.Pos.Z();
    fOutput_ChannelD.V0.Pos.Px = sexa.V0.Pos.Px();
    fOutput_ChannelD.V0.Pos.Py = sexa.V0.Pos.Py();
    fOutput_ChannelD.V0.Pos.Pz = sexa.V0.Pos.Pz();
    // -- Rest of info
    fOutput_ChannelD.V0.Pos.Index = sexa.V0.Pos.Index();
    fOutput_ChannelD.V0.Pos.DCAxy = sexa.V0.Pos.DCAxy();
    fOutput_ChannelD.V0.Pos.DCAz = sexa.V0.Pos.DCAz();
    fOutput_ChannelD.V0.Pos.TPCSignal = sexa.V0.Pos.TPCSignal();
    fOutput_ChannelD.V0.Pos.NSigmaPion = sexa.V0.Pos.NSigmaPion();
    fOutput_ChannelD.V0.Pos.NSigmaKaon = sexa.V0.Pos.NSigmaKaon();
    fOutput_ChannelD.V0.Pos.NSigmaProton = sexa.V0.Pos.NSigmaProton();
    // -- V0 Pos @ V0 (`Flat::State_NoE`)
    fOutput_ChannelD.V0.Pos_atV0.X = sexa.V0.Pos_atPCA_X();
    fOutput_ChannelD.V0.Pos_atV0.Y = sexa.V0.Pos_atPCA_Y();
    fOutput_ChannelD.V0.Pos_atV0.Z = sexa.V0.Pos_atPCA_Z();
    fOutput_ChannelD.V0.Pos_atV0.Px = sexa.V0.Pos_atPCA_Px();
    fOutput_ChannelD.V0.Pos_atV0.Py = sexa.V0.Pos_atPCA_Py();
    fOutput_ChannelD.V0.Pos_atV0.Pz = sexa.V0.Pos_atPCA_Pz();
}

void Finder::StoreMC(const View::MC::ChannelD& sexa) {
    // `Flat::MC_ChannelD` //

    // `Flat::MC_Sexaquark`
    // -- Before (`Flat::LV`)
    fOutput_MC_ChannelD.Before.Px = sexa.Px();
    fOutput_MC_ChannelD.Before.Py = sexa.Py();
    fOutput_MC_ChannelD.Before.Pz = sexa.Pz();
    fOutput_MC_ChannelD.Before.Energy = static_cast<float>(Truth::Sexaquark::Energy(sexa, fSettings.SexaquarkMass));
    // -- After (`Flat::LV`)
    fOutput_MC_ChannelD.After.Px = static_cast<float>(Truth::Sexaquark::AfterPx(sexa));
    fOutput_MC_ChannelD.After.Py = static_cast<float>(Truth::Sexaquark::AfterPy(sexa));
    fOutput_MC_ChannelD.After.Pz = static_cast<float>(Truth::Sexaquark::AfterPz(sexa));
    fOutput_MC_ChannelD.After.Energy = static_cast<float>(Truth::Sexaquark::AfterE(sexa));
    // -- Nucleon (`Flat::LV`)
    fOutput_MC_ChannelD.Nucleon.Px = sexa.Nucleon_Px();
    fOutput_MC_ChannelD.Nucleon.Py = sexa.Nucleon_Py();
    fOutput_MC_ChannelD.Nucleon.Pz = sexa.Nucleon_Pz();
    fOutput_MC_ChannelD.Nucleon.Energy =
        static_cast<float>(Truth::Sexaquark::NucleonEnergy(sexa, Const::Particle_Mass[Const::ReactionChannel_NucleonPID[fSettings.ReactionChannel]]));
    // -- PV (`Flat::Coordinates`)
    fOutput_MC_ChannelD.PV = fInput_MC_PV;
    // -- SV (`Flat::Coordinates`)
    fOutput_MC_ChannelD.SV.X = sexa.SV_X();
    fOutput_MC_ChannelD.SV.Y = sexa.SV_Y();
    fOutput_MC_ChannelD.SV.Z = sexa.SV_Z();
    // -- Reaction ID + Flags
    fOutput_MC_ChannelD.ReactionID = Truth::Sexaquark::ReactionID(sexa);
    fOutput_MC_ChannelD.IsSignal = Truth::Sexaquark::IsSignal(sexa);
    fOutput_MC_ChannelD.IsHybrid = Truth::Sexaquark::IsHybrid(sexa);

    // V0 (`Flat::MC_V0`)
    // -- `Flat::MC`
    fOutput_MC_ChannelD.V0.McEntry = static_cast<int>(sexa.V0.Entry);
    fOutput_MC_ChannelD.V0.PdgCode = sexa.V0.PdgCode();
    fOutput_MC_ChannelD.V0.ReactionID = sexa.V0.ReactionID();
    fOutput_MC_ChannelD.V0.IsTrue = sexa.V0.IsTrue();
    fOutput_MC_ChannelD.V0.IsSignal = sexa.V0.IsSignal();
    fOutput_MC_ChannelD.V0.IsSecondary = sexa.V0.IsSecondary();
    // -- `Flat::LV`
    fOutput_MC_ChannelD.V0.Px = sexa.V0.Px();
    fOutput_MC_ChannelD.V0.Py = sexa.V0.Py();
    fOutput_MC_ChannelD.V0.Pz = sexa.V0.Pz();
    fOutput_MC_ChannelD.V0.Energy = sexa.V0.Energy();
    // -- Mother (`Flat::MC_Id`)
    fOutput_MC_ChannelD.V0.Mother.McEntry = sexa.V0.Mother_Entry();
    fOutput_MC_ChannelD.V0.Mother.PdgCode = sexa.V0.Mother_PdgCode();
    // -- @ decay (`Flat::Coordinates`)
    fOutput_MC_ChannelD.V0.AtDecay.X = sexa.V0.Decay_X();
    fOutput_MC_ChannelD.V0.AtDecay.Y = sexa.V0.Decay_Y();
    fOutput_MC_ChannelD.V0.AtDecay.Z = sexa.V0.Decay_Z();
    // -- hybrid flag
    fOutput_MC_ChannelD.V0.IsHybrid = sexa.V0.IsHybrid();

    // V0 Neg
    // -- `Flat::MC`
    fOutput_MC_ChannelD.V0.Neg.McEntry = sexa.V0.Neg_Entry();
    fOutput_MC_ChannelD.V0.Neg.PdgCode = sexa.V0.Neg_PdgCode();
    fOutput_MC_ChannelD.V0.Neg.ReactionID = sexa.V0.Neg_ReactionID();
    fOutput_MC_ChannelD.V0.Neg.IsTrue = sexa.V0.Neg_IsTrue();
    fOutput_MC_ChannelD.V0.Neg.IsSignal = sexa.V0.Neg_IsSignal();
    fOutput_MC_ChannelD.V0.Neg.IsSecondary = sexa.V0.Neg_IsSecondary();
    // -- `Flat::PxPyPz`
    fOutput_MC_ChannelD.V0.Neg_Momentum.Px = sexa.V0.Neg_Px();
    fOutput_MC_ChannelD.V0.Neg_Momentum.Py = sexa.V0.Neg_Py();
    fOutput_MC_ChannelD.V0.Neg_Momentum.Pz = sexa.V0.Neg_Pz();

    // V0 Pos
    // -- `Flat::MC`
    fOutput_MC_ChannelD.V0.Pos.McEntry = sexa.V0.Pos_Entry();
    fOutput_MC_ChannelD.V0.Pos.PdgCode = sexa.V0.Pos_PdgCode();
    fOutput_MC_ChannelD.V0.Pos.ReactionID = sexa.V0.Pos_ReactionID();
    fOutput_MC_ChannelD.V0.Pos.IsTrue = sexa.V0.Pos_IsTrue();
    fOutput_MC_ChannelD.V0.Pos.IsSignal = sexa.V0.Pos_IsSignal();
    fOutput_MC_ChannelD.V0.Pos.IsSecondary = sexa.V0.Pos_IsSecondary();
    // -- `Flat::PxPyPz`
    fOutput_MC_ChannelD.V0.Pos_Momentum.Px = sexa.V0.Pos_Px();
    fOutput_MC_ChannelD.V0.Pos_Momentum.Py = sexa.V0.Pos_Py();
    fOutput_MC_ChannelD.V0.Pos_Momentum.Pz = sexa.V0.Pos_Pz();

    // Kaon (`Flat::MC_Track`)
    // -- `Flat::LV`
    fOutput_MC_ChannelD.Kaon.Px = sexa.Kaon.Px();
    fOutput_MC_ChannelD.Kaon.Py = sexa.Kaon.Py();
    fOutput_MC_ChannelD.Kaon.Pz = sexa.Kaon.Pz();
    fOutput_MC_ChannelD.Kaon.Energy = sexa.Kaon.Energy();
    // -- `Flat::MC`
    fOutput_MC_ChannelD.Kaon.McEntry = static_cast<int>(sexa.Kaon.Entry);
    fOutput_MC_ChannelD.Kaon.PdgCode = sexa.Kaon.PdgCode();
    fOutput_MC_ChannelD.Kaon.ReactionID = sexa.Kaon.ReactionID();
    fOutput_MC_ChannelD.Kaon.IsTrue = sexa.Kaon.IsTrue();
    fOutput_MC_ChannelD.Kaon.IsSignal = sexa.Kaon.IsSignal();
    fOutput_MC_ChannelD.Kaon.IsSecondary = sexa.Kaon.IsSecondary();
    // -- Mother (`Flat::MC_Id`)
    fOutput_MC_ChannelD.Kaon.Mother.McEntry = sexa.Kaon.Mother_McEntry();
    fOutput_MC_ChannelD.Kaon.Mother.PdgCode = sexa.Kaon.Mother_PdgCode();
    // -- GrandMother (`Flat::MC_Id`)
    fOutput_MC_ChannelD.Kaon.GrandMother.McEntry = sexa.Kaon.GrandMother_McEntry();
    fOutput_MC_ChannelD.Kaon.GrandMother.PdgCode = sexa.Kaon.GrandMother_PdgCode();
}

void Finder::StoreDummyMC_ChannelD() {
    // `Flat::MC_ChannelD` //

    // `Flat::MC_Sexaquark`
    // -- Before (`Flat::LV`)
    fOutput_MC_ChannelD.Before.Px = Const::DummyFloat;
    fOutput_MC_ChannelD.Before.Py = Const::DummyFloat;
    fOutput_MC_ChannelD.Before.Pz = Const::DummyFloat;
    fOutput_MC_ChannelD.Before.Energy = Const::DummyFloat;
    // -- After (`Flat::LV`)
    fOutput_MC_ChannelD.After.Px = Const::DummyFloat;
    fOutput_MC_ChannelD.After.Py = Const::DummyFloat;
    fOutput_MC_ChannelD.After.Pz = Const::DummyFloat;
    fOutput_MC_ChannelD.After.Energy = Const::DummyFloat;
    // -- Nucleon (`Flat::LV`)
    fOutput_MC_ChannelD.Nucleon.Px = Const::DummyFloat;
    fOutput_MC_ChannelD.Nucleon.Py = Const::DummyFloat;
    fOutput_MC_ChannelD.Nucleon.Pz = Const::DummyFloat;
    fOutput_MC_ChannelD.Nucleon.Energy = Const::DummyFloat;
    // -- PV (`Flat::Coordinates`)
    fOutput_MC_ChannelD.PV.X = Const::DummyFloat;
    fOutput_MC_ChannelD.PV.Y = Const::DummyFloat;
    fOutput_MC_ChannelD.PV.Z = Const::DummyFloat;
    // -- SV (`Flat::Coordinates`)
    fOutput_MC_ChannelD.SV.X = Const::DummyFloat;
    fOutput_MC_ChannelD.SV.Y = Const::DummyFloat;
    fOutput_MC_ChannelD.SV.Z = Const::DummyFloat;
    // -- Reaction ID + Flags
    fOutput_MC_ChannelD.ReactionID = Const::DummyInt;
    fOutput_MC_ChannelD.IsSignal = false;
    fOutput_MC_ChannelD.IsHybrid = false;

    // V0 (`Flat::MC_V0`)
    // -- `Flat::MC`
    fOutput_MC_ChannelD.V0.McEntry = Const::DummyInt;
    fOutput_MC_ChannelD.V0.PdgCode = Const::DummyInt;
    fOutput_MC_ChannelD.V0.ReactionID = Const::DummyInt;
    fOutput_MC_ChannelD.V0.IsTrue = false;
    fOutput_MC_ChannelD.V0.IsSignal = false;
    fOutput_MC_ChannelD.V0.IsSecondary = false;
    // -- `Flat::LV`
    fOutput_MC_ChannelD.V0.Px = Const::DummyFloat;
    fOutput_MC_ChannelD.V0.Py = Const::DummyFloat;
    fOutput_MC_ChannelD.V0.Pz = Const::DummyFloat;
    fOutput_MC_ChannelD.V0.Energy = Const::DummyFloat;
    // -- Mother (`Flat::MC_Id`)
    fOutput_MC_ChannelD.V0.Mother.McEntry = Const::DummyInt;
    fOutput_MC_ChannelD.V0.Mother.PdgCode = Const::DummyInt;
    // -- @ decay (`Flat::Coordinates`)
    fOutput_MC_ChannelD.V0.AtDecay.X = Const::DummyFloat;
    fOutput_MC_ChannelD.V0.AtDecay.Y = Const::DummyFloat;
    fOutput_MC_ChannelD.V0.AtDecay.Z = Const::DummyFloat;
    // -- hybrid flag
    fOutput_MC_ChannelD.V0.IsHybrid = false;

    // V0 Neg
    // -- `Flat::MC`
    fOutput_MC_ChannelD.V0.Neg.McEntry = Const::DummyInt;
    fOutput_MC_ChannelD.V0.Neg.PdgCode = Const::DummyInt;
    fOutput_MC_ChannelD.V0.Neg.ReactionID = Const::DummyInt;
    fOutput_MC_ChannelD.V0.Neg.IsTrue = false;
    fOutput_MC_ChannelD.V0.Neg.IsSignal = false;
    fOutput_MC_ChannelD.V0.Neg.IsSecondary = false;
    // -- `Flat::PxPyPz`
    fOutput_MC_ChannelD.V0.Neg_Momentum.Px = Const::DummyFloat;
    fOutput_MC_ChannelD.V0.Neg_Momentum.Py = Const::DummyFloat;
    fOutput_MC_ChannelD.V0.Neg_Momentum.Pz = Const::DummyFloat;

    // V0 Pos
    // -- `Flat::MC`
    fOutput_MC_ChannelD.V0.Pos.McEntry = Const::DummyInt;
    fOutput_MC_ChannelD.V0.Pos.PdgCode = Const::DummyInt;
    fOutput_MC_ChannelD.V0.Pos.ReactionID = Const::DummyInt;
    fOutput_MC_ChannelD.V0.Pos.IsTrue = false;
    fOutput_MC_ChannelD.V0.Pos.IsSignal = false;
    fOutput_MC_ChannelD.V0.Pos.IsSecondary = false;
    // -- `Flat::PxPyPz`
    fOutput_MC_ChannelD.V0.Pos_Momentum.Px = Const::DummyFloat;
    fOutput_MC_ChannelD.V0.Pos_Momentum.Py = Const::DummyFloat;
    fOutput_MC_ChannelD.V0.Pos_Momentum.Pz = Const::DummyFloat;

    // Kaon (`Flat::MC_Track`)
    // -- `Flat::LV`
    fOutput_MC_ChannelD.Kaon.Px = Const::DummyFloat;
    fOutput_MC_ChannelD.Kaon.Py = Const::DummyFloat;
    fOutput_MC_ChannelD.Kaon.Pz = Const::DummyFloat;
    fOutput_MC_ChannelD.Kaon.Energy = Const::DummyFloat;
    // -- `Flat::MC`
    fOutput_MC_ChannelD.Kaon.McEntry = Const::DummyInt;
    fOutput_MC_ChannelD.Kaon.PdgCode = Const::DummyInt;
    fOutput_MC_ChannelD.Kaon.ReactionID = Const::DummyInt;
    fOutput_MC_ChannelD.Kaon.IsTrue = false;
    fOutput_MC_ChannelD.Kaon.IsSignal = false;
    fOutput_MC_ChannelD.Kaon.IsSecondary = false;
    // -- Mother (`Flat::MC_Id`)
    fOutput_MC_ChannelD.Kaon.Mother.McEntry = Const::DummyInt;
    fOutput_MC_ChannelD.Kaon.Mother.PdgCode = Const::DummyInt;
    // -- GrandMother (`Flat::MC_Id`)
    fOutput_MC_ChannelD.Kaon.GrandMother.McEntry = Const::DummyInt;
    fOutput_MC_ChannelD.Kaon.GrandMother.PdgCode = Const::DummyInt;
}

// ## END OF CYCLES ## //

void Finder::EndOfAnalysis() {

    Logger::Info(__FUNCTION__, "The following objects have been written into TFile {}:", fSettings.PathOutputFile);

    // write trees //

    if (IsMC()) {
        fOutputTree_Injected->Write();
        Logger::Info(__FUNCTION__, "- TTree \"{}\"", fOutputTree_Injected->GetName());
    }
    fOutputTree->Write();
    Logger::Info(__FUNCTION__, "- TTree \"{}\"", fOutputTree->GetName());

    // write histograms //

    // -- event counter
    fHist_EventCounter->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_EventCounter->GetName());
    // -- cut flows
    fHist_CutFlow->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow->GetName());
    fHist_CutFlow_AntiChannel->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow_AntiChannel->GetName());

    // reset branch addresses //

    fInputChain_PackedEvents->ResetBranchAddresses();
    if (IsMC()) fOutputTree_Injected->ResetBranchAddresses();
    fOutputTree->ResetBranchAddresses();

    Logger::Info(__FUNCTION__, "All done.");
}

}  // namespace Tree2Secondaries
