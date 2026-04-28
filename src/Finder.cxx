#include "Finder/Finder.hxx"

#include <cstddef>
#include <filesystem>
#include <memory>
#include <string_view>

#include "App/Logger.hxx"
#include "Finder/FinderCuts.hxx"
#include "KalmanFitter/BaseKalmanFitter.hxx"
#include "KalmanFitter/KalmanFitterChannelA.hxx"
#include "KalmanFitter/KalmanFitterChannelD.hxx"
#include "Math/BaseMath.hxx"
#include "Math/Constants.hxx"
#include "Seeder/SeederHelixLine.hxx"
#include "Seeder/SeederLineLine.hxx"
#include "Storage/ReadWrite/ReadWriteSchemaFlat.hxx"
#include "Storage/ReadWrite/ReadWriteSchemaVector.hxx"
#include "Truth/TruthSexaquark.hxx"
#include "View/MC/ViewMcInjected.hxx"
#include "View/MC/ViewMcPackedTrack.hxx"
#include "View/MC/ViewMcPackedV0.hxx"
#include "View/Reconstructed/ViewRecTrack.hxx"

namespace Tree2Secondaries {

namespace KF = KalmanFitter;

bool Finder::Initialize() {

    fInputChain_PackedEvents = std::make_unique<TChain>(Const::TreeName_PackedEvents.data());
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
        IO::CreateBranches(fOutputTree_Injected.get(), fOutput_Injected);
    }

    Logger::Info(__FUNCTION__, "Finder initialized successfully.");

    return true;
}

// ## INPUT ZONE ## //

void Finder::ReadInputBranches() {
    // by default, turn off all branches
    fInputChain_PackedEvents->SetBranchStatus("*", false);
    // connect input branches to memory
    IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_Event, IsMC());
    if (IsMC()) IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_Injected, true);
    // -- depending on reaction channels
    switch (GetReactionChannel()) {
        case EReactionChannel::A:
            IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_AntiLambdas, Const::V0_Acronym[PID_V0::AntiLambda]);
            IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_Lambdas, Const::V0_Acronym[PID_V0::Lambda]);
            IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_KaonsZeroShort, Const::V0_Acronym[PID_V0::KaonZeroShort]);
            if (IsMC()) {
                IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_MC_AntiLambdas, Const::V0_Acronym[PID_V0::AntiLambda]);
                IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_MC_Lambdas, Const::V0_Acronym[PID_V0::Lambda]);
                IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_MC_KaonsZeroShort, Const::V0_Acronym[PID_V0::KaonZeroShort]);
            }
            break;
        case EReactionChannel::D:
            IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_AntiLambdas, Const::V0_Acronym[PID_V0::AntiLambda]);
            IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_Lambdas, Const::V0_Acronym[PID_V0::Lambda]);
            IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_NegKaons, false, Const::Particle_Acronym[PID_StableParticle::NegKaon]);
            IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_PosKaons, false, Const::Particle_Acronym[PID_StableParticle::PosKaon]);
            if (IsMC()) {
                IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_MC_AntiLambdas, Const::V0_Acronym[PID_V0::AntiLambda]);
                IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_MC_Lambdas, Const::V0_Acronym[PID_V0::Lambda]);
                IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_MC_NegKaons, true, Const::Particle_Acronym[PID_StableParticle::NegKaon]);
                IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_MC_PosKaons, true, Const::Particle_Acronym[PID_StableParticle::PosKaon]);
            }
            break;
        case EReactionChannel::E:
            IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_AntiLambdas, Const::V0_Acronym[PID_V0::AntiLambda]);
            IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_Lambdas, Const::V0_Acronym[PID_V0::Lambda]);
            IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_NegKaons, false, Const::Particle_Acronym[PID_StableParticle::NegKaon]);
            IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_PosKaons, false, Const::Particle_Acronym[PID_StableParticle::PosKaon]);
            IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_PiMinus, false, Const::Particle_Acronym[PID_StableParticle::PiMinus]);
            IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_PiPlus, false, Const::Particle_Acronym[PID_StableParticle::PiPlus]);
            if (IsMC()) {
                IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_MC_AntiLambdas, Const::V0_Acronym[PID_V0::AntiLambda]);
                IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_MC_Lambdas, Const::V0_Acronym[PID_V0::Lambda]);
                IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_MC_NegKaons, true, Const::Particle_Acronym[PID_StableParticle::NegKaon]);
                IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_MC_PosKaons, true, Const::Particle_Acronym[PID_StableParticle::PosKaon]);
                IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_MC_PiMinus, true, Const::Particle_Acronym[PID_StableParticle::PiMinus]);
                IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_MC_PiPlus, true, Const::Particle_Acronym[PID_StableParticle::PiPlus]);
            }
            break;
        case EReactionChannel::H:
            IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_NegKaons, false, Const::Particle_Acronym[PID_StableParticle::NegKaon]);
            IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_PosKaons, false, Const::Particle_Acronym[PID_StableParticle::PosKaon]);
            if (IsMC()) {
                IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_MC_NegKaons, true, Const::Particle_Acronym[PID_StableParticle::NegKaon]);
                IO::ReadBranches(fInputChain_PackedEvents.get(), fInput_MC_PosKaons, true, Const::Particle_Acronym[PID_StableParticle::PosKaon]);
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
    constexpr int x_nbins = 20;
    constexpr float x_min = 0.;
    constexpr float x_max = 20.;
    std::string_view hist_title = ";Cut N;N Passed Cut";
    fHist_CutFlow = std::make_unique<TH1D>("CutFlow", hist_title.data(), x_nbins, x_min, x_max);
    fHist_CutFlow_AntiChannel = std::make_unique<TH1D>("CutFlow_Anti", hist_title.data(), x_nbins, x_min, x_max);
}

bool Finder::PrepareOutputTree() {

    auto tree_name = std::format("FoundChannel{}", Const::ReactionChannel_Char[fSettings.ReactionChannel]);

    fOutputTree = std::make_unique<TTree>(tree_name.data(), "");
    if (!fOutputTree) {
        Logger::Error(__FUNCTION__, "Couldn't create TTree \"{}\"", tree_name);
        return false;
    }

    return true;
}

void Finder::CreateOutputBranches() {
    switch (GetReactionChannel()) {
        case EReactionChannel::A:
            IO::CreateBranches(fOutputTree.get(), fOutput_ChannelA, IsMC());
            if (IsMC()) IO::CreateBranches(fOutputTree.get(), fOutput_MC_ChannelA);
            break;
        case EReactionChannel::D:
            IO::CreateBranches(fOutputTree.get(), fOutput_ChannelD, IsMC());
            if (IsMC()) IO::CreateBranches(fOutputTree.get(), fOutput_MC_ChannelD);
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

    fOutputTree_Injected = std::make_unique<TTree>(Const::TreeName_Injected.data(), "");
    if (fOutputTree_Injected == nullptr) {
        Logger::Error(__FUNCTION__, "Couldn't create TTree \"{}\"", Const::TreeName_Injected);
        return false;
    }

    return true;
}

void Finder::ProcessInjected() {
    auto n_injected = static_cast<int>(NumberInjected());
    for (int entry = 0; entry < n_injected; ++entry) {

        View::MC::Injected sexa{&fInput_Injected, entry};  // NOTE: unguarded because contiguous access

        StoreInjected(sexa);

        fOutputTree_Injected->Fill();
    }
}

void Finder::StoreInjected(const View::MC::Injected& sexa) {
    double nucleon_mass = Const::Particle_Mass[Const::ReactionChannel_NucleonPID[fSettings.ReactionChannel]];  // PENDING: can prettify?

    fOutput_Injected.lv.px = sexa.Px();
    fOutput_Injected.lv.py = sexa.Py();
    fOutput_Injected.lv.pz = sexa.Pz();
    fOutput_Injected.lv.energy = static_cast<float>(Truth::Sexaquark::Energy(sexa, fSettings.SexaquarkMass));  // PENDING: can prettify?

    fOutput_Injected.sv.x = sexa.SV_X();
    fOutput_Injected.sv.y = sexa.SV_Y();
    fOutput_Injected.sv.z = sexa.SV_Z();

    fOutput_Injected.lv_nucleon.px = sexa.Nucleon_Px();
    fOutput_Injected.lv_nucleon.py = sexa.Nucleon_Py();
    fOutput_Injected.lv_nucleon.pz = sexa.Nucleon_Pz();
    fOutput_Injected.lv_nucleon.energy = static_cast<float>(Truth::Sexaquark::NucleonEnergy(sexa, nucleon_mass));  // PENDING: can prettify?

    fOutput_Injected.run_number = fInput_Event.run_number;
    fOutput_Injected.dir_number = fInput_Event.dir_number;
    fOutput_Injected.event_number = fInput_Event.event_number;
    fOutput_Injected.reaction_id = static_cast<unsigned int>(Truth::Sexaquark::ReactionID(sexa));  // PENDING: can prettify?
}

// ## Helpers ## //

void Finder::Assign(Schema::Flat::Track& out, const View::Rec::Track& t) {
    out.esd_entry = t.EsdEntry();
    out.charge = t.Charge<char>();
    out.state = {t.X(), t.Y(), t.Z(), t.Px(), t.Py(), t.Pz()};
    out.pre_dca_xy = t.PreDCAxy();
    out.pre_dca_z = t.PreDCAz();
    out.tpc_signal = t.TPC_Signal();
    out.n_sigma_pion = t.NSigmaPion();
    out.n_sigma_kaon = t.NSigmaKaon();
    out.n_sigma_proton = t.NSigmaProton();
}

void Finder::Assign(Schema::Flat::V0& out, const View::Rec::V0& v) {
    out.decay = {v.X(), v.Y(), v.Z()};
    out.lv = {v.Px(), v.Py(), v.Pz(), v.Energy()};
    out.chi2ndf = v.Chi2NDF();
    Assign(out.neg, v.Neg);
    out.neg_pca_v0 = {v.Neg_PCAwrtV0_X(), v.Neg_PCAwrtV0_Y(), v.Neg_PCAwrtV0_Z(), v.Neg_PCAwrtV0_Px(), v.Neg_PCAwrtV0_Py(), v.Neg_PCAwrtV0_Pz()};
    Assign(out.pos, v.Pos);
    out.pos_pca_v0 = {v.Pos_PCAwrtV0_X(), v.Pos_PCAwrtV0_Y(), v.Pos_PCAwrtV0_Z(), v.Pos_PCAwrtV0_Px(), v.Pos_PCAwrtV0_Py(), v.Pos_PCAwrtV0_Pz()};
}

void Finder::Assign(Schema::Flat::MC_Track& out, const View::MC::PackedTrack& t, bool ascendants_info) {
    out.mc_entry = t.McEntry();
    out.pdg_code = t.PdgCode();
    out.lv = {t.Px(), t.Py(), t.Pz(), t.Energy()};
    out.is_true = t.IsTrue();
    out.is_secondary = t.IsSecondary();
    out.is_signal = t.IsSignal();
    out.reaction_id = t.ReactionID();
    if (!ascendants_info) return;
    out.mother_mc_entry = t.Mother_McEntry();
    out.mother_pdg_code = t.Mother_PdgCode();
    out.gm_mc_entry = t.GrandMother_McEntry();
    out.gm_pdg_code = t.GrandMother_PdgCode();
}

void Finder::Assign(Schema::Flat::MC_V0& out, const View::MC::PackedV0& v) {
    out.lv = {v.Px(), v.Py(), v.Pz(), v.Energy()};
    out.origin = {v.Origin_X(), v.Origin_Y(), v.Origin_Z()};
    out.decay = {v.Decay_X(), v.Decay_Y(), v.Decay_Z()};
    out.mc_entry = v.McEntry();
    out.pdg_code = v.PdgCode();
    out.reaction_id = v.ReactionID();
    out.is_true = v.IsTrue();
    out.is_signal = v.IsSignal();
    out.is_secondary = v.IsSecondary();
    out.is_hybrid = v.IsHybrid();
    out.mother_mc_entry = v.Mother_McEntry();
    out.mother_pdg_code = v.Mother_PdgCode();
    // Assign(out.neg, v.Neg); // PENDING
    // Assign(out.pos, v.Pos); // PENDING
}

// ## Channel A ZONE ## //

void Finder::FindSexaquarks_ChannelA(bool anti_channel) {

    // determine properties based on anti-channel or not
    // -- v0a
    const Schema::Vector::V0s* Packed_V0A = &fInput_AntiLambdas;
    const Schema::Vector::MC_V0s* MC_V0A = &fInput_MC_AntiLambdas;
    if (anti_channel) {
        Packed_V0A = &fInput_Lambdas;
        MC_V0A = &fInput_MC_Lambdas;
    }
    // -- v0b
    const Schema::Vector::V0s* Packed_V0B = &fInput_KaonsZeroShort;
    const Schema::Vector::MC_V0s* MC_V0B = &fInput_MC_KaonsZeroShort;
    // -- cut flow hist
    TH1D* hist = anti_channel ? fHist_CutFlow_AntiChannel.get() : fHist_CutFlow.get();

    // loop over all possible pairs of (anti)lambda + K0S //
    std::size_t n_v0a = Packed_V0A->decay.x->size();
    std::size_t n_v0b = Packed_V0B->decay.x->size();
    for (unsigned int v0a_entry = 0; v0a_entry < n_v0a; ++v0a_entry) {
        for (unsigned int v0b_entry = 0; v0b_entry < n_v0b; ++v0b_entry) {

            // get views //
            View::Rec::V0 v0a{Packed_V0A, v0a_entry};
            View::Rec::V0 v0b{Packed_V0B, v0b_entry};

            // sanity check //
            if (v0a.Neg.EsdEntry() == v0b.Neg.EsdEntry() || v0a.Neg.EsdEntry() == v0b.Pos.EsdEntry() || v0a.Pos.EsdEntry() == v0b.Neg.EsdEntry() ||
                v0a.Pos.EsdEntry() == v0b.Pos.EsdEntry()) {
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
                    fOutput_MC_ChannelA = {};  // NOTE: dummify
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

    // if (sexa.CPA_V0A_SV() < Cuts::ChannelA::Min_La_CPAwrtSV) return false;
    // cut_flow_hist->Fill(7.);

    // if (sexa.CPA_V0B_SV() < Cuts::ChannelA::Min_K0S_CPAwrtSV) return false;
    // cut_flow_hist->Fill(8.);

    return true;
}

void Finder::Store(const KF::ChannelA& sexa, bool anti_channel) {
    fOutput_ChannelA.event = fInput_Event;

    fOutput_ChannelA.sv.x = static_cast<float>(sexa.X());
    fOutput_ChannelA.sv.y = static_cast<float>(sexa.Y());
    fOutput_ChannelA.sv.z = static_cast<float>(sexa.Z());
    fOutput_ChannelA.lv.px = static_cast<float>(sexa.Px());
    fOutput_ChannelA.lv.py = static_cast<float>(sexa.Py());
    fOutput_ChannelA.lv.pz = static_cast<float>(sexa.Pz());
    fOutput_ChannelA.lv.energy = static_cast<float>(sexa.E());
    fOutput_ChannelA.chi2ndf = static_cast<float>(sexa.Chi2NDF());
    fOutput_ChannelA.e_minus_nucleon = static_cast<float>(sexa.E_MinusNucleon());
    fOutput_ChannelA.is_anti = anti_channel;

    Assign(fOutput_ChannelA.v0a, sexa.V0A);
    Assign(fOutput_ChannelA.v0b, sexa.V0B);

    fOutput_ChannelA.v0a_pca_sv.x = static_cast<float>(sexa.V0A_PCAwrtSV.X());
    fOutput_ChannelA.v0a_pca_sv.y = static_cast<float>(sexa.V0A_PCAwrtSV.Y());
    fOutput_ChannelA.v0a_pca_sv.z = static_cast<float>(sexa.V0A_PCAwrtSV.Z());
    fOutput_ChannelA.v0a_pca_sv.px = static_cast<float>(sexa.V0A_PCAwrtSV.Px());
    fOutput_ChannelA.v0a_pca_sv.py = static_cast<float>(sexa.V0A_PCAwrtSV.Py());
    fOutput_ChannelA.v0a_pca_sv.pz = static_cast<float>(sexa.V0A_PCAwrtSV.Pz());

    fOutput_ChannelA.v0b_pca_sv.x = static_cast<float>(sexa.V0B_PCAwrtSV.X());
    fOutput_ChannelA.v0b_pca_sv.y = static_cast<float>(sexa.V0B_PCAwrtSV.Y());
    fOutput_ChannelA.v0b_pca_sv.z = static_cast<float>(sexa.V0B_PCAwrtSV.Z());
    fOutput_ChannelA.v0b_pca_sv.px = static_cast<float>(sexa.V0B_PCAwrtSV.Px());
    fOutput_ChannelA.v0b_pca_sv.py = static_cast<float>(sexa.V0B_PCAwrtSV.Py());
    fOutput_ChannelA.v0b_pca_sv.pz = static_cast<float>(sexa.V0B_PCAwrtSV.Pz());
}

void Finder::StoreMC(const View::MC::ChannelA& mc_sexa) {

    fOutput_MC_ChannelA.before.px = mc_sexa.Px();
    fOutput_MC_ChannelA.before.py = mc_sexa.Py();
    fOutput_MC_ChannelA.before.pz = mc_sexa.Pz();
    fOutput_MC_ChannelA.before.energy = static_cast<float>(Truth::Sexaquark::Energy(mc_sexa, fSettings.SexaquarkMass));

    fOutput_MC_ChannelA.after.px = static_cast<float>(Truth::Sexaquark::AfterPx(mc_sexa));
    fOutput_MC_ChannelA.after.py = static_cast<float>(Truth::Sexaquark::AfterPy(mc_sexa));
    fOutput_MC_ChannelA.after.pz = static_cast<float>(Truth::Sexaquark::AfterPz(mc_sexa));
    fOutput_MC_ChannelA.after.energy = static_cast<float>(Truth::Sexaquark::AfterE(mc_sexa));

    fOutput_MC_ChannelA.nucleon.px = mc_sexa.Nucleon_Px();
    fOutput_MC_ChannelA.nucleon.py = mc_sexa.Nucleon_Py();
    fOutput_MC_ChannelA.nucleon.pz = mc_sexa.Nucleon_Pz();
    fOutput_MC_ChannelA.nucleon.energy = static_cast<float>(
        Truth::Sexaquark::NucleonEnergy(mc_sexa, Const::Particle_Mass[Const::ReactionChannel_NucleonPID[fSettings.ReactionChannel]]));

    fOutput_MC_ChannelA.sv.x = mc_sexa.SV_X();
    fOutput_MC_ChannelA.sv.y = mc_sexa.SV_Y();
    fOutput_MC_ChannelA.sv.z = mc_sexa.SV_Z();

    fOutput_MC_ChannelA.reaction_id = Truth::Sexaquark::ReactionID(mc_sexa);
    fOutput_MC_ChannelA.is_signal = Truth::Sexaquark::IsSignal(mc_sexa);
    fOutput_MC_ChannelA.is_hybrid = Truth::Sexaquark::IsHybrid(mc_sexa);

    Assign(fOutput_MC_ChannelA.v0a, mc_sexa.V0A);
    Assign(fOutput_MC_ChannelA.v0b, mc_sexa.V0B);
}

// ## Channel D ZONE ## //

void Finder::FindSexaquarks_ChannelD(bool anti_channel) {

    // determine properties based on anti-channel or not
    // -- kaon
    const Schema::Vector::Tracks* Packed_Kaons = &fInput_PosKaons;
    const Schema::Vector::MC_Tracks* MC_Kaons = &fInput_MC_PosKaons;
    const double ka_mass = Const::Particle_Mass[PID_StableParticle::NegKaon];  // NOTE: same as mass of pos. kaon
    if (anti_channel) {
        Packed_Kaons = &fInput_NegKaons;
        MC_Kaons = &fInput_MC_NegKaons;
    }
    // -- v0
    const Schema::Vector::V0s* Packed_V0s = &fInput_AntiLambdas;
    const Schema::Vector::MC_V0s* MC_V0s = &fInput_MC_AntiLambdas;
    if (anti_channel) {
        Packed_V0s = &fInput_Lambdas;
        MC_V0s = &fInput_MC_Lambdas;
    }
    // -- cut flow hist
    TH1D* hist{anti_channel ? fHist_CutFlow_AntiChannel.get() : fHist_CutFlow.get()};

    // loop over all possible pairs of (anti)lambda + (pos/neg)kaon //
    std::size_t n_kaons = Packed_Kaons->state.x->size();
    std::size_t n_v0 = Packed_V0s->decay.x->size();
    for (unsigned int ka_entry = 0; ka_entry < n_kaons; ++ka_entry) {
        for (unsigned int v0_entry = 0; v0_entry < n_v0; ++v0_entry) {

            // get views //
            View::Rec::Track ka{Packed_Kaons, ka_entry};
            View::Rec::V0 v0{Packed_V0s, v0_entry};

            // sanity check //
            if (v0.Neg.EsdEntry() == ka.EsdEntry() || v0.Pos.EsdEntry() == ka.EsdEntry()) continue;

            // PCAs (1) //
            auto [seed_ka, seed_v0, pca_cache] = Seeder::HelixLine::FastCorrectPCAs(ka, v0, fInput_Event.magnetic_field);

            // apply cuts (1) //
            if (!FastCuts_ChannelD(seed_ka, seed_v0, hist)) continue;

            // PCAs derivatives //
            auto [deriv_v0, deriv_ka] = Seeder::HelixLine::ComputeDerivatives(seed_ka, seed_v0, pca_cache);

            // fit vertex //
            auto fit = KF::FitVertex(ka, v0, ka_mass,                           //
                                     {seed_ka, deriv_ka}, {seed_v0, deriv_v0},  //
                                     fInput_Event.magnetic_field);
            KF::ChannelD sexa(fit, seed_ka.pca, seed_v0.pca, ka, v0);

            // apply cuts (2) //
            if (!SlowCuts(sexa, hist)) continue;

            // store //
            Store(sexa, anti_channel);

            if (IsMC()) {
                View::MC::PackedTrack mc_kaon(MC_Kaons, ka_entry);
                View::MC::PackedV0 mc_v0(MC_V0s, v0_entry);
                View::MC::ChannelD mc_sexa(&fInput_Injected, mc_kaon, mc_v0);
                if (View::IsValid(mc_sexa)) {
                    StoreMC(mc_sexa);
                } else {
                    fOutput_MC_ChannelD = {};  // NOTE: dummify
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
    fOutput_ChannelD.event = fInput_Event;
    fOutput_ChannelD.sv.x = static_cast<float>(sexa.X());
    fOutput_ChannelD.sv.y = static_cast<float>(sexa.Y());
    fOutput_ChannelD.sv.z = static_cast<float>(sexa.Z());
    fOutput_ChannelD.lv.px = static_cast<float>(sexa.Px());
    fOutput_ChannelD.lv.py = static_cast<float>(sexa.Py());
    fOutput_ChannelD.lv.pz = static_cast<float>(sexa.Pz());
    fOutput_ChannelD.lv.energy = static_cast<float>(sexa.E());
    fOutput_ChannelD.chi2ndf = static_cast<float>(sexa.Chi2NDF());
    fOutput_ChannelD.e_minus_nucleon = static_cast<float>(sexa.E_MinusNucleon());
    fOutput_ChannelD.is_anti = anti_channel;

    Assign(fOutput_ChannelD.kaon, sexa.Kaon);
    fOutput_ChannelD.kaon_pca_sv.x = static_cast<float>(sexa.Kaon_PCAwrtSV.X());
    fOutput_ChannelD.kaon_pca_sv.y = static_cast<float>(sexa.Kaon_PCAwrtSV.Y());
    fOutput_ChannelD.kaon_pca_sv.z = static_cast<float>(sexa.Kaon_PCAwrtSV.Z());
    fOutput_ChannelD.kaon_pca_sv.px = static_cast<float>(sexa.Kaon_PCAwrtSV.Px());
    fOutput_ChannelD.kaon_pca_sv.py = static_cast<float>(sexa.Kaon_PCAwrtSV.Py());
    fOutput_ChannelD.kaon_pca_sv.pz = static_cast<float>(sexa.Kaon_PCAwrtSV.Pz());

    // V0
    Assign(fOutput_ChannelD.v0, sexa.V0);
    fOutput_ChannelD.v0_pca_sv.x = static_cast<float>(sexa.V0_PCAwrtSV.X());
    fOutput_ChannelD.v0_pca_sv.y = static_cast<float>(sexa.V0_PCAwrtSV.Y());
    fOutput_ChannelD.v0_pca_sv.z = static_cast<float>(sexa.V0_PCAwrtSV.Z());
    fOutput_ChannelD.v0_pca_sv.px = static_cast<float>(sexa.V0_PCAwrtSV.Px());
    fOutput_ChannelD.v0_pca_sv.py = static_cast<float>(sexa.V0_PCAwrtSV.Py());
    fOutput_ChannelD.v0_pca_sv.pz = static_cast<float>(sexa.V0_PCAwrtSV.Pz());
}

void Finder::StoreMC(const View::MC::ChannelD& sexa) {
    fOutput_MC_ChannelD.before.px = sexa.Px();
    fOutput_MC_ChannelD.before.py = sexa.Py();
    fOutput_MC_ChannelD.before.pz = sexa.Pz();
    fOutput_MC_ChannelD.before.energy = static_cast<float>(Truth::Sexaquark::Energy(sexa, fSettings.SexaquarkMass));
    fOutput_MC_ChannelD.after.px = static_cast<float>(Truth::Sexaquark::AfterPx(sexa));
    fOutput_MC_ChannelD.after.py = static_cast<float>(Truth::Sexaquark::AfterPy(sexa));
    fOutput_MC_ChannelD.after.pz = static_cast<float>(Truth::Sexaquark::AfterPz(sexa));
    fOutput_MC_ChannelD.after.energy = static_cast<float>(Truth::Sexaquark::AfterE(sexa));
    fOutput_MC_ChannelD.nucleon.px = sexa.Nucleon_Px();
    fOutput_MC_ChannelD.nucleon.py = sexa.Nucleon_Py();
    fOutput_MC_ChannelD.nucleon.pz = sexa.Nucleon_Pz();
    fOutput_MC_ChannelD.nucleon.energy =
        static_cast<float>(Truth::Sexaquark::NucleonEnergy(sexa, Const::Particle_Mass[Const::ReactionChannel_NucleonPID[fSettings.ReactionChannel]]));
    fOutput_MC_ChannelD.sv.x = sexa.SV_X();
    fOutput_MC_ChannelD.sv.y = sexa.SV_Y();
    fOutput_MC_ChannelD.sv.z = sexa.SV_Z();
    fOutput_MC_ChannelD.reaction_id = Truth::Sexaquark::ReactionID(sexa);
    fOutput_MC_ChannelD.is_signal = Truth::Sexaquark::IsSignal(sexa);
    fOutput_MC_ChannelD.is_hybrid = Truth::Sexaquark::IsHybrid(sexa);

    Assign(fOutput_MC_ChannelD.v0, sexa.V0);
    Assign(fOutput_MC_ChannelD.kaon, sexa.Kaon, true);
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
