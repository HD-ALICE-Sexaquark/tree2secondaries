#include "Finder/Finder.hxx"

#include <cstddef>
#include <filesystem>
#include <memory>
#include <optional>
#include <string_view>

#include "common/Constants.hpp"
#include "common/DB_Particles.hpp"
#include "common/DB_ReactionChannels.hpp"
#include "common/FL_InjectedSexa.hpp"
#include "common/FL_Sexaquark.hpp"
#include "common/Math.hpp"
#include "common/R2DS_Cuts.hpp"
#include "common/VC_InjectedSexa.hpp"
#include "common/VC_InjectedSexaView.hpp"

#include "App/Logger.hxx"
#include "KalmanFitter/BaseKalmanFitter.hxx"
#include "KalmanFitter/KalmanFitterChannelA.hxx"
#include "KalmanFitter/KalmanFitterChannelD.hxx"
#include "Seeder/SeederHelixHelix.hxx"
#include "Seeder/SeederHelixLine.hxx"
#include "Seeder/SeederLineLine.hxx"

namespace R2DS {

namespace KF = KalmanFitter;

bool Finder::Initialize() {

    // Prepare Input Reader //
    // PENDING: could refactor into a function

    fInput_Model = ROOT::RNTupleModel::Create();

    fInput_Event.AddFieldsTo(fInput_Model.get(), fSettings.IsMC);

    if (fSettings.IsMC) {
        fInput_InjectedSexa.AddFieldsTo(fInput_Model.get(), true);
    }

    switch (fSettings.ReactionChannel.name) {
        case 'A': {
            fInput_AntiLambda.AddFieldsTo(fInput_Model.get(), fSettings.IsMC, "AL");
            fInput_Lambda.AddFieldsTo(fInput_Model.get(), fSettings.IsMC, "L");
            fInput_KaonZeroShort.AddFieldsTo(fInput_Model.get(), fSettings.IsMC, "K0S");
            break;
        }
        case 'D': {
            fInput_AntiLambda.AddFieldsTo(fInput_Model.get(), fSettings.IsMC, "AL");
            fInput_Lambda.AddFieldsTo(fInput_Model.get(), fSettings.IsMC, "L");
            break;
        }
        case 'H': {
            fInput_NegKaon.AddFieldsTo(fInput_Model.get(), fSettings.IsMC, fSettings.IsMC, "NK");
            fInput_PosKaon.AddFieldsTo(fInput_Model.get(), fSettings.IsMC, fSettings.IsMC, "PK");
            break;
        }
        default:
            break;
    }

    fReader = ROOT::RNTupleReader::Open(std::move(fInput_Model), R2DS::Name_PackedRNT,
                                        fSettings.PathInputFiles.front());  // PENDING: handle multiple files?

    // Prepare Output //

    fName_FoundRNT = std::format("Found_Channel{}", fSettings.ReactionChannel.name);

    if (!PrepareOutputFile()) return false;

    PrepareOutputHistograms();

    // Prepare Output Writers //
    // PENDING: could refactor into a function

    // -- reconstructed
    fOutput_Model = ROOT::RNTupleModel::Create();
    if (fSettings.ReactionChannel.name == 'A') {
        fOutput_ChannelA.AddFieldsTo(fOutput_Model.get(), fSettings.IsMC);
    }
    if (fSettings.ReactionChannel.name == 'D') {
        fOutput_ChannelD.AddFieldsTo(fOutput_Model.get(), fSettings.IsMC);
    }
    if (fSettings.ReactionChannel.name == 'H') {
        fOutput_ChannelH.AddFieldsTo(fOutput_Model.get(), fSettings.IsMC);
    }
    fWriter = ROOT::RNTupleWriter::Append(std::move(fOutput_Model), fName_FoundRNT, *fOutput_File);
    // -- mc
    if (fSettings.IsMC) {
        fOutput_Model_InjectedSexa = ROOT::RNTupleModel::Create();
        fOutput_InjectedSexa.AddFieldsTo(fOutput_Model_InjectedSexa.get(), false);
        fWriter_InjectedSexa = ROOT::RNTupleWriter::Append(std::move(fOutput_Model_InjectedSexa), R2DS::Name_InjectedSexaRNT, *fOutput_File);
    }

    Logger::Info(__FUNCTION__, "Finder initialized successfully.");

    return true;
}

// ## OUTPUT ZONE ## //

bool Finder::PrepareOutputFile() {

    const std::filesystem::path output_path(fSettings.PathOutputFile);
    if (output_path.has_parent_path()) std::filesystem::create_directories(output_path.parent_path());

    fOutput_File = std::unique_ptr<TFile>(TFile::Open(fSettings.PathOutputFile.c_str(), "RECREATE"));
    if (fOutput_File == nullptr || fOutput_File->IsZombie()) {
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
    constexpr const char* hist_title = ";Cut N;N Passed Cut";
    fHist_CutFlow = std::make_unique<TH1D>("CutFlow", hist_title, x_nbins, x_min, x_max);
    fHist_CutFlow_AntiChannel = std::make_unique<TH1D>("CutFlow_Anti", hist_title, x_nbins, x_min, x_max);
}

// ## Event ZONE ## //

void Finder::ProcessEvent() {  //
    fHist_EventCounter->Fill(0.);
    fPrimaryVertex.SetCoordinates(*fInput_Event.PV_X, *fInput_Event.PV_Y, *fInput_Event.PV_Z);
}

// ## Injected ZONE ## //

void Finder::ProcessInjected() {
    const double n_mass = DB::Particles::FindParticle(fSettings.ReactionChannel.nucleon_pdg)->mass;

    Vector::InjectedSexaView injected_view(&fInput_InjectedSexa);

    for (std::size_t entry_inj = 0; entry_inj < injected_view.Size(); ++entry_inj) {
        injected_view.CacheCalculations(entry_inj, fSettings.SexaquarkMass, n_mass);  // PENDING: why cache, should just read...
        // assign values //
        Assign_InjectedSexa(injected_view, fOutput_InjectedSexa, false);
        // fill rntuple //
        fWriter_InjectedSexa->Fill();
    }
}

// ## Helpers ## //

void Finder::Assign_Event(Flat::Sexaquark& out) {
    *out.RunNumber = *fInput_Event.RunNumber;
    *out.DirNumber = *fInput_Event.DirNumber;
    if (!fSettings.IsMC) *out.DirNumberB = *fInput_Event.DirNumberB;
    *out.EventNumber = *fInput_Event.EventNumber;
    *out.Centrality = *fInput_Event.Centrality;
    *out.MagneticField = *fInput_Event.MagneticField;
    *out.PV_X = *fInput_Event.PV_X;
    *out.PV_Y = *fInput_Event.PV_Y;
    *out.PV_Z = *fInput_Event.PV_Z;
    if (fSettings.IsMC) {
        *out.MC_PV_X = *fInput_Event.MC_PV_X;
        *out.MC_PV_Y = *fInput_Event.MC_PV_Y;
        *out.MC_PV_Z = *fInput_Event.MC_PV_Z;
    }
}

void Finder::Assign_Candidate(const KF::Particle& sexa, Flat::Sexaquark& out, bool anti_channel) {
    *out.SV_X = static_cast<float>(sexa.X());
    *out.SV_Y = static_cast<float>(sexa.Y());
    *out.SV_Z = static_cast<float>(sexa.Z());
    *out.Px = static_cast<float>(sexa.Px());
    *out.Py = static_cast<float>(sexa.Py());
    *out.Pz = static_cast<float>(sexa.Pz());
    *out.Energy = static_cast<float>(sexa.E());
    *out.Chi2NDF = static_cast<float>(sexa.Chi2NDF());
    *out.E_MinusNucleon = static_cast<float>(sexa.E() - DB::Particles::FindParticle(fSettings.ReactionChannel.nucleon_pdg)->mass);
    *out.AntiChannel = anti_channel;
}

void Finder::Assign_InjectedSexa(const std::optional<Vector::InjectedSexaView>& injected, Flat::InjectedSexa& out, bool embedded_to_rec) {
    if (!embedded_to_rec) {
        *out.RunNumber = *fInput_Event.RunNumber;
        *out.DirNumber = *fInput_Event.DirNumber;
        *out.EventNumber = *fInput_Event.EventNumber;
    }
    if (!injected.has_value()) {
        *out.ReactionID = static_cast<int>(injected->ReactionID());
        *out.SV_X = injected->SV_X();
        *out.SV_Y = injected->SV_Y();
        *out.SV_Z = injected->SV_Z();
        *out.Px = injected->Px();
        *out.Py = injected->Py();
        *out.Pz = injected->Pz();
        *out.Energy = static_cast<float>(injected->Energy());
        *out.Nucleon_Px = injected->Nucleon_Px();
        *out.Nucleon_Py = injected->Nucleon_Py();
        *out.Nucleon_Pz = injected->Nucleon_Pz();
        *out.Nucleon_Energy = static_cast<float>(injected->Nucleon_Energy());
    } else {
        *out.ReactionID = Common::DummyInt;
        *out.SV_X = Common::DummyFloat;
        *out.SV_Y = Common::DummyFloat;
        *out.SV_Z = Common::DummyFloat;
        *out.Px = Common::DummyFloat;
        *out.Py = Common::DummyFloat;
        *out.Pz = Common::DummyFloat;
        *out.Energy = Common::DummyFloat;
        *out.Nucleon_Px = Common::DummyFloat;
        *out.Nucleon_Py = Common::DummyFloat;
        *out.Nucleon_Pz = Common::DummyFloat;
        *out.Nucleon_Energy = Common::DummyFloat;
    }
}

void Finder::Assign(const Vector::TrackView& track, Flat::Track& out, bool include_gm) {
    *out.EsdEntry = track.EsdEntry();
    *out.X = track.X();
    *out.Y = track.Y();
    *out.Z = track.Z();
    *out.Px = track.Px();
    *out.Py = track.Py();
    *out.Pz = track.Pz();
    *out.PreDCAxy = track.PreDCAxy();
    *out.PreDCAz = track.PreDCAz();
    *out.TPC_Signal = track.TPC_Signal();
    *out.NSigmaPion = track.NSigmaPion();
    *out.NSigmaKaon = track.NSigmaKaon();
    *out.NSigmaProton = track.NSigmaProton();
    if (fSettings.IsMC) *out.McEntry = track.McEntry();
    // MC zone
    if (!fSettings.IsMC) return;
    *out.PdgCode = track.PdgCode();
    *out.Origin_X = track.Origin_X();
    *out.Origin_Y = track.Origin_Y();
    *out.Origin_Z = track.Origin_Z();
    *out.MC_Px = track.MC_Px();
    *out.MC_Py = track.MC_Py();
    *out.MC_Pz = track.MC_Pz();
    *out.MC_Energy = track.MC_Energy();
    *out.IsTrue = track.IsTrue();
    *out.IsSecondary = track.IsSecondary();
    *out.IsSignal = track.IsSignal();
    *out.ReactionID = track.ReactionID();
    *out.Mother_McEntry = track.Mother_McEntry();
    *out.Mother_PdgCode = track.Mother_PdgCode();
    if (!include_gm) return;
    *out.GM_McEntry = track.GM_McEntry();
    *out.GM_PdgCode = track.GM_PdgCode();
}

void Finder::Assign(const Vector::V0View& v0, Flat::V0& out) {
    *out.Decay_X = v0.Decay_X();
    *out.Decay_Y = v0.Decay_Y();
    *out.Decay_Z = v0.Decay_Z();
    *out.Px = v0.Px();
    *out.Py = v0.Py();
    *out.Pz = v0.Pz();
    *out.Energy = v0.Energy();
    *out.Chi2NDF = v0.Chi2NDF();
    // composition
    Assign(v0.Neg, out.Neg, false);
    Assign(v0.Pos, out.Pos, false);
    // MC zone
    if (!fSettings.IsMC) return;
    *out.McEntry = v0.McEntry();
    *out.PdgCode = v0.PdgCode();
    *out.Origin_X = v0.Origin_X();
    *out.Origin_Y = v0.Origin_Y();
    *out.Origin_Z = v0.Origin_Z();
    *out.MC_Decay_X = v0.MC_Decay_X();
    *out.MC_Decay_Y = v0.MC_Decay_Y();
    *out.MC_Decay_Z = v0.MC_Decay_Z();
    *out.MC_Px = v0.MC_Px();
    *out.MC_Py = v0.MC_Py();
    *out.MC_Pz = v0.MC_Pz();
    *out.MC_Energy = v0.MC_Energy();
    *out.IsTrue = v0.IsTrue();
    *out.IsSecondary = v0.IsSecondary();
    *out.IsSignal = v0.IsSignal();
    *out.IsHybrid = v0.IsHybrid();
    *out.ReactionID = v0.ReactionID();
    *out.Mother_McEntry = v0.Mother_McEntry();
    *out.Mother_PdgCode = v0.Mother_PdgCode();
}

// ## Channel A ZONE ## //

void Finder::FindSexaquarks_ChannelA(bool anti_channel) {

    // determine anti-channel or not //
    Vector::V0* Input_V0A = anti_channel ? &fInput_Lambda : &fInput_AntiLambda;
    TH1D* hist = anti_channel ? fHist_CutFlow_AntiChannel.get() : fHist_CutFlow.get();

    // create views //
    Vector::V0View v0a(Input_V0A);
    Vector::V0View v0b(&fInput_KaonZeroShort);

    // loop over all possible pairs of (anti)lambda + K0S //
    for (std::size_t entry_v0a = 0; entry_v0a < v0a.Size(); ++entry_v0a) {
        // -- cache V0A
        v0a.CacheCalculations(entry_v0a, fPrimaryVertex);

        for (std::size_t entry_v0b = 0; entry_v0b < v0b.Size(); ++entry_v0b) {
            // -- sanity check
            if (v0a.Neg.EsdEntry() == v0b.Neg.EsdEntry() || v0a.Neg.EsdEntry() == v0b.Pos.EsdEntry() || v0a.Pos.EsdEntry() == v0b.Neg.EsdEntry() ||
                v0a.Pos.EsdEntry() == v0b.Pos.EsdEntry()) {
                continue;
            }
            // -- cache V0B
            v0b.CacheCalculations(entry_v0b, fPrimaryVertex);

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
            Assign(sexa, anti_channel);

            // fill //
            fWriter->Fill();
        }
    }
}

bool Finder::FastCuts_ChannelA(const Seeder::Seed& seed_v0a, const Seeder::Seed& seed_v0b, TH1D* cut_flow_hist) const {

    cut_flow_hist->Fill(0.);

    if (Common::Math::SquaredDistance(seed_v0a.pca.xyz, seed_v0b.pca.xyz) >
        R2DS::Cuts::ChannelA::Max_DCAbtwV0s * R2DS::Cuts::ChannelA::Max_DCAbtwV0s) {
        return false;
    }
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

std::optional<Vector::InjectedSexaView> Finder::LinkInjectedSignal(const KF::ChannelA& sexa) {
    if (!sexa.V0A.IsSignal()) return std::nullopt;
    if (!sexa.V0B.IsSignal()) return std::nullopt;
    if (sexa.V0A.ReactionID() != sexa.V0B.ReactionID()) return std::nullopt;
    return Vector::InjectedSexaView{&fInput_InjectedSexa, static_cast<unsigned int>(sexa.V0A.ReactionID() - E2R::ReactionID_Offset)};
}

void Finder::Assign(const KF::ChannelA& sexa, bool anti_channel) {
    // event info
    Assign_Event(fOutput_ChannelA);
    // candidate info
    Assign_Candidate(sexa, fOutput_ChannelA, anti_channel);
    // -- V0A
    Assign(sexa.V0A, fOutput_ChannelA.V0A);
    *fOutput_ChannelA.V0A_PCAwrtSV_X = static_cast<float>(sexa.V0A_PCAwrtSV.X());
    *fOutput_ChannelA.V0A_PCAwrtSV_Y = static_cast<float>(sexa.V0A_PCAwrtSV.Y());
    *fOutput_ChannelA.V0A_PCAwrtSV_Z = static_cast<float>(sexa.V0A_PCAwrtSV.Z());
    *fOutput_ChannelA.V0A_PCAwrtSV_Px = static_cast<float>(sexa.V0A_PCAwrtSV.Px());
    *fOutput_ChannelA.V0A_PCAwrtSV_Py = static_cast<float>(sexa.V0A_PCAwrtSV.Py());
    *fOutput_ChannelA.V0A_PCAwrtSV_Pz = static_cast<float>(sexa.V0A_PCAwrtSV.Pz());
    // -- V0B
    Assign(sexa.V0B, fOutput_ChannelA.V0B);
    *fOutput_ChannelA.V0B_PCAwrtSV_X = static_cast<float>(sexa.V0B_PCAwrtSV.X());
    *fOutput_ChannelA.V0B_PCAwrtSV_Y = static_cast<float>(sexa.V0B_PCAwrtSV.Y());
    *fOutput_ChannelA.V0B_PCAwrtSV_Z = static_cast<float>(sexa.V0B_PCAwrtSV.Z());
    *fOutput_ChannelA.V0B_PCAwrtSV_Px = static_cast<float>(sexa.V0B_PCAwrtSV.Px());
    *fOutput_ChannelA.V0B_PCAwrtSV_Py = static_cast<float>(sexa.V0B_PCAwrtSV.Py());
    *fOutput_ChannelA.V0B_PCAwrtSV_Pz = static_cast<float>(sexa.V0B_PCAwrtSV.Pz());
    // -- MC zone
    if (!fSettings.IsMC) return;
    Assign_InjectedSexa(LinkInjectedSignal(sexa), fOutput_ChannelA.MC, true);
}

// ## Channel D ZONE ## //

void Finder::FindSexaquarks_ChannelD(bool anti_channel) {

    // determine anti-channel or not //
    Vector::V0* Input_V0s = anti_channel ? &fInput_Lambda : &fInput_AntiLambda;
    Vector::Track* Input_Kaons = anti_channel ? &fInput_PosKaon : &fInput_NegKaon;
    TH1D* hist = anti_channel ? fHist_CutFlow_AntiChannel.get() : fHist_CutFlow.get();

    // create views //
    Vector::V0View v0(Input_V0s);
    Vector::TrackView kaon(Input_Kaons);

    // loop over all possible pairs of (anti)lambda + (pos/neg)kaon //
    for (std::size_t entry_v0 = 0; entry_v0 < v0.Size(); ++entry_v0) {
        // -- cache V0
        v0.CacheCalculations(entry_v0, fPrimaryVertex);

        for (std::size_t entry_kaon = 0; entry_kaon < kaon.Size(); ++entry_kaon) {
            // -- sanity check
            if (v0.Neg.EsdEntry() == kaon.EsdEntry() || v0.Pos.EsdEntry() == kaon.EsdEntry()) continue;
            // -- cache kaon
            kaon.CacheCalculations(entry_kaon, fPrimaryVertex, *fInput_Event.MagneticField);

            // PCAs (1) //
            auto [seed_kaon, seed_v0, pca_cache] = Seeder::HelixLine::FastCorrectPCAs(kaon, v0, *fInput_Event.MagneticField);

            // apply cuts (1) //
            if (!FastCuts_ChannelD(seed_kaon, seed_v0, hist)) continue;

            // PCAs derivatives //
            auto [deriv_v0, deriv_ka] = Seeder::HelixLine::ComputeDerivatives(seed_kaon, seed_v0, pca_cache);

            // fit vertex //
            auto fit = KF::FitVertex(kaon, v0, DB::Particles::Particle("PosKaon").mass,  //
                                     {seed_kaon, deriv_ka}, {seed_v0, deriv_v0},         //
                                     *fInput_Event.MagneticField);
            KF::ChannelD sexa(fit, seed_v0.pca, seed_kaon.pca, v0, kaon);

            // apply cuts (2) //
            if (!SlowCuts(sexa, hist)) continue;

            // store //
            Assign(sexa, anti_channel);

            // fill //
            fWriter->Fill();
        }
    }
}

bool Finder::FastCuts_ChannelD(const Seeder::Seed& seed_ka, const Seeder::Seed& seed_v0, TH1D* cut_flow_hist) const {

    cut_flow_hist->Fill(0.);

    if (Common::Math::SquaredDistance(seed_ka.pca.xyz, seed_v0.pca.xyz) > Cuts::ChannelD::Max_DCAKaLa * Cuts::ChannelD::Max_DCAKaLa) return false;
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

std::optional<Vector::InjectedSexaView> Finder::LinkInjectedSignal(const KF::ChannelD& sexa) {
    if (!sexa.V0.IsSignal()) return std::nullopt;
    if (!sexa.Kaon.IsSignal()) return std::nullopt;
    if (sexa.V0.ReactionID() != sexa.Kaon.ReactionID()) return std::nullopt;
    return Vector::InjectedSexaView{&fInput_InjectedSexa, static_cast<unsigned int>(sexa.V0.ReactionID() - E2R::ReactionID_Offset)};
}

void Finder::Assign(const KF::ChannelD& sexa, bool anti_channel) {
    // -- event info
    Assign_Event(fOutput_ChannelD);
    // -- candidate info
    Assign_Candidate(sexa, fOutput_ChannelD, anti_channel);
    // -- V0
    Assign(sexa.V0, fOutput_ChannelD.V0);
    *fOutput_ChannelD.V0_PCAwrtSV_X = static_cast<float>(sexa.V0_PCAwrtSV.X());
    *fOutput_ChannelD.V0_PCAwrtSV_Y = static_cast<float>(sexa.V0_PCAwrtSV.Y());
    *fOutput_ChannelD.V0_PCAwrtSV_Z = static_cast<float>(sexa.V0_PCAwrtSV.Z());
    *fOutput_ChannelD.V0_PCAwrtSV_Px = static_cast<float>(sexa.V0_PCAwrtSV.Px());
    *fOutput_ChannelD.V0_PCAwrtSV_Py = static_cast<float>(sexa.V0_PCAwrtSV.Py());
    *fOutput_ChannelD.V0_PCAwrtSV_Pz = static_cast<float>(sexa.V0_PCAwrtSV.Pz());
    // -- kaon
    Assign(sexa.Kaon, fOutput_ChannelD.Kaon, true);
    *fOutput_ChannelD.Kaon_PCAwrtSV_X = static_cast<float>(sexa.Kaon_PCAwrtSV.X());
    *fOutput_ChannelD.Kaon_PCAwrtSV_Y = static_cast<float>(sexa.Kaon_PCAwrtSV.Y());
    *fOutput_ChannelD.Kaon_PCAwrtSV_Z = static_cast<float>(sexa.Kaon_PCAwrtSV.Z());
    *fOutput_ChannelD.Kaon_PCAwrtSV_Px = static_cast<float>(sexa.Kaon_PCAwrtSV.Px());
    *fOutput_ChannelD.Kaon_PCAwrtSV_Py = static_cast<float>(sexa.Kaon_PCAwrtSV.Py());
    *fOutput_ChannelD.Kaon_PCAwrtSV_Pz = static_cast<float>(sexa.Kaon_PCAwrtSV.Pz());
    // -- MC zone
    if (!fSettings.IsMC) return;
    Assign_InjectedSexa(LinkInjectedSignal(sexa), fOutput_ChannelD.MC, true);
}

// ## Channel H ZONE ## //

void Finder::FindSexaquarks_ChannelH(bool anti_channel) {

    constexpr double mass_kaon = DB::Particles::Particle("PosKaon").mass;

    // determine anti-channel or not //
    const Vector::Track* Input_Kaon = anti_channel ? &fInput_PosKaon : &fInput_NegKaon;
    TH1D* hist = anti_channel ? fHist_CutFlow_AntiChannel.get() : fHist_CutFlow.get();

    // create views //
    Vector::TrackView kaon1(Input_Kaon);
    Vector::TrackView kaon2(Input_Kaon);

    // loop over all possible pairs of (pos)kaon+(pos)kaon or (neg)kaon+(neg)kaon //
    for (std::size_t entry_kaon1 = 0; entry_kaon1 + 1 < kaon1.Size(); ++entry_kaon1) {
        // cache kaon1 //
        kaon1.CacheCalculations(entry_kaon1, fPrimaryVertex, *fInput_Event.MagneticField);

        for (std::size_t entry_kaon2 = entry_kaon1 + 1; entry_kaon2 < kaon2.Size(); ++entry_kaon2) {
            // sanity check //
            if (kaon1.EsdEntry() == kaon2.EsdEntry()) continue;
            // cache kaon2 //
            kaon2.CacheCalculations(entry_kaon2, fPrimaryVertex, *fInput_Event.MagneticField);

            // PCAs (1) //
            auto [seed_kaon1, seed_kaon2, pca_cache] = Seeder::HelixHelix::FastCorrectPCAs(kaon1, kaon2, *fInput_Event.MagneticField);

            // apply cuts (1) //
            if (!FastCuts_ChannelH(seed_kaon1, seed_kaon2, hist)) continue;

            // PCAs derivatives //
            auto [deriv_kaon1, deriv_kaon2] = Seeder::HelixHelix::ComputeDerivatives(seed_kaon1, seed_kaon2, pca_cache);

            // fit vertex //
            auto fit = KF::FitVertex(kaon1, kaon2, mass_kaon, mass_kaon,  //
                                     {seed_kaon1, deriv_kaon1}, {seed_kaon2, deriv_kaon2}, *fInput_Event.MagneticField);
            KF::ChannelH sexa(fit, seed_kaon1.pca, seed_kaon2.pca, kaon1, kaon2);

            // apply cuts (2) //
            if (!SlowCuts(sexa, hist)) continue;

            // store //
            Assign(sexa, anti_channel);

            // fill //
            fWriter->Fill();
        }
    }
}

bool Finder::FastCuts_ChannelH(const Seeder::Seed& seed_kaon1, const Seeder::Seed& seed_kaon2, TH1D* cut_flow_hist) const {
    //
    // PENDING
    //
    return true;
}

bool Finder::SlowCuts(const KF::ChannelH& sexa, TH1D* cut_flow_hist) const {
    //
    // PENDING
    //
    return true;
}

std::optional<Vector::InjectedSexaView> Finder::LinkInjectedSignal(const KF::ChannelH& sexa) {
    if (!sexa.Kaon1.IsSignal()) return std::nullopt;
    if (!sexa.Kaon2.IsSignal()) return std::nullopt;
    if (sexa.Kaon1.ReactionID() != sexa.Kaon2.ReactionID()) return std::nullopt;
    return Vector::InjectedSexaView{&fInput_InjectedSexa, static_cast<unsigned int>(sexa.Kaon1.ReactionID() - E2R::ReactionID_Offset)};
}

void Finder::Assign(const KF::ChannelH& sexa, bool anti_channel) {
    // -- event info
    Assign_Event(fOutput_ChannelH);
    // -- candidate info
    Assign_Candidate(sexa, fOutput_ChannelH, anti_channel);
    // -- kaon 1
    Assign(sexa.Kaon1, fOutput_ChannelH.Kaon1, true);
    *fOutput_ChannelH.Kaon1_PCAwrtSV_X = static_cast<float>(sexa.Kaon1_PCAwrtSV.X());
    *fOutput_ChannelH.Kaon1_PCAwrtSV_Y = static_cast<float>(sexa.Kaon1_PCAwrtSV.Y());
    *fOutput_ChannelH.Kaon1_PCAwrtSV_Z = static_cast<float>(sexa.Kaon1_PCAwrtSV.Z());
    *fOutput_ChannelH.Kaon1_PCAwrtSV_Px = static_cast<float>(sexa.Kaon1_PCAwrtSV.Px());
    *fOutput_ChannelH.Kaon1_PCAwrtSV_Py = static_cast<float>(sexa.Kaon1_PCAwrtSV.Py());
    *fOutput_ChannelH.Kaon1_PCAwrtSV_Pz = static_cast<float>(sexa.Kaon1_PCAwrtSV.Pz());
    // -- kaon 2
    Assign(sexa.Kaon2, fOutput_ChannelH.Kaon2, true);
    *fOutput_ChannelH.Kaon2_PCAwrtSV_X = static_cast<float>(sexa.Kaon2_PCAwrtSV.X());
    *fOutput_ChannelH.Kaon2_PCAwrtSV_Y = static_cast<float>(sexa.Kaon2_PCAwrtSV.Y());
    *fOutput_ChannelH.Kaon2_PCAwrtSV_Z = static_cast<float>(sexa.Kaon2_PCAwrtSV.Z());
    *fOutput_ChannelH.Kaon2_PCAwrtSV_Px = static_cast<float>(sexa.Kaon2_PCAwrtSV.Px());
    *fOutput_ChannelH.Kaon2_PCAwrtSV_Py = static_cast<float>(sexa.Kaon2_PCAwrtSV.Py());
    *fOutput_ChannelH.Kaon2_PCAwrtSV_Pz = static_cast<float>(sexa.Kaon2_PCAwrtSV.Pz());
    // -- MC zone
    if (!fSettings.IsMC) return;
    Assign_InjectedSexa(LinkInjectedSignal(sexa), fOutput_ChannelH.MC, true);
}

// ## END OF CYCLES ## //

void Finder::EndOfAnalysis() {

    Logger::Info(__FUNCTION__, "The following objects have been written into TFile {}:", fSettings.PathOutputFile);

    if (fSettings.IsMC) {
        Logger::Info(__FUNCTION__, "- RNTuple \"{}\"", R2DS::Name_InjectedSexaRNT);
    }
    Logger::Info(__FUNCTION__, "- RNTuple \"{}\"", fName_FoundRNT);

    // write histograms //

    // -- event counter
    fHist_EventCounter->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_EventCounter->GetName());
    // -- cut flows
    fHist_CutFlow->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow->GetName());
    fHist_CutFlow_AntiChannel->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow_AntiChannel->GetName());

    Logger::Info(__FUNCTION__, "All done.");
}

}  // namespace R2DS
