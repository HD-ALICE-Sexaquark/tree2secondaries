#include "Packager/Packager.hxx"

#include <cstddef>

#include "common/Constants.hpp"
#include "common/DB_Particles.hpp"
#include "common/FL_Event.hpp"
#include "common/R2DS_Cuts.hpp"
#include "common/VC_InjectedSexa.hpp"
#include "common/VC_McParticle.hpp"
#include "common/VC_McParticleView.hpp"
#include "common/VC_Track.hpp"
#include "common/VC_TrackView.hpp"
#include "common/VC_V0.hpp"

#include "App/Logger.hxx"
#include "KalmanFitter/KalmanFitterV0.hxx"
#if !R2DS_LEGACY_KF
#include "KalmanFitter/BaseKalmanFitter.hxx"
#include "Seeder/SeederHelixHelix.hxx"
#else
#include "Legacy/LegacyFitter.hxx"
#include "Legacy/LegacyHelixHelix.hxx"
#endif

namespace R2DS {

namespace KF = KalmanFitter;

bool Packager::Initialize() {

    // Prepare Input Reader //
    // PENDING: could refactor into a function

    fInput_Model = ROOT::RNTupleModel::Create();

    fInput_Event.AddFieldsTo(fInput_Model.get(), fSettings.IsMC);

    if (fSettings.IsMC) {
        fInput_InjectedSexa.AddFieldsTo(fInput_Model.get(), false);
        fInput_McParticle.AddFieldsTo(fInput_Model.get());
    }

    fInput_Track.AddFieldsTo(fInput_Model.get(), fSettings.IsMC, false, "Track");

    fReader =
        ROOT::RNTupleReader::Open(std::move(fInput_Model), E2R::Name_OutputRNT, fSettings.PathInputFiles.front());  // PENDING: handle multiple files?

    // Prepare Output File //

    fOutput_File = std::make_unique<TFile>(fSettings.PathOutputFile.c_str(), "RECREATE");
    if (!fOutput_File || fOutput_File->IsZombie()) {
        Logger::Error(__FUNCTION__, "Couldn't create TFile {}.", fSettings.PathOutputFile.c_str());
        return false;
    }

    PrepareOutputHistograms();

    // Prepare Output Writer //
    // PENDING: could refactor into a function

    fOutput_Model = ROOT::RNTupleModel::Create();

    fOutput_Event.AddFieldsTo(fOutput_Model.get(), fSettings.IsMC);

    fOutput_InjectedSexa.AddFieldsTo(fOutput_Model.get(), true);

    fOutput_NegKaon.AddFieldsTo(fOutput_Model.get(), fSettings.IsMC, fSettings.IsMC, "NK");
    fOutput_PosKaon.AddFieldsTo(fOutput_Model.get(), fSettings.IsMC, fSettings.IsMC, "PK");

    fOutput_AntiLambda.AddFieldsTo(fOutput_Model.get(), fSettings.IsMC, "AL");
    fOutput_Lambda.AddFieldsTo(fOutput_Model.get(), fSettings.IsMC, "L");
    fOutput_KaonZeroShort.AddFieldsTo(fOutput_Model.get(), fSettings.IsMC, "K0S");

    fWriter = ROOT::RNTupleWriter::Append(std::move(fOutput_Model), R2DS::Name_PackedRNT, *fOutput_File);

    Logger::Info(__FUNCTION__, "Packager initialized successfully.");

    return true;
}

// ## OUTPUT ZONE ## //

void Packager::PrepareOutputHistograms() {

    // event counter //
    fHist_EventCounter = std::make_unique<TH1D>("N_Events", ";;N_Events", 1, 0, 1);

    constexpr int x_nbins = 20;
    constexpr float x_min = 0.;
    constexpr float x_max = 20.;
    constexpr const char* hist_title = ";Cut N;N Passed Cut";

    fHist_CutFlow_AntiProton = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", DB::Particles::Particle("AntiProton").acronym).c_str(), hist_title, x_nbins, x_min, x_max);
    fHist_CutFlow_Proton = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", DB::Particles::Particle("Proton").acronym).c_str(), hist_title, x_nbins, x_min, x_max);
    fHist_CutFlow_NegKaon = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", DB::Particles::Particle("NegKaon").acronym).c_str(), hist_title, x_nbins, x_min, x_max);
    fHist_CutFlow_PosKaon = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", DB::Particles::Particle("PosKaon").acronym).c_str(), hist_title, x_nbins, x_min, x_max);
    fHist_CutFlow_PiMinus = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", DB::Particles::Particle("PiMinus").acronym).c_str(), hist_title, x_nbins, x_min, x_max);
    fHist_CutFlow_PiPlus = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", DB::Particles::Particle("PiPlus").acronym).c_str(), hist_title, x_nbins, x_min, x_max);

    fHist_CutFlow_AntiLambda = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", DB::Particles::Particle("AntiLambda").acronym).c_str(), hist_title, x_nbins, x_min, x_max);
    fHist_CutFlow_Lambda = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", DB::Particles::Particle("Lambda").acronym).c_str(), hist_title, x_nbins, x_min, x_max);
    fHist_CutFlow_KaonZeroShort = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", DB::Particles::Particle("KaonZeroShort").acronym).c_str(), hist_title, x_nbins, x_min, x_max);
}

// ## Event ZONE ## //

void Packager::ProcessEvent() {
    fOutput_Event = fInput_Event;
    fHist_EventCounter->Fill(0.);
    // cache pv
    fPrimaryVertex.SetCoordinates(*fInput_Event.PV_X, *fInput_Event.PV_Y, *fInput_Event.PV_Z);
}

// ## MC/Injected ZONE ## //

// Loop over all MC particles.
// Select Primary Particles (particles with no mother), generated via the Sexaquark-Reaction Generator, and with valid Reaction ID;
// and store their origin vertex as the coordinates for this particular secondary vertex.
void Packager::ProcessInjected() {

    // clang-format off
    std::array<float, E2R::NReactionsPerEvent> sv_x; sv_x.fill(Common::DummyFloat);
    std::array<float, E2R::NReactionsPerEvent> sv_y; sv_y.fill(Common::DummyFloat);
    std::array<float, E2R::NReactionsPerEvent> sv_z; sv_z.fill(Common::DummyFloat);
    std::array<char, E2R::NReactionsPerEvent> sv_found{};
    // clang-format on

    Vector::McParticleView mc_view(&fInput_McParticle);
    for (std::size_t entry_mc = 0; entry_mc < mc_view.Size(); ++entry_mc) {
        mc_view.Entry = entry_mc;  // NOTE: avoid calling the full `CacheCalculations()` machinery
        mc_view.CacheAscendantsInfo();
        mc_view.Classify_AsInSexaquarkSimulations(fSettings.ReactionChannel, Common::DummyInt);  // NOTE: no hypothesis

        if (!mc_view.IsFirstGenSignal()) continue;

        auto reaction_idx = static_cast<std::size_t>(mc_view.ReactionID() - E2R::ReactionID_Offset);
        if (sv_found[reaction_idx] == 1) continue;

        sv_x[reaction_idx] = mc_view.Origin_X();
        sv_y[reaction_idx] = mc_view.Origin_Y();
        sv_z[reaction_idx] = mc_view.Origin_Z();
        sv_found[reaction_idx] = 1;
    }

    // store //

    fOutput_InjectedSexa = fInput_InjectedSexa;
    for (std::size_t i = 0; i < E2R::NReactionsPerEvent; ++i) {
        fOutput_InjectedSexa.SV_X->push_back(sv_x[i]);
        fOutput_InjectedSexa.SV_Y->push_back(sv_y[i]);
        fOutput_InjectedSexa.SV_Z->push_back(sv_z[i]);
    }
}

// ## Tracks ZONE ## //

// Filter and group tracks into indices vectors, according to their respective species.
void Packager::ProcessTracks() {
    Vector::TrackView track_view(&fInput_Track);
    for (std::size_t entry_track = 0; entry_track < track_view.Size(); ++entry_track) {
        track_view.CacheCalculations(entry_track, fPrimaryVertex, *fInput_Event.MagneticField);
        // PID and pre-selection //
        if (track_view.Charge() < 0) {
            if (PassesProtonCuts(track_view, fHist_CutFlow_AntiProton.get())) {
                fEntries_AntiProton.push_back(entry_track);
            }
            if (PassesKaonCuts(track_view, fHist_CutFlow_NegKaon.get())) {
                fEntries_NegKaon.push_back(entry_track);
            }
            if (PassesPionCuts(track_view, fHist_CutFlow_PiMinus.get())) {
                fEntries_PiMinus.push_back(entry_track);
            }
        }
        if (track_view.Charge() > 0) {
            if (PassesProtonCuts(track_view, fHist_CutFlow_Proton.get())) {
                fEntries_Proton.push_back(entry_track);
            }
            if (PassesKaonCuts(track_view, fHist_CutFlow_PosKaon.get())) {
                fEntries_PosKaon.push_back(entry_track);
            }
            if (PassesPionCuts(track_view, fHist_CutFlow_PiPlus.get())) {
                fEntries_PiPlus.push_back(entry_track);
            }
        }
    }  // end of loop over tracks

#if R2DS_DEBUG
    Logger::Debug(__FUNCTION__, "n_antiprotons = {}", fEntries_AntiProton.size());
    Logger::Debug(__FUNCTION__, "n_protons     = {}", fEntries_Proton.size());
    Logger::Debug(__FUNCTION__, "n_negkaons    = {}", fEntries_NegKaon.size());
    Logger::Debug(__FUNCTION__, "n_poskaons    = {}", fEntries_PosKaon.size());
    Logger::Debug(__FUNCTION__, "n_piminus     = {}", fEntries_PiMinus.size());
    Logger::Debug(__FUNCTION__, "n_piplus      = {}", fEntries_PiPlus.size());
#endif
}

// NOTE: exclusive to secondary charged kaons.
void Packager::PackTracks(const DB::Particles::Definition& pid) {

    // determine rules based on particle species //
    std::vector<std::size_t>* vec_entries = nullptr;
    Vector::Track* out = nullptr;
    switch (pid.pdg_code) {
        case DB::Particles::Particle("NegKaon").pdg_code: {
            vec_entries = &fEntries_NegKaon;
            out = &fOutput_NegKaon;
            break;
        }
        case DB::Particles::Particle("PosKaon").pdg_code: {
            vec_entries = &fEntries_PosKaon;
            out = &fOutput_PosKaon;
            break;
        }
        default:
            return;
    }

    // prepare views //
    auto track = std::make_unique<Vector::TrackView>(&fInput_Track);
    std::unique_ptr<Vector::McParticleView> linked_mc = nullptr;
    if (fSettings.IsMC) linked_mc = std::make_unique<Vector::McParticleView>(&fInput_McParticle);

    // loop over selected tracks //
    for (const std::size_t& entry_track : *vec_entries) {
        track->CacheCalculations(entry_track, fPrimaryVertex, *fInput_Event.MagneticField);

        if (fSettings.IsMC) {
            auto mc_entry = track->McEntry();  // NOTE: cannot be invalid, by construction
            linked_mc->Entry = mc_entry;       // NOTE: prevent full `CacheCalculations()` machinery
            linked_mc->CacheAscendantsInfo();
            linked_mc->CacheDescendantsInfo();
            linked_mc->Classify_AsInSexaquarkSimulations(fSettings.ReactionChannel, pid.pdg_code);
        }

        Store(track.get(), linked_mc.get(), *out);
    }  // end of loop over selected tracks
}

void Packager::Store(const Vector::TrackView* track, const Vector::McParticleView* linked_mc, Vector::Track& df) {
    df.EsdEntry->push_back(track->EsdEntry());
    df.X->push_back(track->X());
    df.Y->push_back(track->Y());
    df.Z->push_back(track->Z());
    df.Px->push_back(track->Px());
    df.Py->push_back(track->Py());
    df.Pz->push_back(track->Pz());
    df.Charge->push_back(track->Charge());
    df.PreDCAxy->push_back(track->PreDCAxy());
    df.PreDCAz->push_back(track->PreDCAz());
    df.TPC_Signal->push_back(track->TPC_Signal());
    df.NSigmaPion->push_back(track->NSigmaPion());
    df.NSigmaKaon->push_back(track->NSigmaKaon());
    df.NSigmaProton->push_back(track->NSigmaProton());
    df.CovMatrix->push_back(track->CovMatrix());
    // PENDING E2R_EXTRA_INFO
    if (fSettings.IsMC) df.McEntry->push_back(track->McEntry());
    if (linked_mc == nullptr) {
        df.PdgCode->push_back(Common::DummyInt);
        df.Origin_X->push_back(Common::DummyFloat);
        df.Origin_Y->push_back(Common::DummyFloat);
        df.Origin_Z->push_back(Common::DummyFloat);
        df.MC_Px->push_back(Common::DummyFloat);
        df.MC_Py->push_back(Common::DummyFloat);
        df.MC_Pz->push_back(Common::DummyFloat);
        df.MC_Energy->push_back(Common::DummyFloat);
        df.IsTrue->push_back(Common::DummyChar);
        df.IsSecondary->push_back(Common::DummyChar);
        df.IsSignal->push_back(Common::DummyChar);
        df.ReactionID->push_back(Common::DummyInt);
        // mother
        df.Mother_McEntry->push_back(Common::DummyInt);
        df.Mother_PdgCode->push_back(Common::DummyInt);
        // grandmother
        df.GM_McEntry->push_back(Common::DummyInt);
        df.GM_PdgCode->push_back(Common::DummyInt);
    } else {
        df.PdgCode->push_back(linked_mc->PdgCode());
        df.Origin_X->push_back(linked_mc->Origin_X());
        df.Origin_Y->push_back(linked_mc->Origin_Y());
        df.Origin_Z->push_back(linked_mc->Origin_Z());
        df.MC_Px->push_back(linked_mc->Px());
        df.MC_Py->push_back(linked_mc->Py());
        df.MC_Pz->push_back(linked_mc->Pz());
        df.MC_Energy->push_back(linked_mc->Energy());
        df.IsTrue->push_back(static_cast<char>(linked_mc->IsTrue()));
        df.IsSecondary->push_back(static_cast<char>(linked_mc->IsSecondary()));
        df.IsSignal->push_back(static_cast<char>(linked_mc->IsSignal()));
        df.ReactionID->push_back(linked_mc->ReactionID());
        // mother
        df.Mother_McEntry->push_back(linked_mc->Mother_McEntry());
        df.Mother_PdgCode->push_back(linked_mc->Mother_PdgCode());
        // grandmother
        df.GM_McEntry->push_back(linked_mc->GrandMother_McEntry());
        df.GM_PdgCode->push_back(linked_mc->GrandMother_PdgCode());
    }
}

bool Packager::PassesProtonCuts(const Vector::TrackView& track, TH1D* cut_flow_hist) const {

    cut_flow_hist->Fill(0.);
    if (std::abs(track.NSigmaProton()) > R2DS::Cuts::Proton::AbsMax_NSigmaProton) return false;
    cut_flow_hist->Fill(1.);

    return true;
}

bool Packager::PassesKaonCuts(const Vector::TrackView& track, TH1D* cut_flow_hist) const {

    cut_flow_hist->Fill(0.);
    if (std::abs(track.NSigmaKaon()) > R2DS::Cuts::Kaon::AbsMax_NSigmaKaon) return false;
    cut_flow_hist->Fill(1.);

    return true;
}

bool Packager::PassesPionCuts(const Vector::TrackView& track, TH1D* cut_flow_hist) const {

    cut_flow_hist->Fill(0.);
    if (std::abs(track.NSigmaPion()) > R2DS::Cuts::Pion::AbsMax_NSigmaPion) return false;
    cut_flow_hist->Fill(1.);

    return true;
}

// ## V0s ZONE ## //

void Packager::FindV0s(const DB::Particles::Definition& pid) {

    // determine rules based on V0 species //
    Vector::V0* out = nullptr;
    const std::vector<std::size_t>* vec_neg = &fEntries_PiMinus;
    const std::vector<std::size_t>* vec_pos = &fEntries_PiPlus;
    switch (pid.pdg_code) {
        case DB::Particles::Particle("AntiLambda").pdg_code: {
            out = &fOutput_AntiLambda;
            vec_neg = &fEntries_AntiProton;
            break;
        }
        case DB::Particles::Particle("Lambda").pdg_code: {
            out = &fOutput_Lambda;
            vec_pos = &fEntries_Proton;
            break;
        }
        case DB::Particles::Particle("KaonZeroShort").pdg_code: {
            out = &fOutput_KaonZeroShort;
            break;
        }
        default: {
            Logger::Error(__FUNCTION__, "Invalid PID {} for a V0.", pid.name);
            return;
        }
    }
    auto pid_neg = DB::Particles::FindParticle(pid.daughters_pdg[0]);  // PENDING: needs helper, not well done...
    auto pid_pos = DB::Particles::FindParticle(pid.daughters_pdg[1]);

    // prepare views //
    auto neg = std::make_unique<Vector::TrackView>(&fInput_Track);
    auto pos = std::make_unique<Vector::TrackView>(&fInput_Track);
    std::unique_ptr<Vector::McParticleView> mc_neg;
    std::unique_ptr<Vector::McParticleView> mc_pos;
    std::unique_ptr<Vector::McParticleView> mc_v0;
    if (fSettings.IsMC) {
        mc_neg = std::make_unique<Vector::McParticleView>(&fInput_McParticle);
        mc_pos = std::make_unique<Vector::McParticleView>(&fInput_McParticle);
        mc_v0 = std::make_unique<Vector::McParticleView>(&fInput_McParticle);
    }

    // loop over all possible pairs of tracks //
    for (const auto& entry_neg : *vec_neg) {
        neg->CacheCalculations(entry_neg, fPrimaryVertex, *fInput_Event.MagneticField);

        if (fSettings.IsMC) {
            auto mc_entry_neg = neg->McEntry();  // NOTE: cannot be invalid, by construction
            mc_neg->Entry = mc_entry_neg;        // NOTE: prevent full `CacheCalculations()` machinery
            mc_neg->CacheAscendantsInfo();
            mc_neg->CacheDescendantsInfo();
            mc_neg->Classify_AsInSexaquarkSimulations(fSettings.ReactionChannel, pid_neg->pdg_code);
        }

        for (const auto& entry_pos : *vec_pos) {
            if (entry_neg == entry_pos) continue;
            pos->CacheCalculations(entry_pos, fPrimaryVertex, *fInput_Event.MagneticField);

            if (fSettings.IsMC) {
                auto mc_entry_pos = pos->McEntry();  // NOTE: cannot be invalid, by construction
                mc_pos->Entry = mc_entry_pos;        // NOTE: prevent full `CacheCalculations()` machinery
                mc_pos->CacheAscendantsInfo();
                mc_pos->CacheDescendantsInfo();
                mc_pos->Classify_AsInSexaquarkSimulations(fSettings.ReactionChannel, pid_pos->pdg_code);
            }

#if R2DS_LEGACY_KF
            // PCAs //
            auto [seed_neg, seed_pos, deriv_neg, deriv_pos] = Legacy::HelixHelix::FullPCAs(neg, pos, mass_neg, mass_pos, fInput_Event.magnetic_field);

            // fit vertex //
            auto l_fit = Legacy::Fit(neg, pos, mass_neg, mass_pos, fInput_Event.magnetic_field);
            auto fit = KF::Particle::FromLegacy(l_fit);
#else
            // PCAs //
            auto [seed_neg, seed_pos, pca_cache] = Seeder::HelixHelix::FastCorrectPCAs(*neg, *pos, *fInput_Event.MagneticField);

            // apply cuts (1) //
            if (!FastCuts(seed_neg, seed_pos, pid)) continue;

            // PCAs derivatives //
            auto [deriv_neg, deriv_pos] = Seeder::HelixHelix::ComputeDerivatives(seed_neg, seed_pos, pca_cache);

            // fit vertex //
            auto fit =
                KF::FitVertex(*neg, *pos, pid_neg->mass, pid_pos->mass, {seed_neg, deriv_neg}, {seed_pos, deriv_pos}, *fInput_Event.MagneticField);

            // create composite particle //
            KF::V0 v0(fit, seed_neg.pca, seed_pos.pca, *neg, *pos);

            // apply cuts (2) //
            if (!SlowCuts(v0, pid)) continue;
#endif
            // store //
            Store(&v0, mc_neg.get(), mc_pos.get(), mc_v0.get(), *out);
        }  // end of loop over pos
    }  // end of loop over neg
}

bool Packager::FastCuts_Lambda(const Seeder::Seed& seed_neg, const Seeder::Seed& seed_pos, TH1D* cut_flow_hist) const {
    cut_flow_hist->Fill(0.);

    if ((seed_neg.pca.GetXYZ_AsROOT() - seed_pos.pca.GetXYZ_AsROOT()).Mag2() > Cuts::Lambda::Max_DCAbtwDau * Cuts::Lambda::Max_DCAbtwDau)
        return false;
    cut_flow_hist->Fill(1.);

    return true;
}

bool Packager::FastCuts_KaonZeroShort(const Seeder::Seed& seed_neg, const Seeder::Seed& seed_pos, TH1D* cut_flow_hist) const {
    cut_flow_hist->Fill(0.);

    if ((seed_neg.pca.GetXYZ_AsROOT() - seed_pos.pca.GetXYZ_AsROOT()).Mag2() >
        R2DS::Cuts::KaonZeroShort::Max_DCAbtwDau * R2DS::Cuts::KaonZeroShort::Max_DCAbtwDau)
        return false;
    cut_flow_hist->Fill(1.);

    return true;
}

bool Packager::SlowCuts_Lambda(const KF::V0& v0, TH1D* cut_flow_hist) const {

    double mass = v0.Mass().value_or(Common::DummyDouble);  // cached
    if (mass < Cuts::Lambda::Min_Mass || mass > Cuts::Lambda::Max_Mass) return false;
    cut_flow_hist->Fill(2.);

    if (v0.SquaredRadius2D() < Cuts::Lambda::Min_Radius2D * Cuts::Lambda::Min_Radius2D) return false;
    cut_flow_hist->Fill(3.);

    if (v0.SquaredDCA_Neg_V0() > Cuts::Lambda::Max_DCAnegV0 * Cuts::Lambda::Max_DCAnegV0) return false;
    cut_flow_hist->Fill(4.);

    if (v0.SquaredDCA_Pos_V0() > Cuts::Lambda::Max_DCAposV0 * Cuts::Lambda::Max_DCAposV0) return false;
    cut_flow_hist->Fill(5.);

    // if (v0.Pt() < Cuts::Lambda::Min_Pt) return false; // PENDING
    // cut_flow_hist->Fill(6.); // PENDING

    if (std::abs(v0.Rapidity()) > Cuts::Lambda::AbsMax_Rapidity) return false;
    cut_flow_hist->Fill(7.);

    if (v0.AbsArmQtOverAlpha() > Cuts::Lambda::AbsMax_ArmQtOverAlpha) return false;
    cut_flow_hist->Fill(8.);

    double cpa_wrt_pv = v0.CPA_Vertex(*fInput_Event.PV_X, *fInput_Event.PV_Y, *fInput_Event.PV_Z);  // cached
    if (cpa_wrt_pv < Cuts::Lambda::Min_CPAwrtPV || cpa_wrt_pv > Cuts::Lambda::Max_CPAwrtPV) return false;
    cut_flow_hist->Fill(9.);

    if (v0.DCA_Vertex(*fInput_Event.PV_X, *fInput_Event.PV_Y, *fInput_Event.PV_Z) < Cuts::Lambda::Min_DCAwrtPV) return false;
    cut_flow_hist->Fill(10.);

    return true;
}

bool Packager::SlowCuts_KaonZeroShort(const KF::V0& v0, TH1D* cut_flow_hist) const {

    // if (v0.Pt() < Cuts::KaonZeroShort::Min_Pt) return false; // PENDING
    // cut_flow_hist->Fill(2.); // PENDING

    double mass = v0.Mass().value_or(Common::DummyDouble);  // cached
    if (mass < Cuts::KaonZeroShort::Min_Mass || mass > Cuts::KaonZeroShort::Max_Mass) return false;
    cut_flow_hist->Fill(3.);

    if (std::abs(v0.Rapidity()) > Cuts::KaonZeroShort::AbsMax_Rapidity) return false;
    cut_flow_hist->Fill(4.);

    if (v0.SquaredRadius2D() < Cuts::KaonZeroShort::Min_Radius2D * Cuts::KaonZeroShort::Min_Radius2D) return false;
    cut_flow_hist->Fill(5.);

    if (v0.SquaredDCA_Neg_V0() > Cuts::KaonZeroShort::Max_DCAnegV0 * Cuts::KaonZeroShort::Max_DCAnegV0) return false;
    cut_flow_hist->Fill(6.);

    if (v0.SquaredDCA_Pos_V0() > Cuts::KaonZeroShort::Max_DCAposV0 * Cuts::KaonZeroShort::Max_DCAposV0) return false;
    cut_flow_hist->Fill(7.);

    double cpa_wrt_pv = v0.CPA_Vertex(*fInput_Event.PV_X, *fInput_Event.PV_Y, *fInput_Event.PV_Z);  // cached
    if (cpa_wrt_pv < Cuts::KaonZeroShort::Min_CPAwrtPV || cpa_wrt_pv > Cuts::KaonZeroShort::Max_CPAwrtPV) return false;
    cut_flow_hist->Fill(8.);

    if (v0.DCA_Vertex(*fInput_Event.PV_X, *fInput_Event.PV_Y, *fInput_Event.PV_Z) < Cuts::KaonZeroShort::Min_DCAwrtPV) return false;
    cut_flow_hist->Fill(9.);

    return true;
}

void Packager::Store(const KF::V0* v0, const Vector::McParticleView* mc_neg, const Vector::McParticleView* mc_pos,
                     const Vector::McParticleView* mc_v0, Vector::V0& df) {
    df.Decay_X->push_back(static_cast<float>(v0->X()));
    df.Decay_Y->push_back(static_cast<float>(v0->Y()));
    df.Decay_Z->push_back(static_cast<float>(v0->Z()));
    df.Px->push_back(static_cast<float>(v0->Px()));
    df.Py->push_back(static_cast<float>(v0->Py()));
    df.Pz->push_back(static_cast<float>(v0->Pz()));
    df.Energy->push_back(static_cast<float>(v0->E()));
    // df.CovMatrix.Push(v0->Cov<7>()); // PENDING
    df.Chi2NDF->push_back(static_cast<float>(v0->Chi2NDF()));
    // -- neg daughter
    Store(&v0->Neg, mc_neg, df.Neg);
    df.Neg_PCAwrtV0_X->push_back(static_cast<float>(v0->Neg_PCAwrtV0.X()));
    df.Neg_PCAwrtV0_Y->push_back(static_cast<float>(v0->Neg_PCAwrtV0.Y()));
    df.Neg_PCAwrtV0_Z->push_back(static_cast<float>(v0->Neg_PCAwrtV0.Z()));
    df.Neg_PCAwrtV0_Px->push_back(static_cast<float>(v0->Neg_PCAwrtV0.Px()));
    df.Neg_PCAwrtV0_Py->push_back(static_cast<float>(v0->Neg_PCAwrtV0.Py()));
    df.Neg_PCAwrtV0_Pz->push_back(static_cast<float>(v0->Neg_PCAwrtV0.Pz()));
    // -- pos daughter
    Store(&v0->Pos, mc_pos, df.Pos);
    df.Pos_PCAwrtV0_X->push_back(static_cast<float>(v0->Pos_PCAwrtV0.X()));
    df.Pos_PCAwrtV0_Y->push_back(static_cast<float>(v0->Pos_PCAwrtV0.Y()));
    df.Pos_PCAwrtV0_Z->push_back(static_cast<float>(v0->Pos_PCAwrtV0.Z()));
    df.Pos_PCAwrtV0_Px->push_back(static_cast<float>(v0->Pos_PCAwrtV0.Px()));
    df.Pos_PCAwrtV0_Py->push_back(static_cast<float>(v0->Pos_PCAwrtV0.Py()));
    df.Pos_PCAwrtV0_Pz->push_back(static_cast<float>(v0->Pos_PCAwrtV0.Pz()));
    // -- linked mc info
    if (mc_v0 == nullptr) {
        df.McEntry->push_back(Common::DummyInt);
        df.PdgCode->push_back(Common::DummyInt);
        df.MC_Px->push_back(Common::DummyFloat);
        df.MC_Py->push_back(Common::DummyFloat);
        df.MC_Pz->push_back(Common::DummyFloat);
        df.MC_Energy->push_back(Common::DummyFloat);
        df.Origin_X->push_back(Common::DummyFloat);
        df.Origin_Y->push_back(Common::DummyFloat);
        df.Origin_Z->push_back(Common::DummyFloat);
        df.MC_Decay_X->push_back(Common::DummyFloat);
        df.MC_Decay_Y->push_back(Common::DummyFloat);
        df.MC_Decay_Z->push_back(Common::DummyFloat);
        df.IsTrue->push_back(Common::DummyChar);
        df.IsSignal->push_back(Common::DummyChar);
        df.IsSecondary->push_back(Common::DummyChar);
        df.IsHybrid->push_back(Common::DummyChar);
        df.ReactionID->push_back(Common::DummyInt);
        df.Mother_McEntry->push_back(Common::DummyInt);
        df.Mother_PdgCode->push_back(Common::DummyInt);
    } else {
        df.McEntry->push_back(static_cast<int>(mc_v0->Entry));
        df.PdgCode->push_back(mc_v0->PdgCode());
        df.MC_Px->push_back(mc_v0->Px());
        df.MC_Py->push_back(mc_v0->Py());
        df.MC_Pz->push_back(mc_v0->Pz());
        df.MC_Energy->push_back(mc_v0->Energy());
        df.Origin_X->push_back(mc_v0->Origin_X());
        df.Origin_Y->push_back(mc_v0->Origin_Y());
        df.Origin_Z->push_back(mc_v0->Origin_Z());
        df.MC_Decay_X->push_back(static_cast<float>(mc_v0->Decay_X()));
        df.MC_Decay_Y->push_back(static_cast<float>(mc_v0->Decay_Y()));
        df.MC_Decay_Z->push_back(static_cast<float>(mc_v0->Decay_Z()));
        df.IsTrue->push_back(static_cast<char>(mc_v0->IsTrue()));
        df.IsSignal->push_back(static_cast<char>(mc_v0->IsSignal()));
        df.IsSecondary->push_back(static_cast<char>(mc_v0->IsSecondary()));
        // df.IsHybrid->push_back(static_cast<char>(mc_v0->IsHybrid())); // PENDING: not so easy, huh
        df.ReactionID->push_back(mc_v0->ReactionID());
        df.Mother_McEntry->push_back(mc_v0->Mother_McEntry());
        df.Mother_PdgCode->push_back(mc_v0->Mother_PdgCode());
    }
}

// ## END OF CYCLES ## //

void Packager::EndOfEvent() {

    // fill RNTuple //

    fWriter->Fill();

    // clear temporary entries containers //

    fEntries_AntiProton.clear();
    fEntries_Proton.clear();
    fEntries_NegKaon.clear();
    fEntries_PosKaon.clear();
    fEntries_PiMinus.clear();
    fEntries_PiPlus.clear();

    // clear output vector branches //

    if (fSettings.IsMC) fOutput_InjectedSexa.Clear(true);

    fOutput_AntiLambda.Clear(fSettings.IsMC);
    fOutput_Lambda.Clear(fSettings.IsMC);
    fOutput_KaonZeroShort.Clear(fSettings.IsMC);

    fOutput_NegKaon.Clear(fSettings.IsMC, fSettings.IsMC);
    fOutput_PosKaon.Clear(fSettings.IsMC, fSettings.IsMC);
}

void Packager::EndOfAnalysis() {

    Logger::Info(__FUNCTION__, "The following objects have been written into TFile {}:", fSettings.PathOutputFile);

    Logger::Info(__FUNCTION__, "- RNTuple \"{}\"", R2DS::Name_PackedRNT);

    // write histograms //

    // -- event counter
    fHist_EventCounter->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_EventCounter->GetName());

    // -- selected tracks
    fHist_CutFlow_AntiProton->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow_AntiProton->GetName());
    fHist_CutFlow_Proton->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow_Proton->GetName());
    fHist_CutFlow_NegKaon->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow_NegKaon->GetName());
    fHist_CutFlow_PosKaon->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow_PosKaon->GetName());
    fHist_CutFlow_PiMinus->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow_PiMinus->GetName());
    fHist_CutFlow_PiPlus->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow_PiPlus->GetName());

    // -- found v0s
    fHist_CutFlow_AntiLambda->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow_AntiLambda->GetName());
    fHist_CutFlow_Lambda->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow_Lambda->GetName());
    fHist_CutFlow_KaonZeroShort->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow_KaonZeroShort->GetName());

    Logger::Info(__FUNCTION__, "All done.");
}

}  // namespace R2DS
