#include <cstddef>
#include <memory>

#include "common/Cached_V0.hpp"
#include "common/Constants.hpp"
#include "common/DB_Particles.hpp"
#include "common/MC_Helpers.hpp"
#include "common/Math.hpp"
#include "common/POD_Event.hpp"
#include "common/POD_InjectedSexa.hpp"
#include "common/POD_McParticle.hpp"
#include "common/POD_Track.hpp"
#include "common/POD_V0.hpp"
#include "common/T2DS_Cuts.hpp"

#include "App/Logger.hxx"
#include "KalmanFitter/KalmanFitterParticle.hxx"
#if !T2DS_LEGACY_KF
#include "KalmanFitter/BaseKalmanFitter.hxx"
#include "Seeder/BaseSeeder.hxx"
#include "Seeder/SeederHelixHelix.hxx"
#else
#include "Legacy/LegacyFitter.hxx"
#include "Legacy/LegacyHelixHelix.hxx"
#endif

#include "Packager/Packager.hxx"

namespace CMath = Common::Math;

namespace T2DS {

// ## OUTPUT ZONE ## //

void Packager::PrepareOutputHistograms() {

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
    // -- copy event info
    fOutput.Event = static_cast<POD::Event&>(fInput.Event);
    if (fSettings.IsMC) {
        fOutput.MC_Event = fInput.MC_Event;
    }
    // -- update event counter
    fHist_EventCounter->Fill(0.);
    // -- cache pv
    fPrimaryVertex.SetCoordinates(fOutput.Event.PV_X, fOutput.Event.PV_Y, fOutput.Event.PV_Z);
}

// ## MC/Injected ZONE ## //

// Loop over all MC particles.
// Select particles with no mother, generated via the AntiSexaquark-Reaction Generator, and with valid Reaction IDs;
// and store their origin vertex as the coordinates for this particular secondary vertex.
void Packager::ProcessInjected() {
    const double n_mass = DB::Particles::FindParticle(fReactionChannel.nucleon_pdg).mass;

    // clang-format off
    std::array<float, E2T::NSexaReactionsPerEvent> sv_x; sv_x.fill(Common::DummyFloat);
    std::array<float, E2T::NSexaReactionsPerEvent> sv_y; sv_y.fill(Common::DummyFloat);
    std::array<float, E2T::NSexaReactionsPerEvent> sv_z; sv_z.fill(Common::DummyFloat);
    std::array<bool, E2T::NSexaReactionsPerEvent> sv_found{};
    // clang-format on

    // -- loop over all mc particles
    for (const auto& mc : fInput.McParticle) {
        // -- select only first-gen signal products
        if (!MC::SexaquarkRules::IsGen1Signal(mc, fReactionChannel)) continue;
        // -- derive entry
        auto entry_injected = mc.StatusCode - E2T::ReactionID_Offset;
        // -- if sv not filled, fill it
        if (sv_found[entry_injected]) continue;
        sv_x[entry_injected] = mc.Origin_X;
        sv_y[entry_injected] = mc.Origin_Y;
        sv_z[entry_injected] = mc.Origin_Z;
        sv_found[entry_injected] = true;
    }

    // copy input injected sexa info //
    fOutput.InjectedSexa.resize(E2T::NSexaReactionsPerEvent);
    for (std::size_t entry_inj = 0; entry_inj < E2T::NSexaReactionsPerEvent; ++entry_inj) {
        // cache index lookup //
        POD::InjectedSexa& input_inj = fInput.InjectedSexa[entry_inj];
        POD::Extended::InjectedSexa& output_inj = fOutput.InjectedSexa[entry_inj];
        // fill values //
        static_cast<POD::InjectedSexa&>(output_inj) = input_inj;
        output_inj.SV_X = sv_x[entry_inj];
        output_inj.SV_Y = sv_y[entry_inj];
        output_inj.SV_Z = sv_z[entry_inj];
        output_inj.Energy = static_cast<float>(CMath::Hypot4(input_inj.Px, input_inj.Py, input_inj.Pz, fSettings.SexaquarkMass));
        output_inj.Nucleon_Energy = static_cast<float>(CMath::Hypot4(input_inj.Nucleon_Px, input_inj.Nucleon_Py, input_inj.Nucleon_Pz, n_mass));
    }
}

// ## Tracks ZONE ## //

// Filter and group tracks into indices vectors, according to their respective species.
void Packager::ProcessTracks() {

    // loop over all pre-selected tracks //
    for (std::size_t entry_track = 0; entry_track < fInput.Track.size(); ++entry_track) {
        const POD::Track& track = fInput.Track[entry_track];  // cache index lookup

        /* PENDING: cache calculations to speed up cuts! maybe not needed? */

        // PID and pre-selection //
        if (track.Charge < 0) {
            if (PassesProtonCuts(track, fHist_CutFlow_AntiProton.get())) {
                fEntries_AntiProton.push_back(entry_track);
            }
            if (PassesKaonCuts(track, fHist_CutFlow_NegKaon.get())) {
                fEntries_NegKaon.push_back(entry_track);
            }
            if (PassesPionCuts(track, fHist_CutFlow_PiMinus.get())) {
                fEntries_PiMinus.push_back(entry_track);
            }
        }
        if (track.Charge > 0) {
            if (PassesProtonCuts(track, fHist_CutFlow_Proton.get())) {
                fEntries_Proton.push_back(entry_track);
            }
            if (PassesKaonCuts(track, fHist_CutFlow_PosKaon.get())) {
                fEntries_PosKaon.push_back(entry_track);
            }
            if (PassesPionCuts(track, fHist_CutFlow_PiPlus.get())) {
                fEntries_PiPlus.push_back(entry_track);
            }
        }
    }  // end of loop over tracks

#if T2DS_DEBUG
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

    // determine rules and aliases //
    std::vector<std::size_t>* vec_entries = nullptr;
    const auto& input_tracks = fInput.Track;
    std::vector<POD::Track>* output_tracks = nullptr;
    std::vector<POD::Extended::McParticle>* output_mc = nullptr;
    switch (pid.pdg_code) {
        case DB::Particles::Particle("NegKaon").pdg_code: {
            vec_entries = &fEntries_NegKaon;
            output_tracks = &fOutput.NegKaon;
            output_mc = &fOutput.MC_NegKaon;
            break;
        }
        case DB::Particles::Particle("PosKaon").pdg_code: {
            vec_entries = &fEntries_PosKaon;
            output_tracks = &fOutput.PosKaon;
            output_mc = &fOutput.MC_PosKaon;
            break;
        }
        default:
            return;
    }

    // vector preallocation //
    output_tracks->reserve(vec_entries->size());
    if (fSettings.IsMC) output_mc->reserve(vec_entries->size());

    // loop over selected tracks //
    for (const std::size_t& entry_track : *vec_entries) {
        // NOTE: cuts were already applied in `ProcessTracks(...)`
        // -- reconstructed
        output_tracks->emplace_back(input_tracks[entry_track]);
        // -- mc
        if (fSettings.IsMC) {
            output_mc->emplace_back(BuildMcTrack(fInput.Track_McEntry[entry_track], pid.pdg_code, true));
        }
    }  // end of loop over selected tracks
}

bool Packager::PassesProtonCuts(const POD::Track& track, TH1D* cut_flow_hist) const {
    cut_flow_hist->Fill(0.);

    if (std::abs(track.NSigmasProton) > T2DS::Cuts::Proton::AbsMax_NSigmasProton) return false;
    cut_flow_hist->Fill(1.);

    return true;
}

bool Packager::PassesKaonCuts(const POD::Track& track, TH1D* cut_flow_hist) const {
    cut_flow_hist->Fill(0.);

    if (std::abs(track.NSigmasKaon) > T2DS::Cuts::Kaon::AbsMax_NSigmasKaon) return false;
    cut_flow_hist->Fill(1.);

    return true;
}

bool Packager::PassesPionCuts(const POD::Track& track, TH1D* cut_flow_hist) const {
    cut_flow_hist->Fill(0.);

    if (std::abs(track.NSigmasPion) > T2DS::Cuts::Pion::AbsMax_NSigmasPion) return false;
    cut_flow_hist->Fill(1.);

    return true;
}

POD::Extended::McParticle Packager::BuildMcTrack(unsigned int track_mc_entry, int pdg_code_hypothesis, bool include_gm) {
    // copy linked mc info //
    POD::Extended::McParticle new_mc(fInput.McParticle[track_mc_entry]);
    auto c = MC::SexaquarkRules::ClassifyDownstream(new_mc, fInput.McParticle, fReactionChannel, pdg_code_hypothesis, include_gm, false);
    MC::Apply(new_mc, c);
    return new_mc;
}

// ## V0s ZONE ## //

void Packager::FindV0s(const DB::Particles::Definition& pid) {

    // alias input //
    const auto& input_tracks = fInput.Track;

    // determine rules based on V0 species //
    const std::vector<std::size_t>* vec_entries_neg = &fEntries_PiMinus;
    const std::vector<std::size_t>* vec_entries_pos = &fEntries_PiPlus;
    auto pid_neg = DB::Particles::Particle("PiMinus");
    auto pid_pos = DB::Particles::Particle("PiPlus");
    std::vector<POD::V0>* output_vec_v0 = nullptr;
    std::vector<POD::Track>* output_vec_v0_neg = nullptr;
    std::vector<POD::Track>* output_vec_v0_pos = nullptr;
    std::vector<POD::Extended::McParticle>* output_vec_mc_v0 = nullptr;
    std::vector<POD::Extended::McParticle>* output_vec_mc_v0_neg = nullptr;
    std::vector<POD::Extended::McParticle>* output_vec_mc_v0_pos = nullptr;
    switch (pid.pdg_code) {
        case DB::Particles::Particle("AntiLambda").pdg_code: {
            vec_entries_neg = &fEntries_AntiProton;
            pid_neg = DB::Particles::Particle("AntiProton");
            output_vec_v0 = &fOutput.AntiLambda;
            output_vec_v0_neg = &fOutput.AntiLambda_Neg;
            output_vec_v0_pos = &fOutput.AntiLambda_Pos;
            output_vec_mc_v0 = &fOutput.MC_AntiLambda;
            output_vec_mc_v0_neg = &fOutput.MC_AntiLambda_Neg;
            output_vec_mc_v0_pos = &fOutput.MC_AntiLambda_Pos;
            break;
        }
        case DB::Particles::Particle("Lambda").pdg_code: {
            vec_entries_pos = &fEntries_Proton;
            pid_pos = DB::Particles::Particle("Proton");
            output_vec_v0 = &fOutput.Lambda;
            output_vec_v0_neg = &fOutput.Lambda_Neg;
            output_vec_v0_pos = &fOutput.Lambda_Pos;
            output_vec_mc_v0 = &fOutput.MC_Lambda;
            output_vec_mc_v0_neg = &fOutput.MC_Lambda_Neg;
            output_vec_mc_v0_pos = &fOutput.MC_Lambda_Pos;
            break;
        }
        case DB::Particles::Particle("KaonZeroShort").pdg_code: {
            output_vec_v0 = &fOutput.KaonZeroShort;
            output_vec_v0_neg = &fOutput.KaonZeroShort_Neg;
            output_vec_v0_pos = &fOutput.KaonZeroShort_Pos;
            output_vec_mc_v0 = &fOutput.MC_KaonZeroShort;
            output_vec_mc_v0_neg = &fOutput.MC_KaonZeroShort_Neg;
            output_vec_mc_v0_pos = &fOutput.MC_KaonZeroShort_Pos;
            break;
        }
        default: {
            Logger::Error(__FUNCTION__, "Invalid PID {} for a V0.", pid.name);
            return;
        }
    }

    // loop over all possible pairs of tracks //
    for (const auto& entry_neg : *vec_entries_neg) {
        const POD::Track& track_neg = input_tracks[entry_neg];  // cache index lookup

        for (const auto& entry_pos : *vec_entries_pos) {
            if (entry_neg == entry_pos) continue;
            const POD::Track& track_pos = input_tracks[entry_pos];  // cache index lookup

#if T2DS_LEGACY_KF
            // PCAs //
            auto [seed_neg, seed_pos, deriv_neg, deriv_pos] = Legacy::HelixHelix::FullPCAs(neg, pos, mass_neg, mass_pos, fInput_Event.magnetic_field);

            // fit vertex //
            auto l_fit = Legacy::Fit(neg, pos, mass_neg, mass_pos, fOutput.Event.magnetic_field);
            auto fit = KF::Particle::FromLegacy(l_fit);
#else

            // PCAs //
            auto [seed_neg, seed_pos, pca_cache] = Seeder::HelixHelix::FastCorrectPCAs(track_neg, track_pos, fOutput.Event.MagneticField);

            // apply cuts (1) //
            if (!FastCuts(seed_neg.pca, seed_pos.pca, pid)) continue;

            // PCAs derivatives //
            auto [deriv_neg, deriv_pos] = Seeder::HelixHelix::ComputeDerivatives(seed_neg, seed_pos, pca_cache);

            // fit vertex //
            auto fit = KF::FitVertex(track_neg, track_pos, pid_neg.mass, pid_pos.mass, {seed_neg, deriv_neg}, {seed_pos, deriv_pos},
                                     fOutput.Event.MagneticField);

            // create storage+computation units //
            POD::V0 v0 = CreateV0(fit, seed_neg.pca, seed_pos.pca);
            Cached::V0 c_v0(v0, track_neg, track_pos, fPrimaryVertex);

            // apply cuts (2) //
            if (!SlowCuts(c_v0, pid)) continue;
#endif

            // store reconstructed //
            output_vec_v0->emplace_back(v0);
            output_vec_v0_neg->emplace_back(track_neg);
            output_vec_v0_pos->emplace_back(track_pos);

            // store mc //
            if (fSettings.IsMC) {
                // -- neg
                POD::Extended::McParticle new_mc_neg = BuildMcTrack(fInput.Track_McEntry[entry_neg], pid_neg.pdg_code, true);
                output_vec_mc_v0_neg->emplace_back(new_mc_neg);
                // -- pos
                POD::Extended::McParticle new_mc_pos = BuildMcTrack(fInput.Track_McEntry[entry_pos], pid_pos.pdg_code, true);
                output_vec_mc_v0_pos->emplace_back(new_mc_pos);
                // -- v0
                output_vec_mc_v0->emplace_back(BuildMcV0(new_mc_neg, new_mc_pos, pid.pdg_code));
            }
            // -- push
        }  // end of loop over pos
    }  // end of loop over neg
}

bool Packager::FastCuts_Lambda(const Seeder::PCA& pca_neg, const Seeder::PCA& pca_pos, TH1D* cut_flow_hist) const {
    cut_flow_hist->Fill(0.);

    if (CMath::SquaredDistance(pca_neg.xyz, pca_pos.xyz) > Cuts::Lambda::Max_DCAbtwDau * Cuts::Lambda::Max_DCAbtwDau) {
        return false;
    }
    cut_flow_hist->Fill(1.);

    return true;
}

bool Packager::FastCuts_KaonZeroShort(const Seeder::PCA& pca_neg, const Seeder::PCA& pca_pos, TH1D* cut_flow_hist) const {
    cut_flow_hist->Fill(0.);

    if (CMath::SquaredDistance(pca_neg.xyz, pca_pos.xyz) > T2DS::Cuts::KaonZeroShort::Max_DCAbtwDau * T2DS::Cuts::KaonZeroShort::Max_DCAbtwDau) {
        return false;
    }
    cut_flow_hist->Fill(1.);

    return true;
}

bool Packager::SlowCuts_Lambda(const Cached::V0& c_v0, TH1D* cut_flow_hist) const {

    double mass = c_v0.Mass();  // cached
    if (mass < Cuts::Lambda::Min_Mass || mass > Cuts::Lambda::Max_Mass) return false;
    cut_flow_hist->Fill(2.);

    if (c_v0.Decay_SquaredRadius2D() < Cuts::Lambda::Min_Radius2D * Cuts::Lambda::Min_Radius2D) return false;
    cut_flow_hist->Fill(3.);

    if (c_v0.Neg_SquaredDCA_wrt_V0() > Cuts::Lambda::Max_DCAnegV0 * Cuts::Lambda::Max_DCAnegV0) return false;
    cut_flow_hist->Fill(4.);

    if (c_v0.Pos_SquaredDCA_wrt_V0() > Cuts::Lambda::Max_DCAposV0 * Cuts::Lambda::Max_DCAposV0) return false;
    cut_flow_hist->Fill(5.);

    // if (c_v0.Pt() < Cuts::Lambda::Min_Pt) return false; // PENDING
    // cut_flow_hist->Fill(6.); // PENDING

    if (std::abs(c_v0.Rapidity()) > Cuts::Lambda::AbsMax_Rapidity) return false;
    cut_flow_hist->Fill(7.);

    // if (c_v0.AbsArmQtOverAlpha() > Cuts::Lambda::AbsMax_ArmQtOverAlpha) return false;  // PENDING: not really sure if i like this cut, actually
    // cut_flow_hist->Fill(8.);

    if (c_v0.CPA_wrt_PV() < Cuts::Lambda::Min_CPAwrtPV || c_v0.CPA_wrt_PV() > Cuts::Lambda::Max_CPAwrtPV) return false;
    cut_flow_hist->Fill(8.);

    if (c_v0.SquaredDCA_wrt_PV() < Cuts::Lambda::Min_DCAwrtPV * Cuts::Lambda::Min_DCAwrtPV) return false;
    cut_flow_hist->Fill(9.);

    return true;
}

bool Packager::SlowCuts_KaonZeroShort(const Cached::V0& c_v0, TH1D* cut_flow_hist) const {

    // if (c_v0.Pt() < Cuts::KaonZeroShort::Min_Pt) return false; // PENDING
    // cut_flow_hist->Fill(2.); // PENDING

    double mass = c_v0.Mass();  // cached
    if (mass < Cuts::KaonZeroShort::Min_Mass || mass > Cuts::KaonZeroShort::Max_Mass) return false;
    cut_flow_hist->Fill(3.);

    if (std::abs(c_v0.Rapidity()) > Cuts::KaonZeroShort::AbsMax_Rapidity) return false;
    cut_flow_hist->Fill(4.);

    if (c_v0.Decay_SquaredRadius2D() < Cuts::KaonZeroShort::Min_Radius2D * Cuts::KaonZeroShort::Min_Radius2D) return false;
    cut_flow_hist->Fill(5.);

    if (c_v0.Neg_SquaredDCA_wrt_V0() > Cuts::KaonZeroShort::Max_DCAnegV0 * Cuts::KaonZeroShort::Max_DCAnegV0) return false;
    cut_flow_hist->Fill(6.);

    if (c_v0.Pos_SquaredDCA_wrt_V0() > Cuts::KaonZeroShort::Max_DCAposV0 * Cuts::KaonZeroShort::Max_DCAposV0) return false;
    cut_flow_hist->Fill(7.);

    if (c_v0.CPA_wrt_PV() < Cuts::KaonZeroShort::Min_CPAwrtPV || c_v0.CPA_wrt_PV() > Cuts::KaonZeroShort::Max_CPAwrtPV) return false;
    cut_flow_hist->Fill(8.);

    if (c_v0.SquaredDCA_wrt_PV() < Cuts::KaonZeroShort::Min_DCAwrtPV * Cuts::KaonZeroShort::Min_DCAwrtPV) return false;
    cut_flow_hist->Fill(9.);

    return true;
}

POD::Extended::McParticle Packager::BuildMcV0(const POD::Extended::McParticle& mc_neg, const POD::Extended::McParticle& mc_pos,
                                              int pdg_code_hypothesis) {
    POD::Extended::McParticle mc_v0;

    // -- fill hybridness, independetly of no common mother
    mc_v0.IsHybrid = mc_neg.IsTrueSignal != mc_pos.IsTrueSignal || mc_neg.SignalID != mc_pos.SignalID;

    auto mc_entry = MC::FindCommonMotherMcEntry(mc_neg, mc_pos);
    if (!mc_entry.has_value()) return mc_v0;

    // fill values //
    static_cast<POD::McParticle&>(mc_v0) = fInput.McParticle[mc_entry.value()];
    MC::Apply(mc_v0, MC::SexaquarkRules::ClassifyDownstream(mc_v0, fInput.McParticle, fReactionChannel, pdg_code_hypothesis, false, true));

    return mc_v0;
}

POD::V0 Packager::CreateV0(const KF::Particle& fit, const Seeder::PCA& neg_pca_wrt_v0, const Seeder::PCA& pos_pca_wrt_v0) {
    POD::V0 new_v0;  // non-initialized on purpose
    new_v0.Decay_X = static_cast<float>(fit.X());
    new_v0.Decay_Y = static_cast<float>(fit.Y());
    new_v0.Decay_Z = static_cast<float>(fit.Z());
    new_v0.Px = static_cast<float>(fit.Px());
    new_v0.Py = static_cast<float>(fit.Py());
    new_v0.Pz = static_cast<float>(fit.Pz());
    new_v0.Energy = static_cast<float>(fit.E());
    new_v0.CovMatrix = fit.Cov<7>();
    new_v0.Chi2NDF = static_cast<float>(fit.Chi2NDF());
    // -- negative daughter
    new_v0.Neg_PCAwrtV0_X = static_cast<float>(neg_pca_wrt_v0.X());
    new_v0.Neg_PCAwrtV0_Y = static_cast<float>(neg_pca_wrt_v0.Y());
    new_v0.Neg_PCAwrtV0_Z = static_cast<float>(neg_pca_wrt_v0.Z());
    new_v0.Neg_PCAwrtV0_Px = static_cast<float>(neg_pca_wrt_v0.Px());
    new_v0.Neg_PCAwrtV0_Py = static_cast<float>(neg_pca_wrt_v0.Py());
    new_v0.Neg_PCAwrtV0_Pz = static_cast<float>(neg_pca_wrt_v0.Pz());
    // -- positive daughter
    new_v0.Pos_PCAwrtV0_X = static_cast<float>(pos_pca_wrt_v0.X());
    new_v0.Pos_PCAwrtV0_Y = static_cast<float>(pos_pca_wrt_v0.Y());
    new_v0.Pos_PCAwrtV0_Z = static_cast<float>(pos_pca_wrt_v0.Z());
    new_v0.Pos_PCAwrtV0_Px = static_cast<float>(pos_pca_wrt_v0.Px());
    new_v0.Pos_PCAwrtV0_Py = static_cast<float>(pos_pca_wrt_v0.Py());
    new_v0.Pos_PCAwrtV0_Pz = static_cast<float>(pos_pca_wrt_v0.Pz());

    return new_v0;
}

// ## END OF CYCLES ## //

void Packager::EndOfEvent() {
    // fill schema
    fWriter->Fill();

    // clear schema vectors
    fOutput.Clear(fSettings.IsMC);

    // clear temporary entries containers
    fEntries_AntiProton.clear();
    fEntries_Proton.clear();
    fEntries_NegKaon.clear();
    fEntries_PosKaon.clear();
    fEntries_PiMinus.clear();
    fEntries_PiPlus.clear();
}

void Packager::EndOfAnalysis() {

    Logger::Info(__FUNCTION__, "The following objects have been written into TFile \"{}\":", fSettings.PathOutputFile);

    Logger::Info(__FUNCTION__, "- RNTuple \"{}\"", T2DS::Name_PackedRNT);

    // write histograms //

    fOutput_File->cd();

    // -- event counter
    fHist_EventCounter->Write();
    Logger::Info(__FUNCTION__, "- TH1D \"{}\"", fHist_EventCounter->GetName());

    // -- cut flow selected tracks
    fHist_CutFlow_AntiProton->Write();
    Logger::Info(__FUNCTION__, "- TH1D \"{}\"", fHist_CutFlow_AntiProton->GetName());
    fHist_CutFlow_Proton->Write();
    Logger::Info(__FUNCTION__, "- TH1D \"{}\"", fHist_CutFlow_Proton->GetName());
    fHist_CutFlow_NegKaon->Write();
    Logger::Info(__FUNCTION__, "- TH1D \"{}\"", fHist_CutFlow_NegKaon->GetName());
    fHist_CutFlow_PosKaon->Write();
    Logger::Info(__FUNCTION__, "- TH1D \"{}\"", fHist_CutFlow_PosKaon->GetName());
    fHist_CutFlow_PiMinus->Write();
    Logger::Info(__FUNCTION__, "- TH1D \"{}\"", fHist_CutFlow_PiMinus->GetName());
    fHist_CutFlow_PiPlus->Write();
    Logger::Info(__FUNCTION__, "- TH1D \"{}\"", fHist_CutFlow_PiPlus->GetName());

    // -- cut flow found v0s
    fHist_CutFlow_AntiLambda->Write();
    Logger::Info(__FUNCTION__, "- TH1D \"{}\"", fHist_CutFlow_AntiLambda->GetName());
    fHist_CutFlow_Lambda->Write();
    Logger::Info(__FUNCTION__, "- TH1D \"{}\"", fHist_CutFlow_Lambda->GetName());
    fHist_CutFlow_KaonZeroShort->Write();
    Logger::Info(__FUNCTION__, "- TH1D \"{}\"", fHist_CutFlow_KaonZeroShort->GetName());

    Logger::Info(__FUNCTION__, "All done.");
}

}  // namespace T2DS
