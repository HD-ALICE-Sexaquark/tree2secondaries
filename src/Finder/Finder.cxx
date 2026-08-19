#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <vector>

#include "common/Cached_Sexaquark.hpp"
#include "common/Cached_V0.hpp"
#include "common/Cuts_T2DS_Finder.hpp"
#include "common/DB_Particles.hpp"
#include "common/DB_ReactionChannels.hpp"
#include "common/MC_Helpers.hpp"
#include "common/POD_InjectedSexa.hpp"
#include "common/POD_McParticle.hpp"
#include "common/POD_Sexaquark.hpp"
#include "common/POD_Track.hpp"
#include "common/POD_V0.hpp"

#include "common/Math.hpp"
namespace CMath = Common::Math;

#include "KalmanFitter/BaseKalmanFitter.hxx"
#include "Seeder/BaseSeeder.hxx"
#include "Seeder/SeederHelixHelix.hxx"
#include "Seeder/SeederHelixLine.hxx"
#include "Seeder/SeederLineLine.hxx"

#include "Finder/Finder.hxx"

namespace T2DS {

// ## OUTPUT ZONE ## //

void Finder::PrepareOutputHistograms() {
    // event counter
    fHist_EventCounter = std::make_unique<TH1D>("N_Events", ";N_Events;", 1, 0., 1.);

    // cut flows
    constexpr const char* hist_title = ";;N Passed Cut";
    TAxis* x_axis = nullptr;

    // -- secondary (anti)protons
    fHist_CutFlow_AntiProton =
        std::make_unique<TH1D>(std::format("CutFlow_{}", DB::Particles::Particle("AntiProton").acronym).c_str(),  //
                               hist_title, static_cast<int>(EProton::kNProtonCuts), 0., static_cast<double>(EProton::kNProtonCuts));
    fHist_CutFlow_Proton =
        std::make_unique<TH1D>(std::format("CutFlow_{}", DB::Particles::Particle("Proton").acronym).c_str(),  //
                               hist_title, static_cast<int>(EProton::kNProtonCuts), 0., static_cast<double>(EProton::kNProtonCuts));
    //  > define bin labels
    for (auto* hist_proton : {fHist_CutFlow_AntiProton.get(), fHist_CutFlow_Proton.get()}) {
        x_axis = hist_proton->GetXaxis();
        x_axis->SetBinLabel(static_cast<int>(EProton::kAllPossibleProtons) + 1, "AllPossibleProtons");
        x_axis->SetBinLabel(static_cast<int>(EProton::kPasses_NSigmasProtons) + 1, "Passes_NSigmasProtons");
    }

    // -- secondary charged kaons
    fHist_CutFlow_NegKaon = std::make_unique<TH1D>(std::format("CutFlow_{}", DB::Particles::Particle("NegKaon").acronym).c_str(),  //
                                                   hist_title, static_cast<int>(EKaon::kNKaonCuts), 0., static_cast<double>(EKaon::kNKaonCuts));
    fHist_CutFlow_PosKaon = std::make_unique<TH1D>(std::format("CutFlow_{}", DB::Particles::Particle("PosKaon").acronym).c_str(),  //
                                                   hist_title, static_cast<int>(EKaon::kNKaonCuts), 0., static_cast<double>(EKaon::kNKaonCuts));
    //  > define bin labels
    for (auto* hist_kaon : {fHist_CutFlow_NegKaon.get(), fHist_CutFlow_PosKaon.get()}) {
        x_axis = hist_kaon->GetXaxis();
        x_axis->SetBinLabel(static_cast<int>(EKaon::kAllPossibleKaons) + 1, "AllPossibleKaons");
        x_axis->SetBinLabel(static_cast<int>(EKaon::kPasses_NSigmasKaons) + 1, "Passes_NSigmasKaons");
    }

    // -- secondary charged pions
    fHist_CutFlow_PiMinus = std::make_unique<TH1D>(std::format("CutFlow_{}", DB::Particles::Particle("PiMinus").acronym).c_str(),  //
                                                   hist_title, static_cast<int>(EPion::kNPionCuts), 0., static_cast<double>(EPion::kNPionCuts));
    fHist_CutFlow_PiPlus = std::make_unique<TH1D>(std::format("CutFlow_{}", DB::Particles::Particle("PiPlus").acronym).c_str(),  //
                                                  hist_title, static_cast<int>(EPion::kNPionCuts), 0., static_cast<double>(EPion::kNPionCuts));
    //  > define bin labels
    for (auto* hist_pion : {fHist_CutFlow_PiMinus.get(), fHist_CutFlow_PiPlus.get()}) {
        x_axis = hist_pion->GetXaxis();
        x_axis->SetBinLabel(static_cast<int>(EPion::kAllPossiblePions) + 1, "AllPossiblePions");
        x_axis->SetBinLabel(static_cast<int>(EPion::kPasses_NSigmasPions) + 1, "Passes_NSigmasPions");
    }

    // -- for secondary (anti)lambdas
    fHist_CutFlow_AntiLambda =
        std::make_unique<TH1D>(std::format("CutFlow_{}", DB::Particles::Particle("AntiLambda").acronym).c_str(),  //
                               hist_title, static_cast<int>(ELambda::kNLambdaCuts), 0., static_cast<double>(ELambda::kNLambdaCuts));
    fHist_CutFlow_Lambda =
        std::make_unique<TH1D>(std::format("CutFlow_{}", DB::Particles::Particle("Lambda").acronym).c_str(),  //
                               hist_title, static_cast<int>(ELambda::kNLambdaCuts), 0., static_cast<double>(ELambda::kNLambdaCuts));
    //  > define bin labels
    for (auto* hist_lambda : {fHist_CutFlow_AntiLambda.get(), fHist_CutFlow_Lambda.get()}) {
        x_axis = hist_lambda->GetXaxis();
        x_axis->SetBinLabel(static_cast<int>(ELambda::kAllCombinations) + 1, "AllCombinations");
        x_axis->SetBinLabel(static_cast<int>(ELambda::kPasses_DcaBtwDaughters) + 1, "Passes_DcaBtwDaughters");
        // PENDING
    }

    // -- for secondary k0s
    fHist_CutFlow_KaonZeroShort =
        std::make_unique<TH1D>(std::format("CutFlow_{}", DB::Particles::Particle("KaonZeroShort").acronym).c_str(), hist_title,
                               static_cast<int>(EKaonZeroShort::kNKaonZeroShortCuts), 0., static_cast<double>(EKaonZeroShort::kNKaonZeroShortCuts));
    //  > define bin labels
    x_axis = fHist_CutFlow_KaonZeroShort->GetXaxis();
    x_axis->SetBinLabel(static_cast<int>(EKaonZeroShort::kAllCombinations) + 1, "AllCombinations");
    x_axis->SetBinLabel(static_cast<int>(EKaonZeroShort::kPasses_DcaBtwDaughters) + 1, "Passes_DcaBtwDaughters");
    // PENDING

    // -- for channel a
    fHist_CutFlow_ChannelA = std::make_unique<TH1D>("CutFlow_ChannelA", hist_title, static_cast<int>(EChannelA::kNChannelACuts), 0.,
                                                    static_cast<double>(EChannelA::kNChannelACuts));
    fHist_CutFlow_ChannelA_Bkg = std::make_unique<TH1D>("CutFlow_ChannelA_Bkg", hist_title, static_cast<int>(EChannelA::kNChannelACuts), 0.,
                                                        static_cast<double>(EChannelA::kNChannelACuts));
    //  > define bin labels
    // for (auto* hist_sexa : {fHist_CutFlow_ChannelA.get(), fHist_CutFlow_ChannelA_Bkg.get()}) {
    // x_axis = hist_sexa->GetXaxis();
    // PENDING
    // }

    // -- for channel d
    fHist_CutFlow_ChannelD = std::make_unique<TH1D>("CutFlow_ChannelD", hist_title, static_cast<int>(EChannelD::kNChannelDCuts), 0.,
                                                    static_cast<double>(EChannelD::kNChannelDCuts));
    fHist_CutFlow_ChannelD_Bkg = std::make_unique<TH1D>("CutFlow_ChannelD_Bkg", hist_title, static_cast<int>(EChannelD::kNChannelDCuts), 0.,
                                                        static_cast<double>(EChannelD::kNChannelDCuts));
    //  > define bin labels
    // for (auto* hist_sexa : {fHist_CutFlow_ChannelD.get(), fHist_CutFlow_ChannelD_Bkg.get()}) {
    // x_axis = hist_sexa->GetXaxis();
    // PENDING
    // }

    // -- for channel h
    fHist_CutFlow_ChannelH = std::make_unique<TH1D>("CutFlow_ChannelH", hist_title, static_cast<int>(EChannelH::kNChannelHCuts), 0.,
                                                    static_cast<double>(EChannelH::kNChannelHCuts));
    fHist_CutFlow_ChannelH_Bkg = std::make_unique<TH1D>("CutFlow_ChannelH_Bkg", hist_title, static_cast<int>(EChannelH::kNChannelHCuts), 0.,
                                                        static_cast<double>(EChannelH::kNChannelHCuts));
    //  > define bin labels
    // for (auto* hist_sexa : {fHist_CutFlow_ChannelH.get(), fHist_CutFlow_ChannelH_Bkg.get()}) {
    // x_axis = hist_sexa->GetXaxis();
    // PENDING
    // }
}

// ## Event ZONE ## //

void Finder::ProcessEvent() {
    // update event counter
    fHist_EventCounter->Fill(0.);

    // cache pv and magnetic field
    fPrimaryVertex.SetCoordinates(fInput.Event.PV_X, fInput.Event.PV_Y, fInput.Event.PV_Z);
    fMagneticField = fInput.Event.MagneticField;

    // copy event info into every output rntuple
    fOutput.Event = fInput.Event;
    if (fSettings.IsMC) fOutput.MC_Event = fInput.MC_Event;
}

// ## Injected/MC ZONE ## //

// Loop over all MC particles.
// Select particles with no mother, generated via the AntiSexaquark-Reaction Generator, and with valid Reaction IDs;
// and store their origin vertex as the coordinates for this particular secondary vertex.
void Finder::ProcessInjected() {

    if (!fMcSignalChannel.has_value()) {
        fMcSignalChannel = MC::SexaquarkRules::DetectMcSignalChannel(fInput.McParticle);
        Logger::Info(__FUNCTION__, "Injected reaction channel identified as \"Channel {}\".", fMcSignalChannel->name);
    }

    // there's nothing to collect; copy what is already there, if anything
    if (fMcSignalChannel->name == '0') {
        for (const auto& input_inj : fInput.InjectedSexa) fOutput.Injected.emplace_back(POD::Extended::InjectedSexa{input_inj});
        return;
    }

    const std::size_t n_injected = fInput.InjectedSexa.size();
    if (n_injected == 0) return;

    // struck nucleon mass
    const double n_mass = DB::Particles::FindParticle(fMcSignalChannel->nucleon_pdg).mass;

    // accumulate the first-gen products of every reaction //
    std::vector<InjectedReaction> reactions(n_injected);
    for (const auto& mc : fInput.McParticle) {
        // -- select only first-gen signal products
        if (!MC::SexaquarkRules::IsGen1Signal(mc, fMcSignalChannel.value())) continue;
        // -- derive entry
        auto entry_injected = mc.StatusCode - E2T::ReactionID_Offset;  // `IsGen1Signal` guarantees status code to lie within [600,620[
        if (entry_injected >= n_injected) continue;
        InjectedReaction& reaction = reactions[entry_injected];
        // -- add up lorentz vector
        reaction.after_px += static_cast<double>(mc.Px);
        reaction.after_py += static_cast<double>(mc.Py);
        reaction.after_pz += static_cast<double>(mc.Pz);
        reaction.after_e += static_cast<double>(mc.Energy);
        // -- if sv not filled, fill it
        if (reaction.found == 1) continue;
        reaction.sv_x = mc.Origin_X;
        reaction.sv_y = mc.Origin_Y;
        reaction.sv_z = mc.Origin_Z;
        reaction.found = 1;
    }

    // copy input injected sexa info //
    fOutput.Injected.resize(n_injected);
    for (std::size_t entry_inj = 0; entry_inj < n_injected; ++entry_inj) {
        // cache index lookups //
        const POD::InjectedSexa& input_inj = fInput.InjectedSexa[entry_inj];
        const InjectedReaction& reaction = reactions[entry_inj];
        POD::Extended::InjectedSexa& output_inj = fOutput.Injected[entry_inj];
        // fill values //
        static_cast<POD::InjectedSexa&>(output_inj) = input_inj;
        output_inj.Energy = static_cast<float>(CMath::Hypot4(input_inj.Px, input_inj.Py, input_inj.Pz, fSettings.SexaquarkMass));
        output_inj.Nucleon_Energy = static_cast<float>(CMath::Hypot4(input_inj.Nucleon_Px, input_inj.Nucleon_Py, input_inj.Nucleon_Pz, n_mass));
        if (reaction.found == 0) continue;
        output_inj.SV_X = reaction.sv_x;
        output_inj.SV_Y = reaction.sv_y;
        output_inj.SV_Z = reaction.sv_z;
        output_inj.After_Px = static_cast<float>(reaction.after_px);
        output_inj.After_Py = static_cast<float>(reaction.after_py);
        output_inj.After_Pz = static_cast<float>(reaction.after_pz);
        output_inj.After_Energy = static_cast<float>(reaction.after_e);
    }
}

POD::Linked::InjectedSexa Finder::BuildMcSexaquark(const POD::Extended::McParticle& mc_dau1, const POD::Extended::McParticle& mc_dau2) {
    POD::Linked::InjectedSexa mc_sexa;

    // fill hybridness, independently of no common reaction id
    mc_sexa.IsHybrid = mc_dau1.IsTrueSignal != mc_dau2.IsTrueSignal || mc_dau1.SignalID != mc_dau2.SignalID || mc_dau1.IsHybrid || mc_dau2.IsHybrid;

    // find common reaction id
    auto entry_inj = MC::SexaquarkRules::FindCommonReactionID(mc_dau1, mc_dau2);
    if (!entry_inj.has_value() || entry_inj.value() >= fOutput.Injected.size()) {
        return mc_sexa;
    }

    // if found, fill it with matching injected info
    static_cast<POD::Extended::InjectedSexa&>(mc_sexa) = fOutput.Injected[entry_inj.value()];

    return mc_sexa;
}

// ## Tracks ZONE ## //

// Filter and group tracks into per-species vectors, together with their linked mc particles.
// NOTE: a track can enter more than one species, as the pid hypotheses aren't exclusive.
void Finder::ProcessTracks() {

    // vectors preallocation //
    const std::size_t n_total_tracks = fInput.Track.size();
    fTemp_AntiProton.reserve(n_total_tracks);
    fTemp_Proton.reserve(n_total_tracks);
    fTemp_NegKaon.reserve(n_total_tracks);
    fTemp_PosKaon.reserve(n_total_tracks);
    fTemp_PiMinus.reserve(n_total_tracks);
    fTemp_PiPlus.reserve(n_total_tracks);
    if (fSettings.IsMC) {
        fTemp_MC_AntiProton.reserve(n_total_tracks);
        fTemp_MC_Proton.reserve(n_total_tracks);
        fTemp_MC_NegKaon.reserve(n_total_tracks);
        fTemp_MC_PosKaon.reserve(n_total_tracks);
        fTemp_MC_PiMinus.reserve(n_total_tracks);
        fTemp_MC_PiPlus.reserve(n_total_tracks);
    }

    // loop over all pre-selected tracks //
    for (std::size_t entry_track = 0; entry_track < n_total_tracks; ++entry_track) {
        const POD::Track& track = fInput.Track[entry_track];  // cache index lookup

        // PENDING: cache calculations to speed up cuts! maybe not needed? //

        // PID and pre-selection //
        if (track.Charge < 0) {
            if (PassesCuts_Proton(track, fHist_CutFlow_AntiProton.get())) {
                fTemp_AntiProton.emplace_back(track);
                if (fSettings.IsMC) {
                    fTemp_MC_AntiProton.emplace_back(
                        BuildMcTrack(fInput.Track_McEntry[entry_track], DB::Particles::Particle("AntiProton").pdg_code, true));
                }
            }
            if (PassesCuts_Kaon(track, fHist_CutFlow_NegKaon.get())) {
                fTemp_NegKaon.emplace_back(track);
                if (fSettings.IsMC) {
                    fTemp_MC_NegKaon.emplace_back(BuildMcTrack(fInput.Track_McEntry[entry_track], DB::Particles::Particle("NegKaon").pdg_code, true));
                }
            }
            if (PassesCuts_Pion(track, fHist_CutFlow_PiMinus.get())) {
                fTemp_PiMinus.emplace_back(track);
                if (fSettings.IsMC) {
                    fTemp_MC_PiMinus.emplace_back(BuildMcTrack(fInput.Track_McEntry[entry_track], DB::Particles::Particle("PiMinus").pdg_code, true));
                }
            }
        }
        if (track.Charge > 0) {
            if (PassesCuts_Proton(track, fHist_CutFlow_Proton.get())) {
                fTemp_Proton.emplace_back(track);
                if (fSettings.IsMC) {
                    fTemp_MC_Proton.emplace_back(BuildMcTrack(fInput.Track_McEntry[entry_track], DB::Particles::Particle("Proton").pdg_code, true));
                }
            }
            if (PassesCuts_Kaon(track, fHist_CutFlow_PosKaon.get())) {
                fTemp_PosKaon.emplace_back(track);
                if (fSettings.IsMC) {
                    fTemp_MC_PosKaon.emplace_back(BuildMcTrack(fInput.Track_McEntry[entry_track], DB::Particles::Particle("PosKaon").pdg_code, true));
                }
            }
            if (PassesCuts_Pion(track, fHist_CutFlow_PiPlus.get())) {
                fTemp_PiPlus.emplace_back(track);
                if (fSettings.IsMC) {
                    fTemp_MC_PiPlus.emplace_back(BuildMcTrack(fInput.Track_McEntry[entry_track], DB::Particles::Particle("PiPlus").pdg_code, true));
                }
            }
        }
    }  // end of loop over tracks

#if T2DS_DEBUG
    Logger::Debug(__FUNCTION__, "n_antiprotons = {}", fTemp_AntiProton.size());
    Logger::Debug(__FUNCTION__, "n_protons     = {}", fTemp_Proton.size());
    Logger::Debug(__FUNCTION__, "n_negkaons    = {}", fTemp_NegKaon.size());
    Logger::Debug(__FUNCTION__, "n_poskaons    = {}", fTemp_PosKaon.size());
    Logger::Debug(__FUNCTION__, "n_piminus     = {}", fTemp_PiMinus.size());
    Logger::Debug(__FUNCTION__, "n_piplus      = {}", fTemp_PiPlus.size());
#endif
}

bool Finder::PassesCuts_Proton(const POD::Track& track, TH1D* cut_flow_hist) const {
    FillHist(cut_flow_hist, EProton::kAllPossibleProtons);

    if (std::abs(static_cast<double>(track.NSigmasProton)) > T2DS::Cuts::Proton::AbsMax_NSigmasProton) return false;
    FillHist(cut_flow_hist, EProton::kPasses_NSigmasProtons);

    return true;
}

bool Finder::PassesCuts_Kaon(const POD::Track& track, TH1D* cut_flow_hist) const {
    FillHist(cut_flow_hist, EKaon::kAllPossibleKaons);

    if (std::abs(static_cast<double>(track.NSigmasKaon)) > T2DS::Cuts::Kaon::AbsMax_NSigmasKaon) return false;
    FillHist(cut_flow_hist, EKaon::kPasses_NSigmasKaons);

    return true;
}

bool Finder::PassesCuts_Pion(const POD::Track& track, TH1D* cut_flow_hist) const {
    FillHist(cut_flow_hist, EPion::kAllPossiblePions);

    if (std::abs(static_cast<double>(track.NSigmasPion)) > T2DS::Cuts::Pion::AbsMax_NSigmasPion) return false;
    FillHist(cut_flow_hist, EPion::kPasses_NSigmasPions);

    return true;
}

POD::Extended::McParticle Finder::BuildMcTrack(unsigned int track_mc_entry, int pdg_code_hypothesis, bool include_gm) {
    // copy linked mc info //
    POD::Extended::McParticle new_mc(fInput.McParticle[track_mc_entry]);
    auto c = MC::SexaquarkRules::ClassifyDownstream(new_mc, fInput.McParticle, fMcSignalChannel.value(), pdg_code_hypothesis, include_gm, false);
    MC::Apply(new_mc, c);
    return new_mc;
}

// ## V0s ZONE ## //

void Finder::FindV0s(const DB::Particles::Definition& pid) {

    // determine rules based on V0 species //
    const std::vector<POD::Track>* temp_vec_neg = &fTemp_PiMinus;
    const std::vector<POD::Track>* temp_vec_pos = &fTemp_PiPlus;
    const std::vector<POD::Extended::McParticle>* temp_vec_mc_neg = &fTemp_MC_PiMinus;
    const std::vector<POD::Extended::McParticle>* temp_vec_mc_pos = &fTemp_MC_PiPlus;
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
            temp_vec_neg = &fTemp_AntiProton;
            temp_vec_mc_neg = &fTemp_MC_AntiProton;
            pid_neg = DB::Particles::Particle("AntiProton");
            output_vec_v0 = &fTemp_AntiLambda;
            output_vec_v0_neg = &fTemp_AntiLambda_Neg;
            output_vec_v0_pos = &fTemp_AntiLambda_Pos;
            output_vec_mc_v0 = &fTemp_MC_AntiLambda;
            output_vec_mc_v0_neg = &fTemp_MC_AntiLambda_Neg;
            output_vec_mc_v0_pos = &fTemp_MC_AntiLambda_Pos;
            break;
        }
        case DB::Particles::Particle("Lambda").pdg_code: {
            temp_vec_pos = &fTemp_Proton;
            temp_vec_mc_pos = &fTemp_MC_Proton;
            pid_pos = DB::Particles::Particle("Proton");
            output_vec_v0 = &fTemp_Lambda;
            output_vec_v0_neg = &fTemp_Lambda_Neg;
            output_vec_v0_pos = &fTemp_Lambda_Pos;
            output_vec_mc_v0 = &fTemp_MC_Lambda;
            output_vec_mc_v0_neg = &fTemp_MC_Lambda_Neg;
            output_vec_mc_v0_pos = &fTemp_MC_Lambda_Pos;
            break;
        }
        case DB::Particles::Particle("KaonZeroShort").pdg_code: {
            output_vec_v0 = &fTemp_KaonZeroShort;
            output_vec_v0_neg = &fTemp_KaonZeroShort_Neg;
            output_vec_v0_pos = &fTemp_KaonZeroShort_Pos;
            output_vec_mc_v0 = &fTemp_MC_KaonZeroShort;
            output_vec_mc_v0_neg = &fTemp_MC_KaonZeroShort_Neg;
            output_vec_mc_v0_pos = &fTemp_MC_KaonZeroShort_Pos;
            break;
        }
        default: {
            Logger::Error(__FUNCTION__, "Invalid PID {} for a V0.", pid.name);
            return;
        }
    }

    // loop over all possible pairs of tracks //
    // NOTE: negative and positive species never share a track, hence no sanity check is needed
    for (std::size_t entry_neg = 0; entry_neg < temp_vec_neg->size(); ++entry_neg) {
        const POD::Track& track_neg = (*temp_vec_neg)[entry_neg];  // cache index lookup

        for (std::size_t entry_pos = 0; entry_pos < temp_vec_pos->size(); ++entry_pos) {
            const POD::Track& track_pos = (*temp_vec_pos)[entry_pos];  // cache index lookup

            // apply cuts (1) //
            // PENDING: placeholder to remove duplications

            // PCAs //
            auto [seed_neg, seed_pos, pca_cache] = Seeder::HelixHelix::FastCorrectPCAs(track_neg, track_pos, fMagneticField);

            // apply cuts (2) //
            if (!PostSeedCuts(seed_neg.pca, seed_pos.pca, pid)) continue;

            // PCAs derivatives //
            auto [deriv_neg, deriv_pos] = Seeder::HelixHelix::ComputeDerivatives(seed_neg, seed_pos, pca_cache);

            // fit vertex //
            auto fit = KF::FitVertex(track_neg, track_pos, pid_neg, pid_pos, {seed_neg, deriv_neg}, {seed_pos, deriv_pos}, fMagneticField);

            // create storage+computation units //
            POD::V0 v0 = CreateV0(fit, seed_neg.pca, seed_pos.pca);
            Cached::V0 c_v0(v0, fPrimaryVertex);

            // apply cuts (3) //
            if (!PostFitCuts(c_v0, pid)) continue;

            // store reconstructed //
            output_vec_v0->emplace_back(v0);
            output_vec_v0_neg->emplace_back(track_neg);
            output_vec_v0_pos->emplace_back(track_pos);

            // store mc //
            // NOTE: the daughters were already built in `ProcessTracks(...)`, under this very same pid hypothesis
            if (fSettings.IsMC) {
                // -- neg
                const POD::Extended::McParticle& mc_neg = (*temp_vec_mc_neg)[entry_neg];
                output_vec_mc_v0_neg->emplace_back(mc_neg);
                // -- pos
                const POD::Extended::McParticle& mc_pos = (*temp_vec_mc_pos)[entry_pos];
                output_vec_mc_v0_pos->emplace_back(mc_pos);
                // -- v0
                output_vec_mc_v0->emplace_back(BuildMcV0(mc_neg, mc_pos, pid.pdg_code));
            }
        }  // end of loop over pos
    }  // end of loop over neg
}

bool Finder::PostSeedCuts_Lambda(const Seeder::PCA& pca_neg, const Seeder::PCA& pca_pos, TH1D* cut_flow_hist) const {
    FillHist(cut_flow_hist, ELambda::kAllCombinations);

    if (CMath::SquaredDistance(pca_neg.xyz, pca_pos.xyz) > Cuts::Lambda::Max_DCAbtwDau * Cuts::Lambda::Max_DCAbtwDau) {
        return false;
    }
    FillHist(cut_flow_hist, ELambda::kPasses_DcaBtwDaughters);

    return true;
}

bool Finder::PostSeedCuts_KaonZeroShort(const Seeder::PCA& pca_neg, const Seeder::PCA& pca_pos, TH1D* cut_flow_hist) const {
    FillHist(cut_flow_hist, EKaonZeroShort::kAllCombinations);

    if (CMath::SquaredDistance(pca_neg.xyz, pca_pos.xyz) > T2DS::Cuts::KaonZeroShort::Max_DCAbtwDau * T2DS::Cuts::KaonZeroShort::Max_DCAbtwDau) {
        return false;
    }
    FillHist(cut_flow_hist, EKaonZeroShort::kPasses_DcaBtwDaughters);

    return true;
}

bool Finder::PostFitCuts_Lambda(const Cached::V0& c_v0, TH1D* cut_flow_hist) const {

    double mass = c_v0.Mass();  // cached
    if (mass < Cuts::Lambda::Min_Mass || mass > Cuts::Lambda::Max_Mass) return false;
    // FillHist(cut_flow_hist, 2.);  // PENDING

    if (c_v0.Decay_SquaredRadius2D() < Cuts::Lambda::Min_Decay_Radius2D * Cuts::Lambda::Min_Decay_Radius2D) return false;
    // FillHist(cut_flow_hist, 3.);  // PENDING

    if (c_v0.Neg_SquaredDCA_wrt_V0() > Cuts::Lambda::Max_DCAnegV0 * Cuts::Lambda::Max_DCAnegV0) return false;
    // FillHist(cut_flow_hist, 4.);  // PENDING

    if (c_v0.Pos_SquaredDCA_wrt_V0() > Cuts::Lambda::Max_DCAposV0 * Cuts::Lambda::Max_DCAposV0) return false;
    // FillHist(cut_flow_hist, 5.);  // PENDING

    // if (c_v0.Pt() < Cuts::Lambda::Min_Pt) return false; // PENDING
    // FillHist(cut_flow_hist, 6.);  // PENDING

    if (std::abs(c_v0.Rapidity()) > Cuts::Lambda::AbsMax_Rapidity) return false;
    // FillHist(cut_flow_hist, 7.);  // PENDING

    // if (c_v0.AbsArmQtOverAlpha() > Cuts::Lambda::AbsMax_ArmQtOverAlpha) return false;  // PENDING: not really sure if i like this cut, actually
    // FillHist(cut_flow_hist, 8.);  // PENDING

    if (c_v0.CPA_wrt_PV() < Cuts::Lambda::Min_CPAwrtPV || c_v0.CPA_wrt_PV() > Cuts::Lambda::Max_CPAwrtPV) return false;
    // FillHist(cut_flow_hist, 8.);  // PENDING

    if (c_v0.SquaredDCA_wrt_PV() < Cuts::Lambda::Min_DCAwrtPV * Cuts::Lambda::Min_DCAwrtPV) return false;
    // FillHist(cut_flow_hist, 9.);  // PENDING

    return true;
}

bool Finder::PostFitCuts_KaonZeroShort(const Cached::V0& c_v0, TH1D* cut_flow_hist) const {

    // if (c_v0.Pt() < Cuts::KaonZeroShort::Min_Pt) return false; // PENDING
    // FillHist(cut_flow_hist, 2.);  // PENDING // PENDING

    double mass = c_v0.Mass();  // cached
    if (mass < Cuts::KaonZeroShort::Min_Mass || mass > Cuts::KaonZeroShort::Max_Mass) return false;
    // FillHist(cut_flow_hist, 3.); // PENDING

    if (std::abs(c_v0.Rapidity()) > Cuts::KaonZeroShort::AbsMax_Rapidity) return false;
    // FillHist(cut_flow_hist, 4.); // PENDING

    if (c_v0.Decay_SquaredRadius2D() < Cuts::KaonZeroShort::Min_Decay_Radius2D * Cuts::KaonZeroShort::Min_Decay_Radius2D) return false;
    // FillHist(cut_flow_hist, 5.); // PENDING

    if (c_v0.Neg_SquaredDCA_wrt_V0() > Cuts::KaonZeroShort::Max_DCAnegV0 * Cuts::KaonZeroShort::Max_DCAnegV0) return false;
    // FillHist(cut_flow_hist, 6.); // PENDING

    if (c_v0.Pos_SquaredDCA_wrt_V0() > Cuts::KaonZeroShort::Max_DCAposV0 * Cuts::KaonZeroShort::Max_DCAposV0) return false;
    // FillHist(cut_flow_hist, 7.); // PENDING

    if (c_v0.CPA_wrt_PV() < Cuts::KaonZeroShort::Min_CPAwrtPV || c_v0.CPA_wrt_PV() > Cuts::KaonZeroShort::Max_CPAwrtPV) return false;
    // FillHist(cut_flow_hist, 8.); // PENDING

    if (c_v0.SquaredDCA_wrt_PV() < Cuts::KaonZeroShort::Min_DCAwrtPV * Cuts::KaonZeroShort::Min_DCAwrtPV) return false;
    // FillHist(cut_flow_hist, 9.); // PENDING

    return true;
}

POD::Extended::McParticle Finder::BuildMcV0(const POD::Extended::McParticle& mc_neg, const POD::Extended::McParticle& mc_pos,
                                            int pdg_code_hypothesis) {
    POD::Extended::McParticle mc_v0;

    // -- fill hybridness, independetly of no common mother
    mc_v0.IsHybrid = mc_neg.IsTrueSignal != mc_pos.IsTrueSignal || mc_neg.SignalID != mc_pos.SignalID;

    auto mc_entry = MC::FindCommonMotherMcEntry(mc_neg, mc_pos);
    if (!mc_entry.has_value()) return mc_v0;

    // fill values //
    static_cast<POD::McParticle&>(mc_v0) = fInput.McParticle[mc_entry.value()];
    MC::Apply(mc_v0, MC::SexaquarkRules::ClassifyDownstream(mc_v0, fInput.McParticle, fMcSignalChannel.value(), pdg_code_hypothesis, false, true));

    return mc_v0;
}

POD::V0 Finder::CreateV0(const KF::FitResult& fit, const Seeder::PCA& neg_pca_wrt_v0, const Seeder::PCA& pos_pca_wrt_v0) {
    POD::V0 new_v0;  // non-initialized on purpose
    new_v0.Decay_X = static_cast<float>(fit.mother.X());
    new_v0.Decay_Y = static_cast<float>(fit.mother.Y());
    new_v0.Decay_Z = static_cast<float>(fit.mother.Z());
    new_v0.Px = static_cast<float>(fit.mother.Px());
    new_v0.Py = static_cast<float>(fit.mother.Py());
    new_v0.Pz = static_cast<float>(fit.mother.Pz());
    new_v0.Energy = static_cast<float>(fit.mother.E());
    new_v0.CovMatrix = fit.mother.Cov<7>();
    new_v0.Chi2NDF = static_cast<float>(fit.mother.Chi2NDF());
    // -- negative daughter
    new_v0.Neg_PCAwrtV0_X = static_cast<float>(neg_pca_wrt_v0.X());
    new_v0.Neg_PCAwrtV0_Y = static_cast<float>(neg_pca_wrt_v0.Y());
    new_v0.Neg_PCAwrtV0_Z = static_cast<float>(neg_pca_wrt_v0.Z());
    new_v0.Neg_Fit_Px = static_cast<float>(fit.Dau1_Px());
    new_v0.Neg_Fit_Py = static_cast<float>(fit.Dau1_Py());
    new_v0.Neg_Fit_Pz = static_cast<float>(fit.Dau1_Pz());
    new_v0.Neg_Fit_Energy = static_cast<float>(fit.Dau1_E());
    // -- positive daughter
    new_v0.Pos_PCAwrtV0_X = static_cast<float>(pos_pca_wrt_v0.X());
    new_v0.Pos_PCAwrtV0_Y = static_cast<float>(pos_pca_wrt_v0.Y());
    new_v0.Pos_PCAwrtV0_Z = static_cast<float>(pos_pca_wrt_v0.Z());
    new_v0.Pos_Fit_Px = static_cast<float>(fit.Dau2_Px());
    new_v0.Pos_Fit_Py = static_cast<float>(fit.Dau2_Py());
    new_v0.Pos_Fit_Pz = static_cast<float>(fit.Dau2_Pz());
    new_v0.Pos_Fit_Energy = static_cast<float>(fit.Dau2_E());

    return new_v0;
}

// ## Channel A ZONE ## //

void Finder::FindSexaquarks_ChannelA(bool is_bkg_channel) {

    // determine rules and aliases //
    // (anti)lambdas
    // -- rec
    const auto& input_lambdas = is_bkg_channel ? fTemp_Lambda : fTemp_AntiLambda;
    const auto& input_lambdas_neg = is_bkg_channel ? fTemp_Lambda_Neg : fTemp_AntiLambda_Neg;
    const auto& input_lambdas_pos = is_bkg_channel ? fTemp_Lambda_Pos : fTemp_AntiLambda_Pos;
    const std::size_t n_lambdas = input_lambdas.size();
    // -- mc
    const std::vector<POD::Extended::McParticle>* input_mc_lambdas = nullptr;
    const std::vector<POD::Extended::McParticle>* input_mc_lambdas_neg = nullptr;
    const std::vector<POD::Extended::McParticle>* input_mc_lambdas_pos = nullptr;
    if (fSettings.IsMC) {
        input_mc_lambdas = is_bkg_channel ? &fTemp_MC_Lambda : &fTemp_MC_AntiLambda;
        input_mc_lambdas_neg = is_bkg_channel ? &fTemp_MC_Lambda_Neg : &fTemp_MC_AntiLambda_Neg;
        input_mc_lambdas_pos = is_bkg_channel ? &fTemp_MC_Lambda_Pos : &fTemp_MC_AntiLambda_Pos;
    }
    // kaon-zero-short
    // -- rec
    const auto& input_k0s = fTemp_KaonZeroShort;
    const auto& input_k0s_neg = fTemp_KaonZeroShort_Neg;
    const auto& input_k0s_pos = fTemp_KaonZeroShort_Pos;
    const std::size_t n_k0s = input_k0s.size();
    // -- mc
    const std::vector<POD::Extended::McParticle>* input_mc_k0s = nullptr;
    const std::vector<POD::Extended::McParticle>* input_mc_k0s_neg = nullptr;
    const std::vector<POD::Extended::McParticle>* input_mc_k0s_pos = nullptr;
    if (fSettings.IsMC) {
        input_mc_k0s = &fTemp_MC_KaonZeroShort;
        input_mc_k0s_neg = &fTemp_MC_KaonZeroShort_Neg;
        input_mc_k0s_pos = &fTemp_MC_KaonZeroShort_Pos;
    }
    // cut flow hists
    TH1D* hist = is_bkg_channel ? fHist_CutFlow_ChannelA_Bkg.get() : fHist_CutFlow_ChannelA.get();
    // daughter hypotheses, which is where the fit reads their masses from
    const DB::Particles::Definition pid_lambda = is_bkg_channel ? DB::Particles::Particle("Lambda") : DB::Particles::Particle("AntiLambda");
    constexpr DB::Particles::Definition pid_k0s = DB::Particles::Particle("KaonZeroShort");

    // loop over all possible pairs of (anti)lambda + K0S //
    for (std::size_t entry_lambda = 0; entry_lambda < n_lambdas; ++entry_lambda) {
        // cache index lookups //
        const POD::V0& lambda = input_lambdas[entry_lambda];
        const POD::Track& lambda_neg = input_lambdas_neg[entry_lambda];
        const POD::Track& lambda_pos = input_lambdas_pos[entry_lambda];

        for (std::size_t entry_k0s = 0; entry_k0s < n_k0s; ++entry_k0s) {
            // cache index lookups //
            const POD::V0& k0s = input_k0s[entry_k0s];
            const POD::Track& k0s_neg = input_k0s_neg[entry_k0s];
            const POD::Track& k0s_pos = input_k0s_pos[entry_k0s];

            // sanity check //
            if (lambda_neg.EsdEntry == k0s_neg.EsdEntry || lambda_neg.EsdEntry == k0s_pos.EsdEntry || lambda_pos.EsdEntry == k0s_neg.EsdEntry ||
                lambda_pos.EsdEntry == k0s_pos.EsdEntry) {
                continue;
            }

            // PCAs //
            Seeder::LineLine::Cache pca_cache;
            auto [seed_lambda, seed_k0s] = Seeder::LineLine::FastPCAs(lambda, k0s, &pca_cache);

            // apply cuts (1) //
            if (!PostSeedCuts_ChannelA(seed_lambda.pca, seed_k0s.pca, hist)) continue;

            // PCAs derivatives //
            auto [deriv_lambda, deriv_k0s] = Seeder::LineLine::ComputeDerivatives(pca_cache);

            // fit vertex //
            auto fit = KF::FitVertex(lambda, k0s, pid_lambda, pid_k0s, {seed_lambda, deriv_lambda}, {seed_k0s, deriv_k0s});

            // create storage+computation units //
            POD::Sexaquark sexa = Create_ChannelA(fit, seed_lambda.pca, seed_k0s.pca, is_bkg_channel);
            Cached::Sexaquark c_sexa(sexa, fPrimaryVertex);

            // apply cuts (2) //
            if (!PostFitCuts_ChannelA(c_sexa, hist)) continue;

            // store reconstructed //
            fOutput.ChannelA.emplace_back(sexa);
            fOutput.ChannelA_V0A.emplace_back(lambda);
            fOutput.ChannelA_V0A_Neg.emplace_back(lambda_neg);
            fOutput.ChannelA_V0A_Pos.emplace_back(lambda_pos);
            fOutput.ChannelA_V0B.emplace_back(k0s);
            fOutput.ChannelA_V0B_Neg.emplace_back(k0s_neg);
            fOutput.ChannelA_V0B_Pos.emplace_back(k0s_pos);

            // store mc //
            if (fSettings.IsMC) {
                // -- V0A
                const auto& mc_lambda = (*input_mc_lambdas)[entry_lambda];
                fOutput.MC_ChannelA_V0A.emplace_back(mc_lambda);
                fOutput.MC_ChannelA_V0A_Neg.emplace_back((*input_mc_lambdas_neg)[entry_lambda]);
                fOutput.MC_ChannelA_V0A_Pos.emplace_back((*input_mc_lambdas_pos)[entry_lambda]);
                // -- V0B
                const auto& mc_k0s = (*input_mc_k0s)[entry_k0s];
                fOutput.MC_ChannelA_V0B.emplace_back(mc_k0s);
                fOutput.MC_ChannelA_V0B_Neg.emplace_back((*input_mc_k0s_neg)[entry_k0s]);
                fOutput.MC_ChannelA_V0B_Pos.emplace_back((*input_mc_k0s_pos)[entry_k0s]);
                // -- h-dibaryon
                fOutput.MC_ChannelA.emplace_back(BuildMcSexaquark(mc_lambda, mc_k0s));
            }
        }
    }
}

bool Finder::PostSeedCuts_ChannelA(const Seeder::PCA& pca_v0a, const Seeder::PCA& pca_v0b, TH1D* hist_cut_flow) const {
    FillHist(hist_cut_flow, EChannelA::kAllCombinations);

    // if (Common::Math::SquaredDistance(pca_v0a.xyz, pca_v0b.xyz) > T2DS::Cuts::ChannelA::Max_DCAbtwV0s * T2DS::Cuts::ChannelA::Max_DCAbtwV0s) {
    // return false;
    // }
    // FillHist(hist_cut_flow, 1.);

    return true;
}

bool Finder::PostFitCuts_ChannelA(const Cached::Sexaquark& c_sexa, TH1D* hist_cut_flow) const {

    // if (c_sexa.SV_SquaredRadius2D() < Cuts::ChannelA::Min_Radius2D * Cuts::ChannelA::Min_Radius2D) return false;
    // FillHist(hist_cut_flow, 2.); // PENDING

    // if (c_sexa.Dau1_SquaredDCA_wrt_SV() > Cuts::ChannelA::Max_DCALaSV * Cuts::ChannelA::Max_DCALaSV) return false;
    // FillHist(hist_cut_flow, 3.); // PENDING

    // if (c_sexa.Dau2_SquaredDCA_wrt_SV() > Cuts::ChannelA::Max_DCAK0SV * Cuts::ChannelA::Max_DCAK0SV) return false;
    // FillHist(hist_cut_flow, 4.); // PENDING

    // if (c_sexa.CPA_wrt(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) < Cuts::ChannelA::Min_CPAwrtPV) return false;
    // FillHist(hist_cut_flow, 5.); // PENDING

    // if (c_sexa.DCA_wrt(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) > Cuts::ChannelA::Max_DCAwrtPV) return false;
    // FillHist(hist_cut_flow, 6.); // PENDING

    // if (c_sexa.CPA_V0A_SV() < Cuts::ChannelA::Min_La_CPAwrtSV) return false;
    // FillHist(hist_cut_flow, 7.); // PENDING

    // if (c_sexa.CPA_V0B_SV() < Cuts::ChannelA::Min_K0S_CPAwrtSV) return false;
    // FillHist(hist_cut_flow, 8.); // PENDING

    return true;
}

POD::Sexaquark Finder::Create_ChannelA(const KF::FitResult& fit, const Seeder::PCA& pca_v0a, const Seeder::PCA& pca_v0b, bool is_bkg_channel) {
    POD::Sexaquark sexa;  // non-initialized on purpose
    sexa.SV_X = static_cast<float>(fit.mother.X());
    sexa.SV_Y = static_cast<float>(fit.mother.Y());
    sexa.SV_Z = static_cast<float>(fit.mother.Z());
    sexa.Px = static_cast<float>(fit.mother.Px());
    sexa.Py = static_cast<float>(fit.mother.Py());
    sexa.Pz = static_cast<float>(fit.mother.Pz());
    sexa.Energy = static_cast<float>(fit.mother.E());
    sexa.Chi2NDF = static_cast<float>(fit.mother.Chi2NDF());
    sexa.E_MinusNucleon = static_cast<float>(fit.mother.E() - DB::Particles::Particle("Neutron").mass);  // small optimization
    sexa.IsBkgChannel = is_bkg_channel;
    // -- V0A
    sexa.Dau1_PCAwrtSV_X = static_cast<float>(pca_v0a.X());
    sexa.Dau1_PCAwrtSV_Y = static_cast<float>(pca_v0a.Y());
    sexa.Dau1_PCAwrtSV_Z = static_cast<float>(pca_v0a.Z());
    sexa.Dau1_Fit_Px = static_cast<float>(fit.Dau1_Px());
    sexa.Dau1_Fit_Py = static_cast<float>(fit.Dau1_Py());
    sexa.Dau1_Fit_Pz = static_cast<float>(fit.Dau1_Pz());
    sexa.Dau1_Fit_Energy = static_cast<float>(fit.Dau1_E());
    // -- V0B
    sexa.Dau2_PCAwrtSV_X = static_cast<float>(pca_v0b.X());
    sexa.Dau2_PCAwrtSV_Y = static_cast<float>(pca_v0b.Y());
    sexa.Dau2_PCAwrtSV_Z = static_cast<float>(pca_v0b.Z());
    sexa.Dau2_Fit_Px = static_cast<float>(fit.Dau2_Px());
    sexa.Dau2_Fit_Py = static_cast<float>(fit.Dau2_Py());
    sexa.Dau2_Fit_Pz = static_cast<float>(fit.Dau2_Pz());
    sexa.Dau2_Fit_Energy = static_cast<float>(fit.Dau2_E());

    return sexa;
}

// ## Channel D ZONE ## //

void Finder::FindSexaquarks_ChannelD(bool is_bkg_channel) {

    // determine rules and aliases //
    // -- lambda
    //    -- rec
    const auto& input_lambdas = is_bkg_channel ? fTemp_Lambda : fTemp_AntiLambda;
    const auto& input_lambdas_neg = is_bkg_channel ? fTemp_Lambda_Neg : fTemp_AntiLambda_Neg;
    const auto& input_lambdas_pos = is_bkg_channel ? fTemp_Lambda_Pos : fTemp_AntiLambda_Pos;
    const std::size_t n_lambdas = input_lambdas.size();
    //    -- mc
    const std::vector<POD::Extended::McParticle>* input_mc_lambdas = nullptr;
    const std::vector<POD::Extended::McParticle>* input_mc_lambdas_neg = nullptr;
    const std::vector<POD::Extended::McParticle>* input_mc_lambdas_pos = nullptr;
    if (fSettings.IsMC) {
        input_mc_lambdas = is_bkg_channel ? &fTemp_MC_Lambda : &fTemp_MC_AntiLambda;
        input_mc_lambdas_neg = is_bkg_channel ? &fTemp_MC_Lambda_Neg : &fTemp_MC_AntiLambda_Neg;
        input_mc_lambdas_pos = is_bkg_channel ? &fTemp_MC_Lambda_Pos : &fTemp_MC_AntiLambda_Pos;
    }
    // -- kaon
    //    -- rec
    const auto& input_kaons = is_bkg_channel ? fTemp_NegKaon : fTemp_PosKaon;
    const std::size_t n_kaons = input_kaons.size();
    //    -- mc
    const std::vector<POD::Extended::McParticle>* input_mc_kaons = nullptr;
    if (fSettings.IsMC) {
        input_mc_kaons = is_bkg_channel ? &fTemp_MC_NegKaon : &fTemp_MC_PosKaon;
    }
    // -- daughter hypotheses
    const DB::Particles::Definition pid_kaon = is_bkg_channel ? DB::Particles::Particle("NegKaon") : DB::Particles::Particle("PosKaon");
    const DB::Particles::Definition pid_lambda = is_bkg_channel ? DB::Particles::Particle("Lambda") : DB::Particles::Particle("AntiLambda");
    // -- cut flow hist
    TH1D* hist = is_bkg_channel ? fHist_CutFlow_ChannelD_Bkg.get() : fHist_CutFlow_ChannelD.get();

    // loop over all possible pairs of (anti)lambda + (pos/neg)kaon //
    for (std::size_t entry_lambda = 0; entry_lambda < n_lambdas; ++entry_lambda) {
        // cache index lookups //
        const POD::V0& lambda = input_lambdas[entry_lambda];
        const POD::Track& lambda_neg = input_lambdas_neg[entry_lambda];
        const POD::Track& lambda_pos = input_lambdas_pos[entry_lambda];

        for (std::size_t entry_kaon = 0; entry_kaon < n_kaons; ++entry_kaon) {
            // cache index lookup //
            const POD::Track& kaon = input_kaons[entry_kaon];

            // -- sanity check
            if (lambda_neg.EsdEntry == kaon.EsdEntry || lambda_pos.EsdEntry == kaon.EsdEntry) continue;

            // PCAs (1) //
            auto [seed_kaon, seed_v0, pca_cache] = Seeder::HelixLine::FastCorrectPCAs(kaon, lambda, fMagneticField);

            // apply cuts (1) //
            if (!PostSeedCuts_ChannelD(seed_v0.pca, seed_kaon.pca, hist)) continue;

            // PCAs derivatives //
            auto [deriv_ka, deriv_v0] = Seeder::HelixLine::ComputeDerivatives(seed_kaon, seed_v0, pca_cache);

            // fit vertex //
            auto fit = KF::FitVertex(kaon, lambda, pid_kaon, pid_lambda, {seed_kaon, deriv_ka}, {seed_v0, deriv_v0}, fMagneticField);

            // create storage+computation units //
            POD::Sexaquark sexa = Create_ChannelD(fit, seed_v0.pca, seed_kaon.pca, is_bkg_channel);
            Cached::Sexaquark c_sexa(sexa, fPrimaryVertex);

            // apply cuts (2) //
            if (!PostFitCuts_ChannelD(c_sexa, hist)) continue;

            // store reconstructed //
            fOutput.ChannelD.emplace_back(sexa);
            fOutput.ChannelD_V0.emplace_back(lambda);
            fOutput.ChannelD_V0_Neg.emplace_back(lambda_neg);
            fOutput.ChannelD_V0_Pos.emplace_back(lambda_pos);
            fOutput.ChannelD_Kaon.emplace_back(kaon);

            // store mc //
            if (fSettings.IsMC) {
                // -- V0
                const auto& mc_lambda = (*input_mc_lambdas)[entry_lambda];
                fOutput.MC_ChannelD_V0.emplace_back(mc_lambda);
                fOutput.MC_ChannelD_V0_Neg.emplace_back((*input_mc_lambdas_neg)[entry_lambda]);
                fOutput.MC_ChannelD_V0_Pos.emplace_back((*input_mc_lambdas_pos)[entry_lambda]);
                // -- Kaon
                const auto& mc_kaon = (*input_mc_kaons)[entry_kaon];
                fOutput.MC_ChannelD_Kaon.emplace_back(mc_kaon);
                // -- h-dibaryon
                fOutput.MC_ChannelD.emplace_back(BuildMcSexaquark(mc_lambda, mc_kaon));
            }
        }
    }
}

bool Finder::PostSeedCuts_ChannelD(const Seeder::PCA& pca_v0, const Seeder::PCA& pca_ka, TH1D* hist_cut_flow) const {
    FillHist(hist_cut_flow, EChannelD::kAllCombinations);

    // if (Common::Math::SquaredDistance(pca_ka.xyz, pca_v0.xyz) > Cuts::ChannelD::Max_DCAKaLa * Cuts::ChannelD::Max_DCAKaLa) return false;
    // FillHist(hist_cut_flow, 1.); // PENDING

    return true;
}

bool Finder::PostFitCuts_ChannelD(const Cached::Sexaquark& c_sexa, TH1D* hist_cut_flow) const {

    // double sq_radius_2d = c_sexa.SV_SquaredRadius2D();
    // if (sq_radius_2d < Cuts::ChannelD::Min_Radius2D * Cuts::ChannelD::Min_Radius2D ||
    // sq_radius_2d > Cuts::ChannelD::Max_Radius2D * Cuts::ChannelD::Max_Radius2D) {
    // return false;
    // }
    // FillHist(hist_cut_flow, 2.); // PENDING

    // if (sexa.AbsRapidity_MinusNucleon() > Cuts::ChannelD::AbsMax_Rapidity) return false;  // PENDING: kinematics, affected by Fermi motion
    // FillHist(hist_cut_flow, 3.); // PENDING

    // if (sexa.CPA_Vertex(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) < Cuts::ChannelD::Min_CPAwrtPV ||
    // sexa.CPA_Vertex(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) > Cuts::ChannelD::Max_CPAwrtPV) {
    // return false;  // PENDING: kinematics, affected by Fermi motion
    // }
    // FillHist(hist_cut_flow, 3.); // PENDING

    // if (c_sexa.Dau1_SquaredDCA_wrt_SV() > Cuts::ChannelD::Max_DCALaSV * Cuts::ChannelD::Max_DCALaSV) return false;
    // FillHist(hist_cut_flow, 4.); // PENDING

    // if (c_sexa.Dau2_SquaredDCA_wrt_SV() > Cuts::ChannelD::Max_DCAKaSV * Cuts::ChannelD::Max_DCAKaSV) return false;
    // FillHist(hist_cut_flow, 5.); // PENDING

    // if (sexa.DCA_V0Neg_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelD::Max_DCALaNegSV) return false;
    // FillHist(hist_cut_flow, 6.); // PENDING

    // if (sexa.DCA_V0Pos_wrt_SV(fInput_Event.MagneticField) > Cuts::ChannelD::Max_DCALaPosSV) return false;
    // FillHist(hist_cut_flow, 7.); // PENDING

    return true;
}

POD::Sexaquark Finder::Create_ChannelD(const KF::FitResult& fit, const Seeder::PCA& pca_v0, const Seeder::PCA& pca_ka, bool is_bkg_channel) {
    POD::Sexaquark sexa;  // non-initialized on purpose
    sexa.SV_X = static_cast<float>(fit.mother.X());
    sexa.SV_Y = static_cast<float>(fit.mother.Y());
    sexa.SV_Z = static_cast<float>(fit.mother.Z());
    sexa.Px = static_cast<float>(fit.mother.Px());
    sexa.Py = static_cast<float>(fit.mother.Py());
    sexa.Pz = static_cast<float>(fit.mother.Pz());
    sexa.Energy = static_cast<float>(fit.mother.E());
    sexa.Chi2NDF = static_cast<float>(fit.mother.Chi2NDF());
    sexa.E_MinusNucleon = static_cast<float>(fit.mother.E() - DB::Particles::Particle("Proton").mass);  // small optimization
    sexa.IsBkgChannel = is_bkg_channel;
    // V0
    sexa.Dau1_PCAwrtSV_X = static_cast<float>(pca_v0.X());
    sexa.Dau1_PCAwrtSV_Y = static_cast<float>(pca_v0.Y());
    sexa.Dau1_PCAwrtSV_Z = static_cast<float>(pca_v0.Z());
    sexa.Dau1_Fit_Px = static_cast<float>(fit.Dau2_Px());
    sexa.Dau1_Fit_Py = static_cast<float>(fit.Dau2_Py());
    sexa.Dau1_Fit_Pz = static_cast<float>(fit.Dau2_Pz());
    sexa.Dau1_Fit_Energy = static_cast<float>(fit.Dau2_E());
    // kaon
    sexa.Dau2_PCAwrtSV_X = static_cast<float>(pca_ka.X());
    sexa.Dau2_PCAwrtSV_Y = static_cast<float>(pca_ka.Y());
    sexa.Dau2_PCAwrtSV_Z = static_cast<float>(pca_ka.Z());
    sexa.Dau2_Fit_Px = static_cast<float>(fit.Dau1_Px());
    sexa.Dau2_Fit_Py = static_cast<float>(fit.Dau1_Py());
    sexa.Dau2_Fit_Pz = static_cast<float>(fit.Dau1_Pz());
    sexa.Dau2_Fit_Energy = static_cast<float>(fit.Dau1_E());

    return sexa;
}

// ## Channel H ZONE ## //

void Finder::FindSexaquarks_ChannelH(bool is_bkg_channel) {

    // determine rules and aliases //
    // kaons
    // -- rec
    const auto& input_kaons = is_bkg_channel ? fTemp_NegKaon : fTemp_PosKaon;
    const std::size_t n_kaons = input_kaons.size();
    // -- mc
    const std::vector<POD::Extended::McParticle>* input_mc_kaons = nullptr;
    if (fSettings.IsMC) input_mc_kaons = is_bkg_channel ? &fTemp_MC_NegKaon : &fTemp_MC_PosKaon;
    // -- daughter hypothesis, which is where the fit reads their mass from
    const DB::Particles::Definition pid_kaon = is_bkg_channel ? DB::Particles::Particle("NegKaon") : DB::Particles::Particle("PosKaon");
    // cut flow hist
    TH1D* hist = is_bkg_channel ? fHist_CutFlow_ChannelH_Bkg.get() : fHist_CutFlow_ChannelH.get();

    // loop over all possible pairs of (pos)kaon+(pos)kaon or (neg)kaon+(neg)kaon //
    for (std::size_t entry_kaon1 = 0; entry_kaon1 + 1 < n_kaons; ++entry_kaon1) {
        const POD::Track& kaon1 = input_kaons[entry_kaon1];  // cache index lookup

        for (std::size_t entry_kaon2 = entry_kaon1 + 1; entry_kaon2 < n_kaons; ++entry_kaon2) {
            // NOTE: sanity check not needed, because loops don't intersect
            const POD::Track& kaon2 = input_kaons[entry_kaon2];  // cache index lookup

            // PCAs (1) //
            auto [seed_kaon1, seed_kaon2, pca_cache] = Seeder::HelixHelix::FastCorrectPCAs(kaon1, kaon2, fMagneticField);

            // apply cuts (1) //
            if (!PostSeedCuts_ChannelH(seed_kaon1.pca, seed_kaon2.pca, hist)) continue;

            // PCAs derivatives //
            auto [deriv_kaon1, deriv_kaon2] = Seeder::HelixHelix::ComputeDerivatives(seed_kaon1, seed_kaon2, pca_cache);

            // fit vertex //
            auto fit = KF::FitVertex(kaon1, kaon2, pid_kaon, pid_kaon, {seed_kaon1, deriv_kaon1}, {seed_kaon2, deriv_kaon2}, fMagneticField);

            // create storage+computation units //
            POD::Sexaquark sexa = Create_ChannelH(fit, seed_kaon1.pca, seed_kaon2.pca, is_bkg_channel);
            Cached::Sexaquark c_sexa(sexa, fPrimaryVertex);

            // apply cuts (2) //
            if (!PostFitCuts_ChannelH(c_sexa, hist)) continue;

            // store reconstructed //
            fOutput.ChannelH.emplace_back(sexa);
            fOutput.ChannelH_Kaon1.emplace_back(kaon1);
            fOutput.ChannelH_Kaon2.emplace_back(kaon2);

            // store mc //
            if (fSettings.IsMC) {
                // -- Kaon1
                const auto& mc_kaon1 = (*input_mc_kaons)[entry_kaon1];
                fOutput.MC_ChannelH_Kaon1.emplace_back(mc_kaon1);
                // -- Kaon2
                const auto& mc_kaon2 = (*input_mc_kaons)[entry_kaon2];
                fOutput.MC_ChannelH_Kaon2.emplace_back(mc_kaon2);
                // -- h-dibaryon
                fOutput.MC_ChannelH.emplace_back(BuildMcSexaquark(mc_kaon1, mc_kaon2));
            }
        }
    }
}

bool Finder::PostSeedCuts_ChannelH(const Seeder::PCA& pca_kaon1, const Seeder::PCA& pca_kaon2, TH1D* hist_cut_flow) const {
    FillHist(hist_cut_flow, EChannelH::kAllCombinations);

    // PENDING //

    return true;
}

bool Finder::PostFitCuts_ChannelH(const Cached::Sexaquark& c_sexa, TH1D* hist_cut_flow) const {

    // PENDING //

    return true;
}

POD::Sexaquark Finder::Create_ChannelH(const KF::FitResult& fit, const Seeder::PCA& pca_kaon1, const Seeder::PCA& pca_kaon2, bool is_bkg_channel) {
    POD::Sexaquark sexa;  // non-initialized on purpose
    sexa.SV_X = static_cast<float>(fit.mother.X());
    sexa.SV_Y = static_cast<float>(fit.mother.Y());
    sexa.SV_Z = static_cast<float>(fit.mother.Z());
    sexa.Px = static_cast<float>(fit.mother.Px());
    sexa.Py = static_cast<float>(fit.mother.Py());
    sexa.Pz = static_cast<float>(fit.mother.Pz());
    sexa.Energy = static_cast<float>(fit.mother.E());
    sexa.Chi2NDF = static_cast<float>(fit.mother.Chi2NDF());
    sexa.E_MinusNucleon = static_cast<float>(fit.mother.E() - DB::Particles::Particle("Proton").mass);  // small optimization
    sexa.IsBkgChannel = is_bkg_channel;
    // -- Kaon 1
    sexa.Dau1_PCAwrtSV_X = static_cast<float>(pca_kaon1.X());
    sexa.Dau1_PCAwrtSV_Y = static_cast<float>(pca_kaon1.Y());
    sexa.Dau1_PCAwrtSV_Z = static_cast<float>(pca_kaon1.Z());
    sexa.Dau1_Fit_Px = static_cast<float>(fit.Dau1_Px());
    sexa.Dau1_Fit_Py = static_cast<float>(fit.Dau1_Py());
    sexa.Dau1_Fit_Pz = static_cast<float>(fit.Dau1_Pz());
    sexa.Dau1_Fit_Energy = static_cast<float>(fit.Dau1_E());
    // -- Kaon 2
    sexa.Dau2_PCAwrtSV_X = static_cast<float>(pca_kaon2.X());
    sexa.Dau2_PCAwrtSV_Y = static_cast<float>(pca_kaon2.Y());
    sexa.Dau2_PCAwrtSV_Z = static_cast<float>(pca_kaon2.Z());
    sexa.Dau2_Fit_Px = static_cast<float>(fit.Dau2_Px());
    sexa.Dau2_Fit_Py = static_cast<float>(fit.Dau2_Py());
    sexa.Dau2_Fit_Pz = static_cast<float>(fit.Dau2_Pz());
    sexa.Dau2_Fit_Energy = static_cast<float>(fit.Dau2_E());

    return sexa;
}

// ## END OF CYCLES ## //

void Finder::EndOfEvent() {

    // in case of data, don't keep event with no candidates
    const bool has_rec_candidates = !fOutput.ChannelA.empty() || !fOutput.ChannelD.empty() || !fOutput.ChannelH.empty();
    // in case of MC, keep event with injected or reconstructed candidates
    const bool has_injected = fSettings.IsMC && !fOutput.Injected.empty();

    if (has_rec_candidates || has_injected) fWriter->Fill();
    fOutput.Clear(fSettings.IsMC);

    // clear transient tracks //
    fTemp_AntiProton.clear();
    fTemp_Proton.clear();
    fTemp_NegKaon.clear();
    fTemp_PosKaon.clear();
    fTemp_PiMinus.clear();
    fTemp_PiPlus.clear();

    // clear transient v0s //
    fTemp_AntiLambda.clear();
    fTemp_AntiLambda_Neg.clear();
    fTemp_AntiLambda_Pos.clear();
    fTemp_Lambda.clear();
    fTemp_Lambda_Neg.clear();
    fTemp_Lambda_Pos.clear();
    fTemp_KaonZeroShort.clear();
    fTemp_KaonZeroShort_Neg.clear();
    fTemp_KaonZeroShort_Pos.clear();

    if (!fSettings.IsMC) return;

    // clear transient mc //
    fTemp_MC_AntiProton.clear();
    fTemp_MC_Proton.clear();
    fTemp_MC_NegKaon.clear();
    fTemp_MC_PosKaon.clear();
    fTemp_MC_PiMinus.clear();
    fTemp_MC_PiPlus.clear();
    fTemp_MC_AntiLambda.clear();
    fTemp_MC_AntiLambda_Neg.clear();
    fTemp_MC_AntiLambda_Pos.clear();
    fTemp_MC_Lambda.clear();
    fTemp_MC_Lambda_Neg.clear();
    fTemp_MC_Lambda_Pos.clear();
    fTemp_MC_KaonZeroShort.clear();
    fTemp_MC_KaonZeroShort_Neg.clear();
    fTemp_MC_KaonZeroShort_Pos.clear();
}

bool Finder::EndOfAnalysis() {

    Logger::Info(__FUNCTION__, "The following objects have been written into TFile \"{}\":", fSettings.PathOutputFile);

    Logger::Info(__FUNCTION__, "- RNTuple \"{}\"", T2DS::Name_FoundSexaquarkRNT);

    // write histograms //

    fOutput_File->cd();
    for (auto* hist : {
             fHist_EventCounter.get(),
             fHist_CutFlow_AntiProton.get(),
             fHist_CutFlow_Proton.get(),
             fHist_CutFlow_NegKaon.get(),
             fHist_CutFlow_PosKaon.get(),
             fHist_CutFlow_PiMinus.get(),
             fHist_CutFlow_PiPlus.get(),
             fHist_CutFlow_AntiLambda.get(),
             fHist_CutFlow_Lambda.get(),
             fHist_CutFlow_KaonZeroShort.get(),
             fHist_CutFlow_ChannelA.get(),
             fHist_CutFlow_ChannelA_Bkg.get(),
             fHist_CutFlow_ChannelD.get(),
             fHist_CutFlow_ChannelD_Bkg.get(),
             fHist_CutFlow_ChannelH.get(),
             fHist_CutFlow_ChannelH_Bkg.get(),
         }) {
        hist->Write();
        Logger::Info(__FUNCTION__, "- TH1D \"{}\"", hist->GetName());
    }

    Logger::Info(__FUNCTION__, "All done.");

    return fHist_EventCounter->GetEntries() != 0;
}

}  // namespace T2DS
