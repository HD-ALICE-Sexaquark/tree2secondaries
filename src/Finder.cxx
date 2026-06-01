#include <cstddef>
#include <memory>
#include <optional>
#include <string_view>

#include "common/Constants.hpp"
#include "common/DB_Particles.hpp"
#include "common/MC_Helpers.hpp"
#include "common/Math.hpp"
#include "common/POD_ChannelA.hpp"
#include "common/POD_Track.hpp"
#include "common/R2DS_Cuts.hpp"

#include "KalmanFitter/BaseKalmanFitter.hxx"
#include "KalmanFitter/KalmanFitterChannelA.hxx"
#include "KalmanFitter/KalmanFitterChannelD.hxx"
#include "KalmanFitter/KalmanFitterChannelH.hxx"
#include "Seeder/SeederHelixHelix.hxx"
#include "Seeder/SeederHelixLine.hxx"
#include "Seeder/SeederLineLine.hxx"

#include "Finder/Finder.hxx"

namespace R2DS {

// ## OUTPUT ZONE ## //

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

void Finder::ProcessEvent() {
    // -- copy event info
    fOutput.Event = fInput.Event;
    // -- update event counter
    fHist_EventCounter->Fill(0.);
    // -- cache pv
    fPrimaryVertex.SetCoordinates(fOutput.Event.PV_X, fOutput.Event.PV_Y, fOutput.Event.PV_Z);
}

// ## Injected ZONE ## //

void Finder::ProcessInjected() {
    // data kind must be MC //
    if (!fInput.InjectedSexa.has_value()) return;
    // -- copy injected info
    fOutput.InjectedSexa = std::move(fInput.InjectedSexa);
}

// ## Channel A ZONE ## //

void Finder::FindSexaquarks_ChannelA(bool control_channel) {

    // determine rules and aliases //
    // -- lambda
    const auto& input_lambdas = control_channel ? fInput.Lambda : fInput.AntiLambda;
    const std::size_t n_lambdas = input_lambdas.size();
    // -- k0s
    const auto& input_k0s = fInput.KaonZeroShort;
    const std::size_t n_k0s = input_k0s.size();
    // -- cut flow hist
    TH1D* hist = control_channel ? fHist_CutFlow_AntiChannel.get() : fHist_CutFlow.get();

    // loop over all possible pairs of (anti)lambda + K0S //
    for (std::size_t entry_lambda = 0; entry_lambda < n_lambdas; ++entry_lambda) {
        const POD::V0& lambda = input_lambdas[entry_lambda];  // cache index lookup
        /* lambda.CacheCalculations(entry_lambda, fPrimaryVertex); PENDING */

        for (std::size_t entry_k0s = 0; entry_k0s < n_k0s; ++entry_k0s) {
            const POD::V0& k0s = input_k0s[entry_k0s];  // cache index lookup

            // -- sanity check
            if (lambda.Neg.EsdEntry == k0s.Neg.EsdEntry || lambda.Neg.EsdEntry == k0s.Pos.EsdEntry || lambda.Pos.EsdEntry == k0s.Neg.EsdEntry ||
                lambda.Pos.EsdEntry == k0s.Pos.EsdEntry) {
                continue;
            }
            /* k0s.CacheCalculations(entry_k0s, fPrimaryVertex); PENDING */

            // PCAs //
            Seeder::LineLine::Cache pca_cache;
            auto [seed_lambda, seed_k0s] = Seeder::LineLine::FastPCAs(lambda, k0s, &pca_cache);

            // apply cuts (1) //
            if (!FastCuts_ChannelA(seed_lambda, seed_k0s, hist)) continue;

            // PCAs derivatives //
            auto [deriv_lambda, deriv_k0s] = Seeder::LineLine::ComputeDerivatives(pca_cache);

            // fit vertex //
            auto fit = KF::FitVertex(lambda, k0s, {seed_lambda, deriv_lambda}, {seed_k0s, deriv_k0s});
            KF::ChannelA kf_sexa(fit, seed_lambda.pca, seed_k0s.pca, lambda, k0s);

            // apply cuts (2) //
            if (!SlowCuts(kf_sexa, hist)) continue;

            // store //
            POD::ChannelA new_sexa;
            // -- reconstructed info
            BuildRecInfo(new_sexa, kf_sexa, control_channel);
            // -- linked mc info
            BuildMcInfo(new_sexa, control_channel);
            // fill //
            fOutput.ChannelA.emplace_back(new_sexa);
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

bool Finder::SlowCuts(const KF::ChannelA& kf_sexa, TH1D* cut_flow_hist) const {

    if (kf_sexa.SquaredRadius2D() < Cuts::ChannelA::Min_Radius2D * Cuts::ChannelA::Min_Radius2D) return false;
    cut_flow_hist->Fill(2.);

    if (kf_sexa.SquaredDCA_V0A_SV() > Cuts::ChannelA::Max_DCALaSV * Cuts::ChannelA::Max_DCALaSV) return false;
    cut_flow_hist->Fill(3.);

    if (kf_sexa.SquaredDCA_V0B_SV() > Cuts::ChannelA::Max_DCAK0SV * Cuts::ChannelA::Max_DCAK0SV) return false;
    cut_flow_hist->Fill(4.);

    // if (kf_sexa.CPA_wrt(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) < Cuts::ChannelA::Min_CPAwrtPV) return false;
    // cut_flow_hist->Fill(5.);

    // if (kf_sexa.DCA_wrt(fInput_Event.PV.X, fInput_Event.PV.Y, fInput_Event.PV.Z) > Cuts::ChannelA::Max_DCAwrtPV) return false;
    // cut_flow_hist->Fill(6.);

    // if (kf_sexa.CPA_V0A_SV() < Cuts::ChannelA::Min_La_CPAwrtSV) return false;
    // cut_flow_hist->Fill(7.);

    // if (kf_sexa.CPA_V0B_SV() < Cuts::ChannelA::Min_K0S_CPAwrtSV) return false;
    // cut_flow_hist->Fill(8.);

    return true;
}

void Finder::BuildMcInfo(POD::ChannelA& new_sexa, bool control_channel) {
    // data kind must be MC (already filled in `ProcessInjected(...)` //
    if (!fOutput.InjectedSexa.has_value()) return;
    new_sexa.MC = MC::SexaquarkRules::FindInjected(new_sexa.V0A.MC, new_sexa.V0B.MC, fOutput.InjectedSexa.value(), control_channel);
    new_sexa.IsSignal = new_sexa.MC.has_value();
    // clang-format off
    if (new_sexa.IsSignal) { new_sexa.IsHybrid = false; return; }
    // clang-format on
    auto [sig_a, hyb_a] = MC::IsSignalIsHybrid(new_sexa.V0A.MC);
    auto [sig_b, hyb_b] = MC::IsSignalIsHybrid(new_sexa.V0B.MC);
    new_sexa.IsHybrid = hyb_a || hyb_b || (sig_a != sig_b);
}

void Finder::BuildRecInfo(POD::ChannelA& new_sexa, const KF::ChannelA& kf_sexa, bool control_channel) {
    // candidate info
    new_sexa.SV_X = static_cast<float>(kf_sexa.X());
    new_sexa.SV_Y = static_cast<float>(kf_sexa.Y());
    new_sexa.SV_Z = static_cast<float>(kf_sexa.Z());
    new_sexa.Px = static_cast<float>(kf_sexa.Px());
    new_sexa.Py = static_cast<float>(kf_sexa.Py());
    new_sexa.Pz = static_cast<float>(kf_sexa.Pz());
    new_sexa.Energy = static_cast<float>(kf_sexa.E());
    new_sexa.Chi2NDF = static_cast<float>(kf_sexa.Chi2NDF());
    new_sexa.E_MinusNucleon = static_cast<float>(kf_sexa.E() - DB::Particles::Particle("Neutron").mass);
    new_sexa.ControlChannel = control_channel;
    // -- V0A
    new_sexa.V0A = *kf_sexa.V0A;
    new_sexa.V0A_PCAwrtSV_X = static_cast<float>(kf_sexa.V0A_PCAwrtSV.X());
    new_sexa.V0A_PCAwrtSV_Y = static_cast<float>(kf_sexa.V0A_PCAwrtSV.Y());
    new_sexa.V0A_PCAwrtSV_Z = static_cast<float>(kf_sexa.V0A_PCAwrtSV.Z());
    new_sexa.V0A_PCAwrtSV_Px = static_cast<float>(kf_sexa.V0A_PCAwrtSV.Px());
    new_sexa.V0A_PCAwrtSV_Py = static_cast<float>(kf_sexa.V0A_PCAwrtSV.Py());
    new_sexa.V0A_PCAwrtSV_Pz = static_cast<float>(kf_sexa.V0A_PCAwrtSV.Pz());
    // -- V0B
    new_sexa.V0B = *kf_sexa.V0B;
    new_sexa.V0B_PCAwrtSV_X = static_cast<float>(kf_sexa.V0B_PCAwrtSV.X());
    new_sexa.V0B_PCAwrtSV_Y = static_cast<float>(kf_sexa.V0B_PCAwrtSV.Y());
    new_sexa.V0B_PCAwrtSV_Z = static_cast<float>(kf_sexa.V0B_PCAwrtSV.Z());
    new_sexa.V0B_PCAwrtSV_Px = static_cast<float>(kf_sexa.V0B_PCAwrtSV.Px());
    new_sexa.V0B_PCAwrtSV_Py = static_cast<float>(kf_sexa.V0B_PCAwrtSV.Py());
    new_sexa.V0B_PCAwrtSV_Pz = static_cast<float>(kf_sexa.V0B_PCAwrtSV.Pz());
}

// ## Channel D ZONE ## //

void Finder::FindSexaquarks_ChannelD(bool control_channel) {

    // determine rules and aliases //
    // -- lambda
    const auto& input_lambdas = control_channel ? fInput.Lambda : fInput.AntiLambda;
    const std::size_t n_lambdas = input_lambdas.size();
    // -- kaons
    const auto& input_kaons = control_channel ? fInput.PosKaon : fInput.NegKaon;
    const std::size_t n_kaons = input_kaons.size();
    constexpr double mass_kaon = DB::Particles::Particle("PosKaon").mass;
    // -- cut flow hist
    TH1D* hist = control_channel ? fHist_CutFlow_AntiChannel.get() : fHist_CutFlow.get();

    // loop over all possible pairs of (anti)lambda + (pos/neg)kaon //
    for (std::size_t entry_lambda = 0; entry_lambda < n_lambdas; ++entry_lambda) {
        const POD::V0& lambda = input_lambdas[entry_lambda];  // cache index lookup

        /* v0.CacheCalculations(entry_v0, fPrimaryVertex); // PENDING */

        for (std::size_t entry_kaon = 0; entry_kaon < n_kaons; ++entry_kaon) {
            const POD::Track& kaon = input_kaons[entry_kaon];  // cache index lookup

            // -- sanity check
            if (lambda.Neg.EsdEntry == kaon.EsdEntry || lambda.Pos.EsdEntry == kaon.EsdEntry) continue;

            /* kaon.CacheCalculations(entry_kaon, fPrimaryVertex, fOutput.Event.MagneticField); // PENDING */

            // PCAs (1) //
            auto [seed_kaon, seed_v0, pca_cache] = Seeder::HelixLine::FastCorrectPCAs(kaon, lambda, fOutput.Event.MagneticField);

            // apply cuts (1) //
            if (!FastCuts_ChannelD(seed_kaon, seed_v0, hist)) continue;

            // PCAs derivatives //
            auto [deriv_v0, deriv_ka] = Seeder::HelixLine::ComputeDerivatives(seed_kaon, seed_v0, pca_cache);

            // fit vertex //
            auto fit = KF::FitVertex(kaon, lambda, mass_kaon,                     //
                                     {seed_kaon, deriv_ka}, {seed_v0, deriv_v0},  //
                                     fOutput.Event.MagneticField);
            KF::ChannelD kf_sexa(fit, seed_v0.pca, seed_kaon.pca, lambda, kaon);

            // apply cuts (2) //
            if (!SlowCuts(kf_sexa, hist)) continue;

            // store //
            POD::ChannelD new_sexa;
            // -- reconstructed info
            BuildRecInfo(new_sexa, kf_sexa, control_channel);
            // -- linked mc info
            BuildMcInfo(new_sexa, control_channel);
            // fill //
            fOutput.ChannelD.emplace_back(new_sexa);
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

void Finder::BuildMcInfo(POD::ChannelD& new_sexa, bool control_channel) {
    // data kind must be MC (already filled in `ProcessInjected(...)` //
    if (!fOutput.InjectedSexa.has_value()) return;
    new_sexa.MC = MC::SexaquarkRules::FindInjected(new_sexa.V0.MC, new_sexa.Kaon.MC, fOutput.InjectedSexa.value(), control_channel);
    new_sexa.IsSignal = new_sexa.MC.has_value();
    // clang-format off
    if (new_sexa.IsSignal) { new_sexa.IsHybrid = false; return; }
    // clang-format on
    auto [sig_a, hyb_a] = MC::IsSignalIsHybrid(new_sexa.V0.MC);
    auto [sig_b, hyb_b] = MC::IsSignalIsHybrid(new_sexa.Kaon.MC);
    new_sexa.IsHybrid = hyb_a || hyb_b || (sig_a != sig_b);
}

void Finder::BuildRecInfo(POD::ChannelD& new_sexa, const KF::ChannelD& kf_sexa, bool control_channel) {
    // candidate info
    new_sexa.SV_X = static_cast<float>(kf_sexa.X());
    new_sexa.SV_Y = static_cast<float>(kf_sexa.Y());
    new_sexa.SV_Z = static_cast<float>(kf_sexa.Z());
    new_sexa.Px = static_cast<float>(kf_sexa.Px());
    new_sexa.Py = static_cast<float>(kf_sexa.Py());
    new_sexa.Pz = static_cast<float>(kf_sexa.Pz());
    new_sexa.Energy = static_cast<float>(kf_sexa.E());
    new_sexa.Chi2NDF = static_cast<float>(kf_sexa.Chi2NDF());
    new_sexa.E_MinusNucleon = static_cast<float>(kf_sexa.E() - DB::Particles::Particle("Proton").mass);
    new_sexa.ControlChannel = control_channel;
    // -- V0
    new_sexa.V0 = *kf_sexa.V0;
    new_sexa.V0_PCAwrtSV_X = static_cast<float>(kf_sexa.V0_PCAwrtSV.X());
    new_sexa.V0_PCAwrtSV_Y = static_cast<float>(kf_sexa.V0_PCAwrtSV.Y());
    new_sexa.V0_PCAwrtSV_Z = static_cast<float>(kf_sexa.V0_PCAwrtSV.Z());
    new_sexa.V0_PCAwrtSV_Px = static_cast<float>(kf_sexa.V0_PCAwrtSV.Px());
    new_sexa.V0_PCAwrtSV_Py = static_cast<float>(kf_sexa.V0_PCAwrtSV.Py());
    new_sexa.V0_PCAwrtSV_Pz = static_cast<float>(kf_sexa.V0_PCAwrtSV.Pz());
    // -- kaon
    new_sexa.Kaon = *kf_sexa.Kaon;
    new_sexa.Kaon_PCAwrtSV_X = static_cast<float>(kf_sexa.Kaon_PCAwrtSV.X());
    new_sexa.Kaon_PCAwrtSV_Y = static_cast<float>(kf_sexa.Kaon_PCAwrtSV.Y());
    new_sexa.Kaon_PCAwrtSV_Z = static_cast<float>(kf_sexa.Kaon_PCAwrtSV.Z());
    new_sexa.Kaon_PCAwrtSV_Px = static_cast<float>(kf_sexa.Kaon_PCAwrtSV.Px());
    new_sexa.Kaon_PCAwrtSV_Py = static_cast<float>(kf_sexa.Kaon_PCAwrtSV.Py());
    new_sexa.Kaon_PCAwrtSV_Pz = static_cast<float>(kf_sexa.Kaon_PCAwrtSV.Pz());
}

// ## Channel H ZONE ## //

void Finder::FindSexaquarks_ChannelH(bool control_channel) {

    // determine rules and aliases //
    // -- kaons
    const auto& input_kaons = control_channel ? fInput.PosKaon : fInput.NegKaon;
    const std::size_t n_kaons = input_kaons.size();
    // -- cut flow hist
    TH1D* hist = control_channel ? fHist_CutFlow_AntiChannel.get() : fHist_CutFlow.get();
    constexpr double mass_kaon = DB::Particles::Particle("PosKaon").mass;

    // loop over all possible pairs of (pos)kaon+(pos)kaon or (neg)kaon+(neg)kaon //
    for (std::size_t entry_kaon1 = 0; entry_kaon1 + 1 < n_kaons; ++entry_kaon1) {
        const POD::Track& kaon1 = input_kaons[entry_kaon1];  // cache index lookup

        /* kaon1.CacheCalculations(entry_kaon1, fPrimaryVertex, fOutput.Event.MagneticField); // PENDING */

        for (std::size_t entry_kaon2 = entry_kaon1 + 1; entry_kaon2 < n_kaons; ++entry_kaon2) {
            // NOTE: sanity check not needed, because loops don't intersect
            const POD::Track& kaon2 = input_kaons[entry_kaon2];  // cache index lookup

            /* kaon2.CacheCalculations(entry_kaon2, fPrimaryVertex, fOutput.Event.MagneticField); PENDING */

            // PCAs (1) //
            auto [seed_kaon1, seed_kaon2, pca_cache] = Seeder::HelixHelix::FastCorrectPCAs(kaon1, kaon2, fOutput.Event.MagneticField);

            // apply cuts (1) //
            if (!FastCuts_ChannelH(seed_kaon1, seed_kaon2, hist)) continue;

            // PCAs derivatives //
            auto [deriv_kaon1, deriv_kaon2] = Seeder::HelixHelix::ComputeDerivatives(seed_kaon1, seed_kaon2, pca_cache);

            // fit vertex //
            auto fit = KF::FitVertex(kaon1, kaon2, mass_kaon, mass_kaon,  //
                                     {seed_kaon1, deriv_kaon1}, {seed_kaon2, deriv_kaon2}, fOutput.Event.MagneticField);
            KF::ChannelH kf_sexa(fit, seed_kaon1.pca, seed_kaon2.pca, kaon1, kaon2);

            // apply cuts (2) //
            if (!SlowCuts(kf_sexa, hist)) continue;

            // store //
            POD::ChannelH new_sexa;
            // -- reconstructed info
            BuildRecInfo(new_sexa, kf_sexa, control_channel);
            // -- linked mc info
            BuildMcInfo(new_sexa, control_channel);
            // fill //
            fOutput.ChannelH.emplace_back(new_sexa);
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

void Finder::BuildMcInfo(POD::ChannelH& new_sexa, bool control_channel) {
    // data kind must be MC (already filled in `ProcessInjected(...)` //
    if (!fOutput.InjectedSexa.has_value()) return;
    new_sexa.MC = MC::SexaquarkRules::FindInjected(new_sexa.Kaon1.MC, new_sexa.Kaon2.MC, fOutput.InjectedSexa.value(), control_channel);
    new_sexa.IsSignal = new_sexa.MC.has_value();
    // clang-format off
    if (new_sexa.IsSignal) { new_sexa.IsHybrid = false; return; }
    // clang-format on
    auto [sig_a, hyb_a] = MC::IsSignalIsHybrid(new_sexa.Kaon1.MC);
    auto [sig_b, hyb_b] = MC::IsSignalIsHybrid(new_sexa.Kaon2.MC);
    new_sexa.IsHybrid = hyb_a || hyb_b || (sig_a != sig_b);
}

void Finder::BuildRecInfo(POD::ChannelH& new_sexa, const KF::ChannelH& kf_sexa, bool control_channel) {
    // candidate info
    new_sexa.SV_X = static_cast<float>(kf_sexa.X());
    new_sexa.SV_Y = static_cast<float>(kf_sexa.Y());
    new_sexa.SV_Z = static_cast<float>(kf_sexa.Z());
    new_sexa.Px = static_cast<float>(kf_sexa.Px());
    new_sexa.Py = static_cast<float>(kf_sexa.Py());
    new_sexa.Pz = static_cast<float>(kf_sexa.Pz());
    new_sexa.Energy = static_cast<float>(kf_sexa.E());
    new_sexa.Chi2NDF = static_cast<float>(kf_sexa.Chi2NDF());
    new_sexa.E_MinusNucleon = static_cast<float>(kf_sexa.E() - DB::Particles::Particle("Proton").mass);
    new_sexa.ControlChannel = control_channel;
    // -- kaon 1
    new_sexa.Kaon1 = *kf_sexa.Kaon1;
    new_sexa.Kaon1_PCAwrtSV_X = static_cast<float>(kf_sexa.Kaon1_PCAwrtSV.X());
    new_sexa.Kaon1_PCAwrtSV_Y = static_cast<float>(kf_sexa.Kaon1_PCAwrtSV.Y());
    new_sexa.Kaon1_PCAwrtSV_Z = static_cast<float>(kf_sexa.Kaon1_PCAwrtSV.Z());
    new_sexa.Kaon1_PCAwrtSV_Px = static_cast<float>(kf_sexa.Kaon1_PCAwrtSV.Px());
    new_sexa.Kaon1_PCAwrtSV_Py = static_cast<float>(kf_sexa.Kaon1_PCAwrtSV.Py());
    new_sexa.Kaon1_PCAwrtSV_Pz = static_cast<float>(kf_sexa.Kaon1_PCAwrtSV.Pz());
    // -- kaon 2
    new_sexa.Kaon2 = *kf_sexa.Kaon2;
    new_sexa.Kaon2_PCAwrtSV_X = static_cast<float>(kf_sexa.Kaon2_PCAwrtSV.X());
    new_sexa.Kaon2_PCAwrtSV_Y = static_cast<float>(kf_sexa.Kaon2_PCAwrtSV.Y());
    new_sexa.Kaon2_PCAwrtSV_Z = static_cast<float>(kf_sexa.Kaon2_PCAwrtSV.Z());
    new_sexa.Kaon2_PCAwrtSV_Px = static_cast<float>(kf_sexa.Kaon2_PCAwrtSV.Px());
    new_sexa.Kaon2_PCAwrtSV_Py = static_cast<float>(kf_sexa.Kaon2_PCAwrtSV.Py());
    new_sexa.Kaon2_PCAwrtSV_Pz = static_cast<float>(kf_sexa.Kaon2_PCAwrtSV.Pz());
}

void Finder::EndOfEvent() {
    // fill schema
    fWriter.Fill();
    // clear schema vectors
    fOutput.Clear();
}

// ## END OF CYCLES ## //

void Finder::EndOfAnalysis() {

    Logger::Info(__FUNCTION__, "The following objects have been written into TFile {}:", fSettings.PathOutputFile);

    Logger::Info(__FUNCTION__, "- RNTuple \"{}\"", R2DS::Name_FoundSexaRNT);

    // write histograms //

    fOutput_File->cd();

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
