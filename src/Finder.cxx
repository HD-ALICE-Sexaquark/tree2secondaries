#include <cstddef>
#include <optional>
#include <string_view>

#include "common/DB_Particles.hpp"
#include "common/MC_Helpers.hpp"
#include "common/Math.hpp"
#include "common/POD_Sexaquark.hpp"
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
    fOutput_Base->Event = fInput.Event;
    // -- update event counter
    fHist_EventCounter->Fill(0.);
    // -- cache pv
    fPrimaryVertex.SetCoordinates(fOutput_Base->Event.PV_X, fOutput_Base->Event.PV_Y, fOutput_Base->Event.PV_Z);
}

// ## Injected ZONE ## //

void Finder::ProcessInjected() {
    // -- copy injected info
    fOutput_Base->Injected = std::move(fInput.InjectedSexa);
}

POD::Linked::InjectedSexa Finder::BuildMcSexaquark(const POD::Extended::McParticle& mc_dau1, const POD::Extended::McParticle& mc_dau2) {
    POD::Linked::InjectedSexa mc_sexa;
    // -- hybridness must be filled independently of found matching injected sexaquark
    mc_sexa.IsHybrid = mc_dau1.IsTrueSignal != mc_dau2.IsTrueSignal || mc_dau1.SignalID != mc_dau2.SignalID || mc_dau1.IsHybrid || mc_dau2.IsHybrid;
    // -- find common reaction id
    auto entry_inj = MC::SexaquarkRules::FindCommonReactionID(mc_dau1, mc_dau2);
    if (!entry_inj.has_value()) {
        return mc_sexa;
    }
    // -- if found, fill it with matching injected info
    static_cast<POD::Extended::InjectedSexa&>(mc_sexa) = fOutput_Base->Injected[entry_inj.value()];

    return mc_sexa;
}

// ## Channel A ZONE ## //

void Finder::FindSexaquarks_ChannelA(bool anti_channel) {

    // determine rules and aliases //
    // -- lambda
    //    -- rec
    const auto& input_lambdas = anti_channel ? fInput.Lambda : fInput.AntiLambda;
    const auto& input_lambdas_neg = anti_channel ? fInput.Lambda_Neg : fInput.AntiLambda_Neg;
    const auto& input_lambdas_pos = anti_channel ? fInput.Lambda_Pos : fInput.AntiLambda_Pos;
    const std::size_t n_lambdas = input_lambdas.size();
    //    -- mc
    const std::vector<POD::Extended::McParticle>* input_mc_lambdas = nullptr;
    const std::vector<POD::Extended::McParticle>* input_mc_lambdas_neg = nullptr;
    const std::vector<POD::Extended::McParticle>* input_mc_lambdas_pos = nullptr;
    if (fSettings.IsMC) {
        input_mc_lambdas = anti_channel ? &fInput.MC_Lambda : &fInput.MC_AntiLambda;
        input_mc_lambdas_neg = anti_channel ? &fInput.MC_Lambda_Neg : &fInput.MC_AntiLambda_Neg;
        input_mc_lambdas_pos = anti_channel ? &fInput.MC_Lambda_Pos : &fInput.MC_AntiLambda_Pos;
    }
    // -- k0s
    //    -- rec
    const auto& input_k0s = fInput.KaonZeroShort;
    const auto& input_k0s_neg = fInput.KaonZeroShort_Neg;
    const auto& input_k0s_pos = fInput.KaonZeroShort_Pos;
    const std::size_t n_k0s = input_k0s.size();
    //    -- mc
    const std::vector<POD::Extended::McParticle>* input_mc_k0s = nullptr;
    const std::vector<POD::Extended::McParticle>* input_mc_k0s_neg = nullptr;
    const std::vector<POD::Extended::McParticle>* input_mc_k0s_pos = nullptr;
    if (fSettings.IsMC) {
        input_mc_k0s = &fInput.MC_KaonZeroShort;
        input_mc_k0s_neg = &fInput.MC_KaonZeroShort_Neg;
        input_mc_k0s_pos = &fInput.MC_KaonZeroShort_Pos;
    }
    // -- cut flow hist
    TH1D* hist = anti_channel ? fHist_CutFlow_AntiChannel.get() : fHist_CutFlow.get();

    // loop over all possible pairs of (anti)lambda + K0S //
    for (std::size_t entry_lambda = 0; entry_lambda < n_lambdas; ++entry_lambda) {
        // cache index lookups //
        const POD::V0& lambda = input_lambdas[entry_lambda];
        const POD::Track& lambda_neg = input_lambdas_neg[entry_lambda];
        const POD::Track& lambda_pos = input_lambdas_pos[entry_lambda];

        /* lambda.CacheCalculations(entry_lambda, fPrimaryVertex); PENDING */

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
            if (!SlowCuts_ChannelA(kf_sexa, hist)) continue;

            // store reconstructed //
            fOutput_ChannelA.Sexaquark.emplace_back(Create_ChannelA(kf_sexa, anti_channel));
            fOutput_ChannelA.V0A.emplace_back(lambda);
            fOutput_ChannelA.V0A_Neg.emplace_back(lambda_neg);
            fOutput_ChannelA.V0A_Pos.emplace_back(lambda_pos);
            fOutput_ChannelA.V0B.emplace_back(k0s);
            fOutput_ChannelA.V0B_Neg.emplace_back(k0s_neg);
            fOutput_ChannelA.V0B_Pos.emplace_back(k0s_pos);

            // store mc //
            if (fSettings.IsMC) {
                // -- V0A
                const auto& mc_lambda = (*input_mc_lambdas)[entry_lambda];
                fOutput_ChannelA.MC_V0A.emplace_back(mc_lambda);
                fOutput_ChannelA.MC_V0A_Neg.emplace_back((*input_mc_lambdas_neg)[entry_lambda]);
                fOutput_ChannelA.MC_V0A_Pos.emplace_back((*input_mc_lambdas_pos)[entry_lambda]);
                // -- V0B
                const auto& mc_k0s = (*input_mc_k0s)[entry_k0s];
                fOutput_ChannelA.MC_V0B.emplace_back(mc_k0s);
                fOutput_ChannelA.MC_V0B_Neg.emplace_back((*input_mc_k0s_neg)[entry_k0s]);
                fOutput_ChannelA.MC_V0B_Pos.emplace_back((*input_mc_k0s_pos)[entry_k0s]);
                // -- h-dibaryon
                fOutput_ChannelA.MC_Sexaquark.emplace_back(BuildMcSexaquark(mc_lambda, mc_k0s));
            }
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

bool Finder::SlowCuts_ChannelA(const KF::ChannelA& kf_sexa, TH1D* cut_flow_hist) const {

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

POD::Sexaquark Finder::Create_ChannelA(const KF::ChannelA& kf_sexa, bool anti_channel) {
    POD::Sexaquark sexa;  // non-initialized on purpose
    sexa.SV_X = static_cast<float>(kf_sexa.X());
    sexa.SV_Y = static_cast<float>(kf_sexa.Y());
    sexa.SV_Z = static_cast<float>(kf_sexa.Z());
    sexa.Px = static_cast<float>(kf_sexa.Px());
    sexa.Py = static_cast<float>(kf_sexa.Py());
    sexa.Pz = static_cast<float>(kf_sexa.Pz());
    sexa.Energy = static_cast<float>(kf_sexa.E());
    sexa.Chi2NDF = static_cast<float>(kf_sexa.Chi2NDF());
    sexa.E_MinusNucleon = static_cast<float>(kf_sexa.E() - DB::Particles::Particle("Neutron").mass);  // small optimization
    sexa.AntiChannel = anti_channel;
    // -- V0A
    sexa.Dau1_PCAwrtSV_X = static_cast<float>(kf_sexa.V0A_PCAwrtSV.X());
    sexa.Dau1_PCAwrtSV_Y = static_cast<float>(kf_sexa.V0A_PCAwrtSV.Y());
    sexa.Dau1_PCAwrtSV_Z = static_cast<float>(kf_sexa.V0A_PCAwrtSV.Z());
    sexa.Dau1_PCAwrtSV_Px = static_cast<float>(kf_sexa.V0A_PCAwrtSV.Px());
    sexa.Dau1_PCAwrtSV_Py = static_cast<float>(kf_sexa.V0A_PCAwrtSV.Py());
    sexa.Dau1_PCAwrtSV_Pz = static_cast<float>(kf_sexa.V0A_PCAwrtSV.Pz());
    // -- V0B
    sexa.Dau2_PCAwrtSV_X = static_cast<float>(kf_sexa.V0B_PCAwrtSV.X());
    sexa.Dau2_PCAwrtSV_Y = static_cast<float>(kf_sexa.V0B_PCAwrtSV.Y());
    sexa.Dau2_PCAwrtSV_Z = static_cast<float>(kf_sexa.V0B_PCAwrtSV.Z());
    sexa.Dau2_PCAwrtSV_Px = static_cast<float>(kf_sexa.V0B_PCAwrtSV.Px());
    sexa.Dau2_PCAwrtSV_Py = static_cast<float>(kf_sexa.V0B_PCAwrtSV.Py());
    sexa.Dau2_PCAwrtSV_Pz = static_cast<float>(kf_sexa.V0B_PCAwrtSV.Pz());

    return sexa;
}

// ## Channel D ZONE ## //

void Finder::FindSexaquarks_ChannelD(bool anti_channel) {

    // determine rules and aliases //
    // -- lambda
    //    -- rec
    const auto& input_lambdas = anti_channel ? fInput.Lambda : fInput.AntiLambda;
    const auto& input_lambdas_neg = anti_channel ? fInput.Lambda_Neg : fInput.AntiLambda_Neg;
    const auto& input_lambdas_pos = anti_channel ? fInput.Lambda_Pos : fInput.AntiLambda_Pos;
    const std::size_t n_lambdas = input_lambdas.size();
    //    -- mc
    const std::vector<POD::Extended::McParticle>* input_mc_lambdas = nullptr;
    const std::vector<POD::Extended::McParticle>* input_mc_lambdas_neg = nullptr;
    const std::vector<POD::Extended::McParticle>* input_mc_lambdas_pos = nullptr;
    if (fSettings.IsMC) {
        input_mc_lambdas = anti_channel ? &fInput.MC_Lambda : &fInput.MC_AntiLambda;
        input_mc_lambdas_neg = anti_channel ? &fInput.MC_Lambda_Neg : &fInput.MC_AntiLambda_Neg;
        input_mc_lambdas_pos = anti_channel ? &fInput.MC_Lambda_Pos : &fInput.MC_AntiLambda_Pos;
    }
    // -- kaon
    //    -- rec
    const auto& input_kaons = anti_channel ? fInput.NegKaon : fInput.PosKaon;
    const std::size_t n_kaons = input_kaons.size();
    //    -- mc
    const std::vector<POD::Extended::McParticle>* input_mc_kaons = nullptr;
    if (fSettings.IsMC) {
        input_mc_kaons = anti_channel ? &fInput.MC_NegKaon : &fInput.MC_PosKaon;
    }
    constexpr double mass_kaon = DB::Particles::Particle("NegKaon").mass;
    // -- cut flow hist
    TH1D* hist = anti_channel ? fHist_CutFlow_AntiChannel.get() : fHist_CutFlow.get();

    // loop over all possible pairs of (anti)lambda + (pos/neg)kaon //
    for (std::size_t entry_lambda = 0; entry_lambda < n_lambdas; ++entry_lambda) {
        // cache index lookups //
        const POD::V0& lambda = input_lambdas[entry_lambda];
        const POD::Track& lambda_neg = input_lambdas_neg[entry_lambda];
        const POD::Track& lambda_pos = input_lambdas_pos[entry_lambda];

        /* v0.CacheCalculations(entry_v0, fPrimaryVertex); // PENDING */

        for (std::size_t entry_kaon = 0; entry_kaon < n_kaons; ++entry_kaon) {
            // cache index lookup //
            const POD::Track& kaon = input_kaons[entry_kaon];

            // -- sanity check
            if (lambda_neg.EsdEntry == kaon.EsdEntry || lambda_pos.EsdEntry == kaon.EsdEntry) continue;

            /* kaon.CacheCalculations(entry_kaon, fPrimaryVertex, fOutput.Event.MagneticField); // PENDING */

            // PCAs (1) //
            auto [seed_kaon, seed_v0, pca_cache] = Seeder::HelixLine::FastCorrectPCAs(kaon, lambda, fOutput_Base->Event.MagneticField);

            // apply cuts (1) //
            if (!FastCuts_ChannelD(seed_kaon, seed_v0, hist)) continue;

            // PCAs derivatives //
            auto [deriv_v0, deriv_ka] = Seeder::HelixLine::ComputeDerivatives(seed_kaon, seed_v0, pca_cache);

            // fit vertex //
            auto fit = KF::FitVertex(kaon, lambda, mass_kaon, {seed_kaon, deriv_ka}, {seed_v0, deriv_v0}, fOutput_Base->Event.MagneticField);
            KF::ChannelD kf_sexa(fit, seed_v0.pca, seed_kaon.pca, lambda, kaon);

            // apply cuts (2) //
            if (!SlowCuts_ChannelD(kf_sexa, hist)) continue;

            // store reconstructed //
            fOutput_ChannelD.Sexaquark.emplace_back(Create_ChannelD(kf_sexa, anti_channel));
            fOutput_ChannelD.V0.emplace_back(lambda);
            fOutput_ChannelD.V0_Neg.emplace_back(lambda_neg);
            fOutput_ChannelD.V0_Pos.emplace_back(lambda_pos);
            fOutput_ChannelD.Kaon.emplace_back(kaon);

            // store mc //
            if (fSettings.IsMC) {
                // -- V0
                const auto& mc_lambda = (*input_mc_lambdas)[entry_lambda];
                fOutput_ChannelD.MC_V0.emplace_back(mc_lambda);
                fOutput_ChannelD.MC_V0_Neg.emplace_back((*input_mc_lambdas_neg)[entry_lambda]);
                fOutput_ChannelD.MC_V0_Pos.emplace_back((*input_mc_lambdas_pos)[entry_lambda]);
                // -- Kaon
                const auto& mc_kaon = (*input_mc_kaons)[entry_kaon];
                fOutput_ChannelD.MC_Kaon.emplace_back(mc_kaon);
                // -- h-dibaryon
                fOutput_ChannelD.MC_Sexaquark.emplace_back(BuildMcSexaquark(mc_lambda, mc_kaon));
            }
        }
    }
}

bool Finder::FastCuts_ChannelD(const Seeder::Seed& seed_ka, const Seeder::Seed& seed_v0, TH1D* cut_flow_hist) const {

    cut_flow_hist->Fill(0.);

    if (Common::Math::SquaredDistance(seed_ka.pca.xyz, seed_v0.pca.xyz) > Cuts::ChannelD::Max_DCAKaLa * Cuts::ChannelD::Max_DCAKaLa) return false;
    cut_flow_hist->Fill(1.);

    return true;
}

bool Finder::SlowCuts_ChannelD(const KF::ChannelD& sexa, TH1D* cut_flow_hist) const {

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

POD::Sexaquark Finder::Create_ChannelD(const KF::ChannelD& kf_sexa, bool anti_channel) {
    POD::Sexaquark sexa;  // non-initialized on purpose
    sexa.SV_X = static_cast<float>(kf_sexa.X());
    sexa.SV_Y = static_cast<float>(kf_sexa.Y());
    sexa.SV_Z = static_cast<float>(kf_sexa.Z());
    sexa.Px = static_cast<float>(kf_sexa.Px());
    sexa.Py = static_cast<float>(kf_sexa.Py());
    sexa.Pz = static_cast<float>(kf_sexa.Pz());
    sexa.Energy = static_cast<float>(kf_sexa.E());
    sexa.Chi2NDF = static_cast<float>(kf_sexa.Chi2NDF());
    sexa.E_MinusNucleon = static_cast<float>(kf_sexa.E() - DB::Particles::Particle("Proton").mass);  // small optimization
    sexa.AntiChannel = anti_channel;
    // -- V0
    sexa.Dau1_PCAwrtSV_X = static_cast<float>(kf_sexa.V0_PCAwrtSV.X());
    sexa.Dau1_PCAwrtSV_Y = static_cast<float>(kf_sexa.V0_PCAwrtSV.Y());
    sexa.Dau1_PCAwrtSV_Z = static_cast<float>(kf_sexa.V0_PCAwrtSV.Z());
    sexa.Dau1_PCAwrtSV_Px = static_cast<float>(kf_sexa.V0_PCAwrtSV.Px());
    sexa.Dau1_PCAwrtSV_Py = static_cast<float>(kf_sexa.V0_PCAwrtSV.Py());
    sexa.Dau1_PCAwrtSV_Pz = static_cast<float>(kf_sexa.V0_PCAwrtSV.Pz());
    // -- kaon
    sexa.Dau2_PCAwrtSV_X = static_cast<float>(kf_sexa.Kaon_PCAwrtSV.X());
    sexa.Dau2_PCAwrtSV_Y = static_cast<float>(kf_sexa.Kaon_PCAwrtSV.Y());
    sexa.Dau2_PCAwrtSV_Z = static_cast<float>(kf_sexa.Kaon_PCAwrtSV.Z());
    sexa.Dau2_PCAwrtSV_Px = static_cast<float>(kf_sexa.Kaon_PCAwrtSV.Px());
    sexa.Dau2_PCAwrtSV_Py = static_cast<float>(kf_sexa.Kaon_PCAwrtSV.Py());
    sexa.Dau2_PCAwrtSV_Pz = static_cast<float>(kf_sexa.Kaon_PCAwrtSV.Pz());

    return sexa;
}

// ## Channel H ZONE ## //

void Finder::FindSexaquarks_ChannelH(bool anti_channel) {

    // determine rules and aliases //
    // -- kaons
    //    -- rec
    const auto& input_kaons = anti_channel ? fInput.NegKaon : fInput.PosKaon;
    const std::size_t n_kaons = input_kaons.size();
    //    -- mc
    const std::vector<POD::Extended::McParticle>* input_mc_kaons = nullptr;
    if (fSettings.IsMC) input_mc_kaons = anti_channel ? &fInput.MC_NegKaon : &fInput.MC_PosKaon;
    constexpr double mass_kaon = DB::Particles::Particle("NegKaon").mass;
    // -- cut flow hist
    TH1D* hist = anti_channel ? fHist_CutFlow_AntiChannel.get() : fHist_CutFlow.get();

    // loop over all possible pairs of (pos)kaon+(pos)kaon or (neg)kaon+(neg)kaon //
    for (std::size_t entry_kaon1 = 0; entry_kaon1 + 1 < n_kaons; ++entry_kaon1) {
        const POD::Track& kaon1 = input_kaons[entry_kaon1];  // cache index lookup

        /* kaon1.CacheCalculations(entry_kaon1, fPrimaryVertex, fOutput.Event.MagneticField); // PENDING */

        for (std::size_t entry_kaon2 = entry_kaon1 + 1; entry_kaon2 < n_kaons; ++entry_kaon2) {
            // NOTE: sanity check not needed, because loops don't intersect
            const POD::Track& kaon2 = input_kaons[entry_kaon2];  // cache index lookup

            /* kaon2.CacheCalculations(entry_kaon2, fPrimaryVertex, fOutput.Event.MagneticField); // PENDING */

            // PCAs (1) //
            auto [seed_kaon1, seed_kaon2, pca_cache] = Seeder::HelixHelix::FastCorrectPCAs(kaon1, kaon2, fOutput_Base->Event.MagneticField);

            // apply cuts (1) //
            if (!FastCuts_ChannelH(seed_kaon1, seed_kaon2, hist)) continue;

            // PCAs derivatives //
            auto [deriv_kaon1, deriv_kaon2] = Seeder::HelixHelix::ComputeDerivatives(seed_kaon1, seed_kaon2, pca_cache);

            // fit vertex //
            auto fit = KF::FitVertex(kaon1, kaon2, mass_kaon, mass_kaon, {seed_kaon1, deriv_kaon1}, {seed_kaon2, deriv_kaon2},
                                     fOutput_Base->Event.MagneticField);
            KF::ChannelH kf_sexa(fit, seed_kaon1.pca, seed_kaon2.pca, kaon1, kaon2);

            // apply cuts (2) //
            if (!SlowCuts_ChannelH(kf_sexa, hist)) continue;

            // store reconstructed //
            fOutput_ChannelH.Sexaquark.emplace_back(Create_ChannelH(kf_sexa, anti_channel));
            fOutput_ChannelH.Kaon1.emplace_back(kaon1);
            fOutput_ChannelH.Kaon2.emplace_back(kaon2);

            // store mc //
            if (fSettings.IsMC) {
                // -- Kaon1
                const auto& mc_kaon1 = (*input_mc_kaons)[entry_kaon1];
                fOutput_ChannelH.MC_Kaon1.emplace_back(mc_kaon1);
                // -- Kaon2
                const auto& mc_kaon2 = (*input_mc_kaons)[entry_kaon2];
                fOutput_ChannelH.MC_Kaon2.emplace_back(mc_kaon2);
                // -- h-dibaryon
                fOutput_ChannelH.MC_Sexaquark.emplace_back(BuildMcSexaquark(mc_kaon1, mc_kaon2));
            }
        }
    }
}

bool Finder::FastCuts_ChannelH(const Seeder::Seed& seed_kaon1, const Seeder::Seed& seed_kaon2, TH1D* cut_flow_hist) const {
    //
    // PENDING
    //
    return true;
}

bool Finder::SlowCuts_ChannelH(const KF::ChannelH& sexa, TH1D* cut_flow_hist) const {
    //
    // PENDING
    //
    return true;
}

POD::Sexaquark Finder::Create_ChannelH(const KF::ChannelH& kf_sexa, bool anti_channel) {
    POD::Sexaquark sexa;  // non-initialized on purpose
    sexa.SV_X = static_cast<float>(kf_sexa.X());
    sexa.SV_Y = static_cast<float>(kf_sexa.Y());
    sexa.SV_Z = static_cast<float>(kf_sexa.Z());
    sexa.Px = static_cast<float>(kf_sexa.Px());
    sexa.Py = static_cast<float>(kf_sexa.Py());
    sexa.Pz = static_cast<float>(kf_sexa.Pz());
    sexa.Energy = static_cast<float>(kf_sexa.E());
    sexa.Chi2NDF = static_cast<float>(kf_sexa.Chi2NDF());
    sexa.E_MinusNucleon = static_cast<float>(kf_sexa.E() - DB::Particles::Particle("Proton").mass);  // small optimization
    sexa.AntiChannel = anti_channel;
    // -- Kaon 1
    sexa.Dau1_PCAwrtSV_X = static_cast<float>(kf_sexa.Kaon1_PCAwrtSV.X());
    sexa.Dau1_PCAwrtSV_Y = static_cast<float>(kf_sexa.Kaon1_PCAwrtSV.Y());
    sexa.Dau1_PCAwrtSV_Z = static_cast<float>(kf_sexa.Kaon1_PCAwrtSV.Z());
    sexa.Dau1_PCAwrtSV_Px = static_cast<float>(kf_sexa.Kaon1_PCAwrtSV.Px());
    sexa.Dau1_PCAwrtSV_Py = static_cast<float>(kf_sexa.Kaon1_PCAwrtSV.Py());
    sexa.Dau1_PCAwrtSV_Pz = static_cast<float>(kf_sexa.Kaon1_PCAwrtSV.Pz());
    // -- Kaon 2
    sexa.Dau2_PCAwrtSV_X = static_cast<float>(kf_sexa.Kaon2_PCAwrtSV.X());
    sexa.Dau2_PCAwrtSV_Y = static_cast<float>(kf_sexa.Kaon2_PCAwrtSV.Y());
    sexa.Dau2_PCAwrtSV_Z = static_cast<float>(kf_sexa.Kaon2_PCAwrtSV.Z());
    sexa.Dau2_PCAwrtSV_Px = static_cast<float>(kf_sexa.Kaon2_PCAwrtSV.Px());
    sexa.Dau2_PCAwrtSV_Py = static_cast<float>(kf_sexa.Kaon2_PCAwrtSV.Py());
    sexa.Dau2_PCAwrtSV_Pz = static_cast<float>(kf_sexa.Kaon2_PCAwrtSV.Pz());

    return sexa;
}

void Finder::EndOfEvent() {
    // if data, don't keep event with no candidates
    // if mc, keep event with injected or reconstructed candidates
    // NOTE: if they're empty, there's nothing to clear
    if (fReactionChannel.name == 'A') {
        if (!fSettings.IsMC && fOutput_ChannelA.Sexaquark.empty()) return;
        if (fSettings.IsMC && fOutput_ChannelA.Sexaquark.empty() && fOutput_ChannelA.Injected.empty()) return;
    }
    if (fReactionChannel.name == 'D') {
        if (!fSettings.IsMC && fOutput_ChannelD.Sexaquark.empty()) return;
        if (fSettings.IsMC && fOutput_ChannelD.Sexaquark.empty() && fOutput_ChannelD.Injected.empty()) return;
    }
    if (fReactionChannel.name == 'H') {
        if (!fSettings.IsMC && fOutput_ChannelH.Sexaquark.empty()) return;
        if (fSettings.IsMC && fOutput_ChannelH.Sexaquark.empty() && fOutput_ChannelH.Injected.empty()) return;
    }
    // fill schema
    fWriter->Fill();
    // clear schema vectors
    if (fReactionChannel.name == 'A')
        fOutput_ChannelA.Clear(fSettings.IsMC);
    else if (fReactionChannel.name == 'D')
        fOutput_ChannelD.Clear(fSettings.IsMC);
    else if (fReactionChannel.name == 'H')
        fOutput_ChannelH.Clear(fSettings.IsMC);
}

// ## END OF CYCLES ## //

void Finder::EndOfAnalysis() {

    Logger::Info(__FUNCTION__, "The following objects have been written into TFile {}:", fSettings.PathOutputFile);

    Logger::Info(__FUNCTION__, "- RNTuple \"{}\"", fName_FoundRNT);

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
