#pragma once

#include <array>
#include <cstddef>
#include <optional>
#include <span>
#include <string_view>

#include "common/Cached_ChannelA.hpp"
#include "common/Cached_ChannelD.hpp"
#include "common/Cached_ChannelH.hpp"
#include "common/Cached_Hdibaryon.hpp"

namespace Skimmer::Variables {

// ## The Variable Registry ## //

template <typename C>
struct Definition {
    using Extractor = double (*)(const C&);  // plain function pointer

    std::string_view Name;
    std::string_view Unit;
    Extractor Extract;
};

template <typename C>
[[nodiscard]] constexpr std::optional<std::size_t> FindIndex(std::span<const Definition<C>> db, std::string_view name) {
    for (std::size_t i = 0; i < db.size(); ++i)
        if (db[i].Name == name) return i;
    return std::nullopt;
}

// == Channel A : antisexaquark + neutron -> antilambda + K0S == //

inline constexpr auto DB_ChannelA = std::to_array<Definition<Cached::ChannelA>>({
    // -- kinematics
    {"Pt", "GeV/c", [](const Cached::ChannelA& c) { return c.Pt(); }},
    {"P", "GeV/c", [](const Cached::ChannelA& c) { return c.P(); }},
    {"Mass", "GeV/c2", [](const Cached::ChannelA& c) { return c.Mass(); }},
    {"Eta", "", [](const Cached::ChannelA& c) { return c.Eta(); }},
    {"Rapidity", "", [](const Cached::ChannelA& c) { return c.Rapidity(); }},
    {"Mass_MinusNucleon", "GeV/c2", [](const Cached::ChannelA& c) { return c.Mass_MinusNucleon(); }},
    {"Rapidity_MinusNucleon", "", [](const Cached::ChannelA& c) { return c.Rapidity_MinusNucleon(); }},
    // -- topology, w.r.t. the secondary vertex
    {"SV_Radius2D", "cm", [](const Cached::ChannelA& c) { return c.SV_Radius2D(); }},
    {"SV_Radius3D", "cm", [](const Cached::ChannelA& c) { return c.SV_Radius3D(); }},
    {"FlightLength", "cm", [](const Cached::ChannelA& c) { return c.FlightLength(); }},
    {"DCA_btw_Daughters", "cm", [](const Cached::ChannelA& c) { return c.DCA_btw_Daughters(); }},
    {"CosTheta_btw_Daughters", "", [](const Cached::ChannelA& c) { return c.CosTheta_btw_Daughters(); }},
    {"Theta_btw_Daughters", "rad", [](const Cached::ChannelA& c) { return c.Theta_btw_Daughters(); }},
    {"Chi2NDF", "", [](const Cached::ChannelA& c) { return static_cast<double>(c.Chi2NDF); }},
    // -- topology, w.r.t. the primary vertex
    {"Qt_wrt_PV", "GeV/c", [](const Cached::ChannelA& c) { return c.Qt_wrt_PV(); }},
    {"Ql_wrt_PV", "GeV/c", [](const Cached::ChannelA& c) { return c.Ql_wrt_PV(); }},
    {"DCA_wrt_PV", "cm", [](const Cached::ChannelA& c) { return c.DCA_wrt_PV(); }},
    {"CPA_wrt_PV", "", [](const Cached::ChannelA& c) { return c.CPA_wrt_PV(); }},
    // -- the antilambda leg
    {"V0A_DCA_wrt_SV", "cm", [](const Cached::ChannelA& c) { return c.V0A_DCA_wrt_SV(); }},
    {"V0A_CPA_wrt_SV", "", [](const Cached::ChannelA& c) { return c.V0A_CPA_wrt_SV(); }},
    {"V0A_DecayLength", "cm", [](const Cached::ChannelA& c) { return c.V0A_DecayLength(); }},
    {"V0A_Neg_Pt", "GeV/c", [](const Cached::ChannelA& c) { return c.V0A_Neg_Pt(); }},
    {"V0A_Neg_Eta", "", [](const Cached::ChannelA& c) { return c.V0A_Neg_Eta(); }},
    {"V0A_Pos_Pt", "GeV/c", [](const Cached::ChannelA& c) { return c.V0A_Pos_Pt(); }},
    {"V0A_Pos_Eta", "", [](const Cached::ChannelA& c) { return c.V0A_Pos_Eta(); }},
    // -- the K0S leg
    {"V0B_DCA_wrt_SV", "cm", [](const Cached::ChannelA& c) { return c.V0B_DCA_wrt_SV(); }},
    {"V0B_CPA_wrt_SV", "", [](const Cached::ChannelA& c) { return c.V0B_CPA_wrt_SV(); }},
    {"V0B_DecayLength", "cm", [](const Cached::ChannelA& c) { return c.V0B_DecayLength(); }},
    {"V0B_Neg_Pt", "GeV/c", [](const Cached::ChannelA& c) { return c.V0B_Neg_Pt(); }},
    {"V0B_Neg_Eta", "", [](const Cached::ChannelA& c) { return c.V0B_Neg_Eta(); }},
    {"V0B_Pos_Pt", "GeV/c", [](const Cached::ChannelA& c) { return c.V0B_Pos_Pt(); }},
    {"V0B_Pos_Eta", "", [](const Cached::ChannelA& c) { return c.V0B_Pos_Eta(); }},
});

// == Channel D : antisexaquark + proton -> antilambda + K+ == //

inline constexpr auto DB_ChannelD = std::to_array<Definition<Cached::ChannelD>>({
    // -- kinematics
    {"Pt", "GeV/c", [](const Cached::ChannelD& c) { return c.Pt(); }},
    {"P", "GeV/c", [](const Cached::ChannelD& c) { return c.P(); }},
    {"Mass", "GeV/c2", [](const Cached::ChannelD& c) { return c.Mass(); }},
    {"Eta", "", [](const Cached::ChannelD& c) { return c.Eta(); }},
    {"Rapidity", "", [](const Cached::ChannelD& c) { return c.Rapidity(); }},
    {"Mass_MinusNucleon", "GeV/c2", [](const Cached::ChannelD& c) { return c.Mass_MinusNucleon(); }},
    {"Rapidity_MinusNucleon", "", [](const Cached::ChannelD& c) { return c.Rapidity_MinusNucleon(); }},
    // -- topology, w.r.t. the secondary vertex
    {"SV_Radius2D", "cm", [](const Cached::ChannelD& c) { return c.SV_Radius2D(); }},
    {"SV_Radius3D", "cm", [](const Cached::ChannelD& c) { return c.SV_Radius3D(); }},
    {"FlightLength", "cm", [](const Cached::ChannelD& c) { return c.FlightLength(); }},
    {"DCA_btw_Daughters", "cm", [](const Cached::ChannelD& c) { return c.DCA_btw_Daughters(); }},
    {"CosTheta_btw_Daughters", "", [](const Cached::ChannelD& c) { return c.CosTheta_btw_Daughters(); }},
    {"Theta_btw_Daughters", "rad", [](const Cached::ChannelD& c) { return c.Theta_btw_Daughters(); }},
    {"Chi2NDF", "", [](const Cached::ChannelD& c) { return static_cast<double>(c.Chi2NDF); }},
    // -- topology, w.r.t. the primary vertex
    {"Qt_wrt_PV", "GeV/c", [](const Cached::ChannelD& c) { return c.Qt_wrt_PV(); }},
    {"Ql_wrt_PV", "GeV/c", [](const Cached::ChannelD& c) { return c.Ql_wrt_PV(); }},
    {"DCA_wrt_PV", "cm", [](const Cached::ChannelD& c) { return c.DCA_wrt_PV(); }},
    {"CPA_wrt_PV", "", [](const Cached::ChannelD& c) { return c.CPA_wrt_PV(); }},
    // -- the antilambda leg
    {"V0_DCA_wrt_SV", "cm", [](const Cached::ChannelD& c) { return c.V0_DCA_wrt_SV(); }},
    {"V0_CPA_wrt_SV", "", [](const Cached::ChannelD& c) { return c.V0_CPA_wrt_SV(); }},
    {"V0_DecayLength", "cm", [](const Cached::ChannelD& c) { return c.V0_DecayLength(); }},
    {"V0_Neg_Pt", "GeV/c", [](const Cached::ChannelD& c) { return c.V0_Neg_Pt(); }},
    {"V0_Neg_Eta", "", [](const Cached::ChannelD& c) { return c.V0_Neg_Eta(); }},
    {"V0_Pos_Pt", "GeV/c", [](const Cached::ChannelD& c) { return c.V0_Pos_Pt(); }},
    {"V0_Pos_Eta", "", [](const Cached::ChannelD& c) { return c.V0_Pos_Eta(); }},
    // -- the charged kaon leg
    {"Kaon_DCA_wrt_SV", "cm", [](const Cached::ChannelD& c) { return c.Kaon_DCA_wrt_SV(); }},
});

// == Channel H : antisexaquark + proton -> K+ + K+ + X == //

inline constexpr auto DB_ChannelH = std::to_array<Definition<Cached::ChannelH>>({
    // -- kinematics
    {"Pt", "GeV/c", [](const Cached::ChannelH& c) { return c.Pt(); }},
    {"P", "GeV/c", [](const Cached::ChannelH& c) { return c.P(); }},
    {"Mass", "GeV/c2", [](const Cached::ChannelH& c) { return c.Mass(); }},
    {"Eta", "", [](const Cached::ChannelH& c) { return c.Eta(); }},
    {"Rapidity", "", [](const Cached::ChannelH& c) { return c.Rapidity(); }},
    {"Mass_MinusNucleon", "GeV/c2", [](const Cached::ChannelH& c) { return c.Mass_MinusNucleon(); }},
    {"Rapidity_MinusNucleon", "", [](const Cached::ChannelH& c) { return c.Rapidity_MinusNucleon(); }},
    // -- topology, w.r.t. the secondary vertex
    {"SV_Radius2D", "cm", [](const Cached::ChannelH& c) { return c.SV_Radius2D(); }},
    {"SV_Radius3D", "cm", [](const Cached::ChannelH& c) { return c.SV_Radius3D(); }},
    {"FlightLength", "cm", [](const Cached::ChannelH& c) { return c.FlightLength(); }},
    {"DCA_btw_Daughters", "cm", [](const Cached::ChannelH& c) { return c.DCA_btw_Daughters(); }},
    {"CosTheta_btw_Daughters", "", [](const Cached::ChannelH& c) { return c.CosTheta_btw_Daughters(); }},
    {"Theta_btw_Daughters", "rad", [](const Cached::ChannelH& c) { return c.Theta_btw_Daughters(); }},
    {"Chi2NDF", "", [](const Cached::ChannelH& c) { return static_cast<double>(c.Chi2NDF); }},
    // -- topology, w.r.t. the primary vertex
    {"Qt_wrt_PV", "GeV/c", [](const Cached::ChannelH& c) { return c.Qt_wrt_PV(); }},
    {"Ql_wrt_PV", "GeV/c", [](const Cached::ChannelH& c) { return c.Ql_wrt_PV(); }},
    {"DCA_wrt_PV", "cm", [](const Cached::ChannelH& c) { return c.DCA_wrt_PV(); }},
    {"CPA_wrt_PV", "", [](const Cached::ChannelH& c) { return c.CPA_wrt_PV(); }},
    // -- the two charged kaon legs
    {"Kaon1_DCA_wrt_SV", "cm", [](const Cached::ChannelH& c) { return c.Kaon1_DCA_wrt_SV(); }},
    {"Kaon2_DCA_wrt_SV", "cm", [](const Cached::ChannelH& c) { return c.Kaon2_DCA_wrt_SV(); }},
});

// == H-dibaryon : H -> Lambda + Lambda == //

// NOTE: every `CV_*` getter returns `Common::DummyDouble` (-999) when the production-vertex constraint
//       did not converge -- see `Cached::Hdibaryon::HasCV()`. Cut on them only alongside `Chi2CV`, or the
//       dummy population lands inside a lower-bound cut and inflates whatever it inflates.

inline constexpr auto DB_Hdibaryon = std::to_array<Definition<Cached::Hdibaryon>>({
    // -- kinematics
    {"Pt", "GeV/c", [](const Cached::Hdibaryon& c) { return c.Pt(); }},
    {"P", "GeV/c", [](const Cached::Hdibaryon& c) { return c.P(); }},
    {"Mass", "GeV/c2", [](const Cached::Hdibaryon& c) { return c.Mass(); }},
    {"Eta", "", [](const Cached::Hdibaryon& c) { return c.Eta(); }},
    {"Rapidity", "", [](const Cached::Hdibaryon& c) { return c.Rapidity(); }},
    // -- topology, w.r.t. the decay vertex
    {"Decay_Radius2D", "cm", [](const Cached::Hdibaryon& c) { return c.Decay_Radius2D(); }},
    {"Decay_Radius3D", "cm", [](const Cached::Hdibaryon& c) { return c.Decay_Radius3D(); }},
    {"DecayLength", "cm", [](const Cached::Hdibaryon& c) { return c.DecayLength(); }},
    {"DCA_btw_Lambdas", "cm", [](const Cached::Hdibaryon& c) { return c.DCA_btw_Lambdas(); }},
    {"Dist_btw_LambdaDVs", "cm", [](const Cached::Hdibaryon& c) { return c.Dist_btw_LambdaDVs(); }},
    {"Chi2NDF", "", [](const Cached::Hdibaryon& c) { return static_cast<double>(c.Chi2NDF); }},
    // -- topology, w.r.t. the primary vertex
    {"DCAxy_wrt_PV", "cm", [](const Cached::Hdibaryon& c) { return c.DCAxy_wrt_PV(); }},
    {"DCAz_wrt_PV", "cm", [](const Cached::Hdibaryon& c) { return c.DCAz_wrt_PV(); }},
    {"DCA_wrt_PV", "cm", [](const Cached::Hdibaryon& c) { return c.DCA_wrt_PV(); }},
    {"CPA_wrt_PV", "", [](const Cached::Hdibaryon& c) { return c.CPA_wrt_PV(); }},
    // -- production-vertex-constrained fit
    {"Chi2CV", "", [](const Cached::Hdibaryon& c) { return static_cast<double>(c.Chi2CV); }},
    {"CV_Mass", "GeV/c2", [](const Cached::Hdibaryon& c) { return c.CV_Mass(); }},
    {"CV_Pt", "GeV/c", [](const Cached::Hdibaryon& c) { return c.CV_Pt(); }},
    {"CV_Radius2D", "cm", [](const Cached::Hdibaryon& c) { return c.CV_Radius2D(); }},
    {"CV_DecayLength_DivByErr", "", [](const Cached::Hdibaryon& c) { return c.CV_DecayLength_DivByErr(); }},
    // -- the first lambda leg
    {"Lambda1_DecayLength", "cm", [](const Cached::Hdibaryon& c) { return c.Lambda1_DecayLength(); }},
    {"Lambda1_CPA_wrt_DV", "", [](const Cached::Hdibaryon& c) { return c.Lambda1_CPA_wrt_DV(); }},
    {"Lambda1_DCA_wrt_DV", "cm", [](const Cached::Hdibaryon& c) { return c.Lambda1_DCA_wrt_DV(); }},
    {"Proton1_Pt", "GeV/c", [](const Cached::Hdibaryon& c) { return c.Proton1_Pt(); }},
    {"Proton1_Eta", "", [](const Cached::Hdibaryon& c) { return c.Proton1_Eta(); }},
    {"Pion1_Pt", "GeV/c", [](const Cached::Hdibaryon& c) { return c.Pion1_Pt(); }},
    {"Pion1_Eta", "", [](const Cached::Hdibaryon& c) { return c.Pion1_Eta(); }},
    // -- the second lambda leg
    {"Lambda2_DecayLength", "cm", [](const Cached::Hdibaryon& c) { return c.Lambda2_DecayLength(); }},
    {"Lambda2_CPA_wrt_DV", "", [](const Cached::Hdibaryon& c) { return c.Lambda2_CPA_wrt_DV(); }},
    {"Lambda2_DCA_wrt_DV", "cm", [](const Cached::Hdibaryon& c) { return c.Lambda2_DCA_wrt_DV(); }},
    {"Proton2_Pt", "GeV/c", [](const Cached::Hdibaryon& c) { return c.Proton2_Pt(); }},
    {"Proton2_Eta", "", [](const Cached::Hdibaryon& c) { return c.Proton2_Eta(); }},
    {"Pion2_Pt", "GeV/c", [](const Cached::Hdibaryon& c) { return c.Pion2_Pt(); }},
    {"Pion2_Eta", "", [](const Cached::Hdibaryon& c) { return c.Pion2_Eta(); }},
    // -- lambda-lambda correlations
    {"DeltaDecayLength_LL", "cm", [](const Cached::Hdibaryon& c) { return c.DeltaDecayLength_LL(); }},
    {"CosTheta_LL", "", [](const Cached::Hdibaryon& c) { return c.CosTheta_LL(); }},
    {"Theta_LL", "rad", [](const Cached::Hdibaryon& c) { return c.Theta_LL(); }},
    {"DeltaEta_LL", "", [](const Cached::Hdibaryon& c) { return c.DeltaEta_LL(); }},
    {"DeltaPhi_LL", "rad", [](const Cached::Hdibaryon& c) { return c.DeltaPhi_LL(); }},
    {"DeltaR_LL", "", [](const Cached::Hdibaryon& c) { return c.DeltaR_LL(); }},
    {"Asym_LL", "", [](const Cached::Hdibaryon& c) { return c.Asym_LL(); }},
    // -- proton-proton correlations
    {"CosTheta_pp", "", [](const Cached::Hdibaryon& c) { return c.CosTheta_pp(); }},
    {"Qrel", "GeV/c", [](const Cached::Hdibaryon& c) { return c.Qrel(); }},
    {"CosThetaStar_L1", "", [](const Cached::Hdibaryon& c) { return c.CosThetaStar_L1(); }},
    {"CosThetaStar_L2", "", [](const Cached::Hdibaryon& c) { return c.CosThetaStar_L2(); }},
    {"CosThetaStar_Pr1", "", [](const Cached::Hdibaryon& c) { return c.CosThetaStar_Pr1(); }},
    {"CosThetaStar_Pr2", "", [](const Cached::Hdibaryon& c) { return c.CosThetaStar_Pr2(); }},
});

}  // namespace Skimmer::Variables
