#pragma once

#include <array>
#include <optional>
#include <span>
#include <string_view>

namespace Tree2Secondaries::Particles {

inline constexpr std::array DecayProducts_PiZero = std::to_array<int>({22, 22});
inline constexpr std::array DecayProducts_KaonZeroShort = std::to_array<int>({-211, 211});
inline constexpr std::array DecayProducts_AntiLambda = std::to_array<int>({-2212, 211});
inline constexpr std::array DecayProducts_Lambda = std::to_array<int>({2212, -211});
inline constexpr std::array DecayProducts_XiPlus = std::to_array<int>({-3122, 211});  // Xi+ -> AntiLambda Pi+

struct Definition {
    std::string_view name;               // canonical name
    std::string_view acronym;            // short name
    int pdg_code;                        // PDG MC numbering scheme
    int charge;                          // in units of e
    double mass;                         // (GeV/c^2)
    double ctau;                         // (cm)
    std::span<const int> daughters_pdg;  // NOTE: only detectable charged modes
};

inline constexpr std::array DB = std::to_array<Definition>({
    {"Photon", "G", 22, 0, 0.0, 1e20, {}},
    //
    {"Electron", "EL", 11, -1, 0.00051100, 1e20, {}},
    {"Positron", "PO", -11, +1, 0.00051100, 1e20, {}},
    //
    {"PiMinus", "PM", -211, -1, 0.13957040, 780.45, {}},
    {"PiPlus", "PP", 211, +1, 0.13957040, 780.45, {}},
    {"PiZero", "P0", 111, 0, 0.13497680, 2.5e-6, DecayProducts_PiZero},
    //
    {"NegKaon", "NK", -321, -1, 0.49367700, 371.1, {}},
    {"PosKaon", "PK", 321, +1, 0.49367700, 371.1, {}},
    {"KaonZeroShort", "K0S", 310, 0, 0.49761100, 2.6844, DecayProducts_KaonZeroShort},
    //
    {"AntiProton", "AP", -2212, -1, 0.93827210, 1e20, {}},
    {"Proton", "P", 2212, +1, 0.93827210, 1e20, {}},
    //
    {"AntiNeutron", "AN", -2112, 0, 0.93956540, 1e20, {}},
    {"Neutron", "N", 2112, 0, 0.93956540, 1e20, {}},
    //
    {"AntiLambda", "AL", -3122, 0, 1.1156830, 7.89, DecayProducts_AntiLambda},
    {"Lambda", "L", 3122, 0, 1.1156830, 7.89, DecayProducts_Lambda},
    //
    {"XiPlus", "XP", -3312, +1, 1.3217100, 4.91, DecayProducts_XiPlus},
});

consteval Definition Particle(std::string_view name) {
    for (const auto& i : DB)
        if (i.name == name) return i;
    throw "unknown particle name";
}

constexpr std::optional<Definition> FindParticle(int pdg_code) {
    for (const auto& i : DB)
        if (i.pdg_code == pdg_code) return i;
    return std::nullopt;
}

}  // namespace Tree2Secondaries::Particles
