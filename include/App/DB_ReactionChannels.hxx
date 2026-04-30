#pragma once

#include <array>
#include <cstddef>
#include <optional>
#include <span>

namespace Tree2Secondaries::ReactionChannels {

struct Definition {
    char name;                          // 'A'..'H'
    int nucleon_pdg;                    // struck nucleon
    std::span<const int> products_pdg;  // first gen. products
};

inline constexpr std::array kFirstGen_A = std::to_array<int>({-3122, 310});
inline constexpr std::array kFirstGen_B = std::to_array<int>({-3122, 310, -211, 211});
inline constexpr std::array kFirstGen_C = std::to_array<int>({-2212, 310, 310, 211});
inline constexpr std::array kFirstGen_D = std::to_array<int>({-3122, 321});
inline constexpr std::array kFirstGen_E = std::to_array<int>({-3122, 321, -211, 211});
inline constexpr std::array kFirstGen_F = std::to_array<int>({-2212, 321, 310, 211});
inline constexpr std::array kFirstGen_G = std::to_array<int>({-3312, -211});
inline constexpr std::array kFirstGen_H = std::to_array<int>({-2212, 321, 321, 111});

inline constexpr std::array DB = std::to_array<Definition>({
    {'A', 2112, kFirstGen_A},
    {'B', 2112, kFirstGen_B},
    {'C', 2112, kFirstGen_C},
    {'D', 2212, kFirstGen_D},
    {'E', 2212, kFirstGen_E},
    {'F', 2212, kFirstGen_F},
    {'G', 2112, kFirstGen_G},
    {'H', 2212, kFirstGen_H},
});

consteval std::size_t IndexOf(char c) {
    for (std::size_t rr = 0; rr < DB.size(); ++rr)
        if (DB[rr].name == c) return rr;
    throw "unknown reaction channel";
}

constexpr std::optional<Definition> Find(char c) {
    for (const auto& r : DB)
        if (r.name == c) return r;
    return std::nullopt;
}

constexpr std::optional<std::size_t> FindIndex(char c) {
    for (std::size_t rr = 0; rr < DB.size(); ++rr)
        if (DB[rr].name == c) return rr;
    return std::nullopt;
}

}  // namespace Tree2Secondaries::ReactionChannels
