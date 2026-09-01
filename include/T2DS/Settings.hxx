#pragma once

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace T2DS {

enum EProgramMode : std::uint8_t { FINDER, VERIFIER };
inline constexpr std::array<const char *, 2> Name_ProgramMode{"FINDER", "VERIFIER"};

struct Settings {
    void Print() const;

    std::string PathOutputFile;
    std::vector<std::string> PathInputFiles;
    std::optional<unsigned long> LimitToNEvents;
    double SexaquarkMass{};
    EProgramMode Mode{EProgramMode::FINDER};
    bool IsMC{false};
};

}  // namespace T2DS
