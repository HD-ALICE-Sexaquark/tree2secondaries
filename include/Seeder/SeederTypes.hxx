#pragma once

#include <array>
#include <utility>

#include "common/POD_PreFoundLambda.hpp"
#include "common/POD_Track.hpp"
#include "common/POD_V0.hpp"

namespace T2DS::Seeder {

struct PCA {
    std::array<double, 3> xyz{};
    std::array<double, 3> mom{};

    [[nodiscard]] double X() const noexcept { return xyz[0]; }
    [[nodiscard]] double Y() const noexcept { return xyz[1]; }
    [[nodiscard]] double Z() const noexcept { return xyz[2]; }
    [[nodiscard]] double Px() const noexcept { return mom[0]; }
    [[nodiscard]] double Py() const noexcept { return mom[1]; }
    [[nodiscard]] double Pz() const noexcept { return mom[2]; }
};

struct Seed {
    PCA pca;
    double theta{0.};
    double sin{0.};
    double cos{0.};
    double sB{0.};
    double cB{0.};
    double ds{0.};
};
using SeedsPair = std::pair<Seed, Seed>;

struct Deriv {
    std::array<double, 6> ds_dr{};
    std::array<double, 6> ds_dr1{};
};

struct Result {
    Seed seed;
    Deriv deriv;
};

struct State {
    double x{}, y{}, z{};
    double px{}, py{}, pz{};
    int charge{};

    static State FromTrack(const POD::Track& t) {
        return {static_cast<double>(t.X),
                static_cast<double>(t.Y),
                static_cast<double>(t.Z),
                static_cast<double>(t.Px),
                static_cast<double>(t.Py),
                static_cast<double>(t.Pz),
                t.Charge};
    }
    static State FromV0(const POD::V0& v) {
        return {static_cast<double>(v.Decay_X),
                static_cast<double>(v.Decay_Y),
                static_cast<double>(v.Decay_Z),
                static_cast<double>(v.Px),
                static_cast<double>(v.Py),
                static_cast<double>(v.Pz),
                0};
    }
    static State FromPreFoundLambda(const POD::Extended::PreFoundLambda& l) {
        return {static_cast<double>(l.Decay_X),
                static_cast<double>(l.Decay_Y),
                static_cast<double>(l.Decay_Z),
                static_cast<double>(l.Px),
                static_cast<double>(l.Py),
                static_cast<double>(l.Pz),
                0};
    }
};

}  // namespace T2DS::Seeder
