#pragma once

#include <array>

#include <Math/Point3D.h>
#include <Math/Vector3D.h>

namespace R2DS::Seeder {

// # Structs # //

struct PCA {
    [[nodiscard]] double X() const noexcept { return xyz[0]; }
    [[nodiscard]] double Y() const noexcept { return xyz[1]; }
    [[nodiscard]] double Z() const noexcept { return xyz[2]; }
    [[nodiscard]] double Px() const noexcept { return mom[0]; }
    [[nodiscard]] double Py() const noexcept { return mom[1]; }
    [[nodiscard]] double Pz() const noexcept { return mom[2]; }

    [[nodiscard]] ROOT::Math::XYZPoint GetXYZ_AsROOT() const { return {X(), Y(), Z()}; }
    [[nodiscard]] ROOT::Math::XYZVector GetPxPyPz_AsROOT() const { return {Px(), Py(), Pz()}; }

    std::array<double, 3> xyz{};
    std::array<double, 3> mom{};
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

struct Deriv {
    std::array<double, 6> ds_dr{};
    std::array<double, 6> ds_dr1{};
};

struct Result {
    Seed seed;
    Deriv deriv;
};

}  // namespace R2DS::Seeder
