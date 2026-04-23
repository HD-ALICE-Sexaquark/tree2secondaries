#pragma once

#include <cmath>
#include <cstdlib>
#include <optional>

#include <Math/Point3D.h>
#include <Math/Vector3D.h>

#include "View/Reconstructed/ViewRecTrack.hxx"
#include "View/Reconstructed/ViewRecV0.hxx"

namespace Tree2Secondaries::Legacy {

// Helper //

template <typename Scalar>
static Scalar IJ(Scalar i, Scalar j) {
    return (j <= i) ? i * (i + 1) / 2 + j : j * (j + 1) / 2 + i;
}

struct Particle {

    void Initialize();

    // Constructors //

    Particle() { Initialize(); };
    static Particle FromTrack(const View::Rec::Track &track, double mass);
    static Particle FromV0(const View::Rec::V0 &v0);

    // Named Accesors //

    [[nodiscard]] double X() const noexcept { return fP[0]; }
    [[nodiscard]] double Y() const noexcept { return fP[1]; }
    [[nodiscard]] double Z() const noexcept { return fP[2]; }
    [[nodiscard]] double Px() const noexcept { return fP[3]; }
    [[nodiscard]] double Py() const noexcept { return fP[4]; }
    [[nodiscard]] double Pz() const noexcept { return fP[5]; }
    [[nodiscard]] double E() const noexcept { return fP[6]; }
    [[nodiscard]] double S() const noexcept { return fP[7]; }

    [[nodiscard]] double SigmaX2() const noexcept { return fC[IJ(0, 0)]; }
    [[nodiscard]] double SigmaXY() const noexcept { return fC[IJ(1, 0)]; }
    [[nodiscard]] double SigmaY2() const noexcept { return fC[IJ(1, 1)]; }
    [[nodiscard]] double SigmaXZ() const noexcept { return fC[IJ(2, 0)]; }
    [[nodiscard]] double SigmaYZ() const noexcept { return fC[IJ(2, 1)]; }
    [[nodiscard]] double SigmaZ2() const noexcept { return fC[IJ(2, 2)]; }
    [[nodiscard]] double SigmaXPx() const noexcept { return fC[IJ(3, 0)]; }
    [[nodiscard]] double SigmaYPx() const noexcept { return fC[IJ(3, 1)]; }
    [[nodiscard]] double SigmaZPx() const noexcept { return fC[IJ(3, 2)]; }
    [[nodiscard]] double SigmaPx2() const noexcept { return fC[IJ(3, 3)]; }
    [[nodiscard]] double SigmaXPy() const noexcept { return fC[IJ(4, 0)]; }
    [[nodiscard]] double SigmaYPy() const noexcept { return fC[IJ(4, 1)]; }
    [[nodiscard]] double SigmaZPy() const noexcept { return fC[IJ(4, 2)]; }
    [[nodiscard]] double SigmaPxPy() const noexcept { return fC[IJ(4, 3)]; }
    [[nodiscard]] double SigmaPy2() const noexcept { return fC[IJ(4, 4)]; }
    [[nodiscard]] double SigmaXPz() const noexcept { return fC[IJ(5, 0)]; }
    [[nodiscard]] double SigmaYPz() const noexcept { return fC[IJ(5, 1)]; }
    [[nodiscard]] double SigmaZPz() const noexcept { return fC[IJ(5, 2)]; }
    [[nodiscard]] double SigmaPxPz() const noexcept { return fC[IJ(5, 3)]; }
    [[nodiscard]] double SigmaPyPz() const noexcept { return fC[IJ(5, 4)]; }
    [[nodiscard]] double SigmaPz2() const noexcept { return fC[IJ(5, 5)]; }

    [[nodiscard]] double SigmaXE() const noexcept { return fC[IJ(6, 0)]; }
    [[nodiscard]] double SigmaYE() const noexcept { return fC[IJ(6, 1)]; }
    [[nodiscard]] double SigmaZE() const noexcept { return fC[IJ(6, 2)]; }
    [[nodiscard]] double SigmaPxE() const noexcept { return fC[IJ(6, 3)]; }
    [[nodiscard]] double SigmaPyE() const noexcept { return fC[IJ(6, 4)]; }
    [[nodiscard]] double SigmaPzE() const noexcept { return fC[IJ(6, 5)]; }
    [[nodiscard]] double SigmaE2() const noexcept { return fC[IJ(6, 6)]; }
    [[nodiscard]] double SigmaXS() const noexcept { return fC[IJ(7, 0)]; }
    [[nodiscard]] double SigmaYS() const noexcept { return fC[IJ(7, 1)]; }
    [[nodiscard]] double SigmaZS() const noexcept { return fC[IJ(7, 2)]; }
    [[nodiscard]] double SigmaPxS() const noexcept { return fC[IJ(7, 3)]; }
    [[nodiscard]] double SigmaPyS() const noexcept { return fC[IJ(7, 4)]; }
    [[nodiscard]] double SigmaPzS() const noexcept { return fC[IJ(7, 5)]; }
    [[nodiscard]] double SigmaES() const noexcept { return fC[IJ(7, 6)]; }
    [[nodiscard]] double SigmaS2() const noexcept { return fC[IJ(7, 7)]; }

    [[nodiscard]] int Charge() const noexcept { return fQ; }

    // Derived Quantities //

    [[nodiscard]] std::array<double, 3> GetXYZ() const { return {X(), Y(), Z()}; }
    [[nodiscard]] std::array<double, 3> GetPxPyPz() const { return {Px(), Py(), Pz()}; }
    [[nodiscard]] ROOT::Math::XYZPoint GetXYZ_AsROOT() const { return {X(), Y(), Z()}; }
    [[nodiscard]] ROOT::Math::XYZVector GetPxPyPz_AsROOT() const { return {Px(), Py(), Pz()}; }

    [[nodiscard]] double SquaredRadius2D() const { return X() * X() + Y() * Y(); }
    [[nodiscard]] double SquaredRadius3D() const { return SquaredRadius2D() + Z() * Z(); }
    [[nodiscard]] double SquaredPt() const noexcept { return Px() * Px() + Py() * Py(); }
    [[nodiscard]] double SquaredMomentum() const noexcept { return SquaredPt() + Pz() * Pz(); }

    [[nodiscard]] double Chi2NDF() const noexcept { return fChi2 / static_cast<double>(fNDF); }

    // NOTE: not guarded on purpose
    [[nodiscard]] double Eta() const { return std::atanh(Pz() / Momentum()); }

    [[nodiscard]] double Pseudorapidity() const { return Eta(); }

    // NOTE: not guarded on purpose
    [[nodiscard]] double Rapidity() const { return std::log((E() + Pz()) / (E() - Pz())) / 2.; }

    [[nodiscard]] double Pt() const noexcept { return std::sqrt(SquaredPt()); }
    [[nodiscard]] double Momentum() const noexcept { return std::sqrt(SquaredMomentum()); }
    [[nodiscard]] double Radius2D() const { return std::sqrt(SquaredRadius2D()); }
    [[nodiscard]] double Radius3D() const { return std::sqrt(SquaredRadius3D()); }
    [[nodiscard]] std::optional<double> Mass() const {
        double m_squared = (E() - Momentum()) * (E() + Momentum());
        if (m_squared < 0.) return std::nullopt;  // protection
        return std::sqrt(m_squared);
    }

    [[nodiscard]] double AbsZ() const { return std::abs(Z()); }
    [[nodiscard]] double AbsEta() const { return std::abs(Eta()); }

    // Member Variables //

    double fC[36];
    double fP[8];
    double fChi2;
    int fNDF;
    int fQ;
};

}  // namespace Tree2Secondaries::Legacy
