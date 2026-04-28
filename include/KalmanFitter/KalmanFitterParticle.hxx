#pragma once

#include <cmath>
#include <format>
#include <optional>

#include <Math/Point3D.h>
#include <Math/Vector3D.h>

#include <Eigen/Eigen>

#include "Math/Constants.hxx"
#include "View/Reconstructed/ViewRecTrack.hxx"
#include "View/Reconstructed/ViewRecV0.hxx"
#if T2S_LEGACY_KF
#include "Legacy/LegacyParticle.hxx"
#endif

namespace Tree2Secondaries::KalmanFitter {

// # KF::Particle # //

struct alignas(T2S_SIMD_ALIGN) Particle {

    // Constructors //

    Particle() = default;
    static Particle FromTrack(const View::Rec::Track &track, double mass);
    static Particle FromV0(const View::Rec::V0 &v0);
#if T2S_LEGACY_KF
    static Particle FromLegacy(const Legacy::Particle &part);
#endif

    // Named Accesors //

    [[nodiscard]] double X() const noexcept { return fP(0); }
    [[nodiscard]] double Y() const noexcept { return fP(1); }
    [[nodiscard]] double Z() const noexcept { return fP(2); }
    [[nodiscard]] double Px() const noexcept { return fP(3); }
    [[nodiscard]] double Py() const noexcept { return fP(4); }
    [[nodiscard]] double Pz() const noexcept { return fP(5); }
    [[nodiscard]] double E() const noexcept { return fP(6); }
    [[nodiscard]] double S() const noexcept { return fP(7); }

    [[nodiscard]] double CovX2() const noexcept { return fC(0, 0); }
    [[nodiscard]] double CovXY() const noexcept { return fC(1, 0); }
    [[nodiscard]] double CovY2() const noexcept { return fC(1, 1); }
    [[nodiscard]] double CovXZ() const noexcept { return fC(2, 0); }
    [[nodiscard]] double CovYZ() const noexcept { return fC(2, 1); }
    [[nodiscard]] double CovZ2() const noexcept { return fC(2, 2); }
    [[nodiscard]] double CovXPx() const noexcept { return fC(3, 0); }
    [[nodiscard]] double CovYPx() const noexcept { return fC(3, 1); }
    [[nodiscard]] double CovZPx() const noexcept { return fC(3, 2); }
    [[nodiscard]] double CovPx2() const noexcept { return fC(3, 3); }
    [[nodiscard]] double CovXPy() const noexcept { return fC(4, 0); }
    [[nodiscard]] double CovYPy() const noexcept { return fC(4, 1); }
    [[nodiscard]] double CovZPy() const noexcept { return fC(4, 2); }
    [[nodiscard]] double CovPxPy() const noexcept { return fC(4, 3); }
    [[nodiscard]] double CovPy2() const noexcept { return fC(4, 4); }
    [[nodiscard]] double CovXPz() const noexcept { return fC(5, 0); }
    [[nodiscard]] double CovYPz() const noexcept { return fC(5, 1); }
    [[nodiscard]] double CovZPz() const noexcept { return fC(5, 2); }
    [[nodiscard]] double CovPxPz() const noexcept { return fC(5, 3); }
    [[nodiscard]] double CovPyPz() const noexcept { return fC(5, 4); }
    [[nodiscard]] double CovPz2() const noexcept { return fC(5, 5); }

    [[nodiscard]] double CovXE() const noexcept { return fC(6, 0); }
    [[nodiscard]] double CovYE() const noexcept { return fC(6, 1); }
    [[nodiscard]] double CovZE() const noexcept { return fC(6, 2); }
    [[nodiscard]] double CovPxE() const noexcept { return fC(6, 3); }
    [[nodiscard]] double CovPyE() const noexcept { return fC(6, 4); }
    [[nodiscard]] double CovPzE() const noexcept { return fC(6, 5); }
    [[nodiscard]] double CovE2() const noexcept { return fC(6, 6); }
    [[nodiscard]] double CovXS() const noexcept { return fC(7, 0); }
    [[nodiscard]] double CovYS() const noexcept { return fC(7, 1); }
    [[nodiscard]] double CovZS() const noexcept { return fC(7, 2); }
    [[nodiscard]] double CovPxS() const noexcept { return fC(7, 3); }
    [[nodiscard]] double CovPyS() const noexcept { return fC(7, 4); }
    [[nodiscard]] double CovPzS() const noexcept { return fC(7, 5); }
    [[nodiscard]] double CovES() const noexcept { return fC(7, 6); }
    [[nodiscard]] double CovS2() const noexcept { return fC(7, 7); }

    [[nodiscard]] double VarX() const noexcept { return CovX2(); }
    [[nodiscard]] double VarY() const noexcept { return CovY2(); }
    [[nodiscard]] double VarZ() const noexcept { return CovZ2(); }
    [[nodiscard]] double VarPx() const noexcept { return CovPx2(); }
    [[nodiscard]] double VarPy() const noexcept { return CovPy2(); }
    [[nodiscard]] double VarPz() const noexcept { return CovPz2(); }
    [[nodiscard]] double VarE() const noexcept { return CovE2(); }
    [[nodiscard]] double VarS() const noexcept { return CovS2(); }

    template <unsigned int N>
    [[nodiscard]] std::array<float, (N * (N + 1)) / 2> Cov() const {
        std::array<float, (N * (N + 1)) / 2> out{};
        for (unsigned int i = 0, k = 0; i < N; ++i) {
            for (unsigned int j = 0; j <= i; ++j, ++k) {
                out[k] = static_cast<float>(fC(i, j));
            }
        }
        return out;
    }

    template <unsigned int N>
    void AppendCov(std::vector<float> &out) const {
        for (unsigned int i = 0; i < N; ++i) {
            for (unsigned int j = 0; j <= i; ++j) {
                out.push_back(static_cast<float>(fC(i, j)));
            }
        }
    }

    [[nodiscard]] double Chi2() const noexcept { return fChi2; }
    [[nodiscard]] int NDF() const noexcept { return fNDF; }
    [[nodiscard]] int Charge() const noexcept { return fQ; }

    // Derived Quantities //

    [[nodiscard]] std::array<double, 3> GetXYZ() const { return {X(), Y(), Z()}; }
    [[nodiscard]] std::array<double, 3> GetPxPyPz() const { return {Px(), Py(), Pz()}; }
    [[nodiscard]] ROOT::Math::XYZPoint GetXYZ_AsROOT() const { return {X(), Y(), Z()}; }
    [[nodiscard]] ROOT::Math::XYZVector GetPxPyPz_AsROOT() const { return {Px(), Py(), Pz()}; }

    [[nodiscard]] double SquaredRadius2D() const noexcept { return X() * X() + Y() * Y(); }
    [[nodiscard]] double SquaredRadius3D() const noexcept { return SquaredRadius2D() + Z() * Z(); }
    [[nodiscard]] double SquaredPt() const noexcept { return Px() * Px() + Py() * Py(); }
    [[nodiscard]] double SquaredMomentum() const noexcept { return SquaredPt() + Pz() * Pz(); }

    [[nodiscard]] double Chi2NDF() const noexcept { return fChi2 / static_cast<double>(fNDF); }

    // NOTE: unguarded on purpose
    [[nodiscard]] double Eta() const { return std::atanh(Pz() / Momentum()); }
    [[nodiscard]] double Pseudorapidity() const { return Eta(); }

    // NOTE: unguarded on purpose
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

    Eigen::Matrix<double, 8, 8> fC{Eigen::Matrix<double, 8, 8>::Zero()};  // symmetric, only lower triangle is used
    Eigen::Vector<double, 8> fP{Eigen::Vector<double, 8>::Zero()};
    double fChi2{};
    int fNDF{Const::DummyInt};
    int fQ{};
};

}  // namespace Tree2Secondaries::KalmanFitter

// # KF::Particle Print Formatter # //

template <>
struct std::formatter<Tree2Secondaries::KalmanFitter::Particle> {
    constexpr auto parse(std::format_parse_context &ctx) { return ctx.begin(); }
    auto format(const Tree2Secondaries::KalmanFitter::Particle &p, std::format_context &ctx) const {
        auto out = ctx.out();
        out = std::format_to(out, "\n");
        out = std::format_to(out, "(X,Y,Z,S)    = ({:13.6e}, {:13.6e}, {:13.6e}, {:13.6e})\n", p.X(), p.Y(), p.Z(), p.S());
        out = std::format_to(out, "(Px,Py,Pz,E) = ({:13.6e}, {:13.6e}, {:13.6e}, {:13.6e})\n", p.Px(), p.Py(), p.Pz(), p.E());
        out = std::format_to(out, "Mass         = {:13.6e}\n", p.Mass().value_or(Tree2Secondaries::Const::DummyInt));
        out = std::format_to(out, "Radius2D     = {:13.6e}\n", p.Radius2D());
        out = std::format_to(out, "Chi2/NDF     = {:13.6e} / {} = {:13.6e}\n", p.Chi2(), p.NDF(), p.Chi2NDF());
        // out = std::format_to(out, "fC           = {}\n", p.fC); // PENDING
        out = std::format_to(out, "\n");
        return out;
    }
};
