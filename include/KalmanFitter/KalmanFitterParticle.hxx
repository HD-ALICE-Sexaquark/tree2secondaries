#pragma once

#include <cmath>
#include <format>
#include <optional>

#include <Math/Point3D.h>
#include <Math/Vector3D.h>

#include <Eigen/Eigen>

#include "common/Constants.hpp"
#include "common/DB_Particles.hpp"
#include "common/POD_PreFoundLambda.hpp"
#include "common/POD_Track.hpp"
#include "common/POD_V0.hpp"

#include "App/Utilities.hxx"  // NOTE: don't remove! print formatter below needs it

namespace T2DS::KF {

// ## Constants ## //

static constexpr double Initial_Css = 1.;
static constexpr int Initial_NDF = -1;

// -- mass constraint's Newton solver
static constexpr int MassConstraint_MaxIter = 10;  // PENDING: maybe it's already done in 2-3 steps
static constexpr double MassConstraint_Tolerance = 1.E-10;
static constexpr double MassConstraint_MinDenom = 1.E-10;
static constexpr double MassConstraint_MinVariance = 1.E-20;

// ## KF::Particle ## //

struct Particle {

    // Constructors //

    // The mass always comes from the particle's identity, never from a parallel argument -- so the two can't disagree.
    // `on_shell` re-asserts an exact mass hypothesis on a composite, which is what a mother mass constraint applied
    // during a previous fit amounts to. It has to be opt-in and explicit because the hypothesis does not survive the
    // POD round-trip: `POD::V0` and `POD::Extended::PreFoundLambda` store (Px,Py,Pz,E) but no mass bookkeeping.

    Particle() = default;
    static Particle FromTrack(const POD::Track &v, const DB::Particles::Definition &pid);
    static Particle FromV0(const POD::V0 &v, const DB::Particles::Definition &pid, bool on_shell);
    static Particle FromPreFoundLambda(const POD::Extended::PreFoundLambda &l, const DB::Particles::Definition &pid, bool on_shell);

    // Modifier //

    void SetMassBookkeeping(const DB::Particles::Definition &pid, bool on_shell);

    // Named Accesors //

    [[nodiscard]] double X() const { return fP(0); }
    [[nodiscard]] double Y() const { return fP(1); }
    [[nodiscard]] double Z() const { return fP(2); }
    [[nodiscard]] double Px() const { return fP(3); }
    [[nodiscard]] double Py() const { return fP(4); }
    [[nodiscard]] double Pz() const { return fP(5); }
    [[nodiscard]] double E() const { return fP(6); }
    [[nodiscard]] double S() const { return fP(7); }

    [[nodiscard]] double CovX2() const { return fC(0, 0); }
    [[nodiscard]] double CovXY() const { return fC(1, 0); }
    [[nodiscard]] double CovY2() const { return fC(1, 1); }
    [[nodiscard]] double CovXZ() const { return fC(2, 0); }
    [[nodiscard]] double CovYZ() const { return fC(2, 1); }
    [[nodiscard]] double CovZ2() const { return fC(2, 2); }
    [[nodiscard]] double CovPxX() const { return fC(3, 0); }
    [[nodiscard]] double CovPxY() const { return fC(3, 1); }
    [[nodiscard]] double CovPxZ() const { return fC(3, 2); }
    [[nodiscard]] double CovPx2() const { return fC(3, 3); }
    [[nodiscard]] double CovPyX() const { return fC(4, 0); }
    [[nodiscard]] double CovPyY() const { return fC(4, 1); }
    [[nodiscard]] double CovPyZ() const { return fC(4, 2); }
    [[nodiscard]] double CovPyPx() const { return fC(4, 3); }
    [[nodiscard]] double CovPy2() const { return fC(4, 4); }
    [[nodiscard]] double CovPzX() const { return fC(5, 0); }
    [[nodiscard]] double CovPzY() const { return fC(5, 1); }
    [[nodiscard]] double CovPzZ() const { return fC(5, 2); }
    [[nodiscard]] double CovPzPx() const { return fC(5, 3); }
    [[nodiscard]] double CovPzPy() const { return fC(5, 4); }
    [[nodiscard]] double CovPz2() const { return fC(5, 5); }

    [[nodiscard]] double CovEX() const { return fC(6, 0); }
    [[nodiscard]] double CovEY() const { return fC(6, 1); }
    [[nodiscard]] double CovEZ() const { return fC(6, 2); }
    [[nodiscard]] double CovEPx() const { return fC(6, 3); }
    [[nodiscard]] double CovEPy() const { return fC(6, 4); }
    [[nodiscard]] double CovEPz() const { return fC(6, 5); }
    [[nodiscard]] double CovE2() const { return fC(6, 6); }

    [[nodiscard]] double CovSX() const { return fC(7, 0); }
    [[nodiscard]] double CovSY() const { return fC(7, 1); }
    [[nodiscard]] double CovSZ() const { return fC(7, 2); }
    [[nodiscard]] double CovSPx() const { return fC(7, 3); }
    [[nodiscard]] double CovSPy() const { return fC(7, 4); }
    [[nodiscard]] double CovSPz() const { return fC(7, 5); }
    [[nodiscard]] double CovSE() const { return fC(7, 6); }
    [[nodiscard]] double CovS2() const { return fC(7, 7); }

    [[nodiscard]] double VarX() const { return CovX2(); }
    [[nodiscard]] double VarY() const { return CovY2(); }
    [[nodiscard]] double VarZ() const { return CovZ2(); }
    [[nodiscard]] double VarPx() const { return CovPx2(); }
    [[nodiscard]] double VarPy() const { return CovPy2(); }
    [[nodiscard]] double VarPz() const { return CovPz2(); }
    [[nodiscard]] double VarE() const { return CovE2(); }
    [[nodiscard]] double VarS() const { return CovS2(); }

    template <unsigned int N>
    [[nodiscard]] std::array<float, N> State() const {
        std::array<float, N> out{};
        for (unsigned int i = 0; i < N; ++i) out[i] = static_cast<float>(fP(i));
        return out;
    }

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

    [[nodiscard]] double Chi2() const { return fChi2; }
    [[nodiscard]] int NDF() const { return fNDF; }
    [[nodiscard]] int Charge() const { return fQ; }

    // Derived Quantities //

    [[nodiscard]] std::array<double, 3> GetXYZ() const { return {X(), Y(), Z()}; }
    [[nodiscard]] std::array<double, 3> GetPxPyPz() const { return {Px(), Py(), Pz()}; }
    [[nodiscard]] ROOT::Math::XYZPoint GetXYZ_AsROOT() const { return {X(), Y(), Z()}; }
    [[nodiscard]] ROOT::Math::XYZVector GetPxPyPz_AsROOT() const { return {Px(), Py(), Pz()}; }

    [[nodiscard]] double SquaredRadius2D() const { return X() * X() + Y() * Y(); }
    [[nodiscard]] double SquaredRadius3D() const { return SquaredRadius2D() + Z() * Z(); }
    [[nodiscard]] double SquaredPt() const { return Px() * Px() + Py() * Py(); }
    [[nodiscard]] double SquaredMomentum() const { return SquaredPt() + Pz() * Pz(); }

    [[nodiscard]] double Chi2NDF() const { return fChi2 / static_cast<double>(fNDF); }

    // NOTE: unguarded on purpose
    [[nodiscard]] double Eta() const { return std::atanh(Pz() / Momentum()); }
    [[nodiscard]] double Pseudorapidity() const { return Eta(); }

    // NOTE: unguarded on purpose
    [[nodiscard]] double Rapidity() const { return std::log((E() + Pz()) / (E() - Pz())) / 2.; }

    [[nodiscard]] double Pt() const { return std::sqrt(SquaredPt()); }
    [[nodiscard]] double Momentum() const { return std::sqrt(SquaredMomentum()); }
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

    Eigen::Matrix<double, 8, 8> fC{Eigen::Matrix<double, 8, 8>::Zero()};  // full symmetric
    Eigen::Vector<double, 8> fP{Eigen::Vector<double, 8>::Zero()};
    std::optional<double> fMassHypo;  // exact mass hypothesis, if the particle has one
    double fChi2{};
    double fSumDaughterMass{0.};  // sum of the daughters' masses, i.e. the lowest physical mass
    int fNDF{Initial_NDF};
    int fQ{};
};

}  // namespace T2DS::KF

// # KF::Particle Print Formatter # //

template <>
struct std::formatter<T2DS::KF::Particle> {
    constexpr auto parse(std::format_parse_context &ctx) { return ctx.begin(); }
    auto format(const T2DS::KF::Particle &p, std::format_context &ctx) const {
        auto out = ctx.out();
        out = std::format_to(out, "\n");
        out = std::format_to(out, "(X,Y,Z,S)    = ({:13.6e}, {:13.6e}, {:13.6e}, {:13.6e})\n", p.X(), p.Y(), p.Z(), p.S());
        out = std::format_to(out, "(Px,Py,Pz,E) = ({:13.6e}, {:13.6e}, {:13.6e}, {:13.6e})\n", p.Px(), p.Py(), p.Pz(), p.E());
        out = std::format_to(out, "Mass         = {:13.6e}\n", p.Mass().value_or(Common::DummyInt));
        out = std::format_to(out, "Radius2D     = {:13.6e}\n", p.Radius2D());
        out = std::format_to(out, "Chi2/NDF     = {:13.6e} / {} = {:13.6e}\n", p.Chi2(), p.NDF(), p.Chi2NDF());
        out = std::format_to(out, "fC           = {}\n", p.fC);
        out = std::format_to(out, "\n");
        return out;
    }
};
