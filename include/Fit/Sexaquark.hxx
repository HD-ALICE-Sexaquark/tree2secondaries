#pragma once

#include <KFParticle.hxx>
#include <KFParticle_Math.hxx>

#include "Math/Constants.hxx"
#include "Math/Math.hxx"

namespace Tree2Secondaries::Fit {

struct alignas(T2S_SIMD_ALIGN) Sexaquark : KF::Particle {
    // constructors //
    Sexaquark() = delete;
    explicit Sexaquark(double nucleon_mass) : Nucleon_Mass{nucleon_mass} {};

    // utilities //
    [[nodiscard]] double E_MinusNucleon() const { return E() - Nucleon_Mass; };
    [[nodiscard]] double Mass_MinusNucleon() const {
        double mass2{E_MinusNucleon() * E_MinusNucleon() - P2()};
        if (mass2 > 0.) return std::sqrt(mass2);
        return Const::DummyDouble;
    };

    // cuts //
    [[nodiscard]] double AbsRapidity_MinusNucleon() const { return std::abs(std::log((E_MinusNucleon() + Pz()) / (E_MinusNucleon() - Pz())) / 2.); };
    [[nodiscard]] double CPA_Point(float x, float y, float z) const {
        return Math::CosinePointingAngle({Px(), Py(), Pz()}, {X(), Y(), Z()}, {x, y, z});
    }

    // member vars //
    double Nucleon_Mass{};
};

}  // namespace Tree2Secondaries::Fit
