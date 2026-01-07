#pragma once

#include <Math/Point3D.h>
#include <Math/Vector3D.h>

#include <KFParticle.hxx>
#include <KFParticle_Math.hxx>

#include "Math/BaseMath.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::Fit {

struct alignas(T2S_SIMD_ALIGN) Sexaquark : KF::Particle {
    // constructors //
    Sexaquark() = delete;
    explicit Sexaquark(double nucleon_mass) : Nucleon_Mass{nucleon_mass} {};

    // utilities //
    [[nodiscard]] ROOT::Math::XYZPoint GetXYZ_AsROOT() const { return {X(), Y(), Z()}; }
    [[nodiscard]] ROOT::Math::XYZVector GetPxPyPz_AsROOT() const { return {Px(), Py(), Pz()}; }
    [[nodiscard]] double E_MinusNucleon() const { return E() - Nucleon_Mass; };
    [[nodiscard]] double Mass_MinusNucleon() const {
        double mass2{E_MinusNucleon() * E_MinusNucleon() - P2()};
        if (mass2 > 0.) return std::sqrt(mass2);
        return Const::DummyDouble;
    };
    [[nodiscard]] double AbsRapidity_MinusNucleon() const { return std::abs(std::log((E_MinusNucleon() + Pz()) / (E_MinusNucleon() - Pz())) / 2.); };

    // cuts //
    [[nodiscard]] double CPA_wrt(double ref_x, double ref_y, double ref_z) const {
        return Tree2Secondaries::Math::CosinePointingAngle(GetPxPyPz_AsROOT(), GetXYZ_AsROOT(), {ref_x, ref_y, ref_z});
    }
    [[nodiscard]] double DCA_wrt(double ref_x, double ref_y, double ref_z) const {
        return Tree2Secondaries::Math::FastDCALineVertex(GetPxPyPz_AsROOT(), GetXYZ_AsROOT(), {ref_x, ref_y, ref_z});
    }

    // member vars //
    double Nucleon_Mass{};
};

}  // namespace Tree2Secondaries::Fit
