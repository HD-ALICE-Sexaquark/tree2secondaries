#ifndef T2S_SECONDARY_SEXAQUARK_HXX
#define T2S_SECONDARY_SEXAQUARK_HXX

#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

#include "Math/Common.hxx"
#include "Secondary/Reconstructed.hxx"

namespace Tree2Secondaries {

class Sexaquark : public Reconstructed {
   public:
    Sexaquark(size_t idx, const ROOT::Math::XYZPoint& vtx, const ROOT::Math::PxPyPzEVector& momentum, double nucleon_mass)
        : Reconstructed{0,                                     //
                        vtx,                                   //
                        ROOT::Math::PxPyPzMVector{momentum}},  //
          fIndex{idx},
          fNoFermiEnergy{momentum.E() + nucleon_mass} {}
    ~Sexaquark() override = default;

    [[nodiscard]] size_t Index() const override { return fIndex; }

    [[nodiscard]] ROOT::Math::XYZPoint SecondaryVertex() const { return fMeasuredVtx; }
    [[nodiscard]] double SecondaryX() const { return fMeasuredVtx.X(); }
    [[nodiscard]] double SecondaryY() const { return fMeasuredVtx.Y(); }
    [[nodiscard]] double SecondaryZ() const { return fMeasuredVtx.Z(); }
    [[nodiscard]] double Radius() const { return fMeasuredVtx.Rho(); }

    // Line Parametrization: //
    //   X(s)  = X0 + Px s   //
    //   Y(s)  = Y0 + Py s   //
    //   Z(s)  = Z0 + Pz s   //
    [[nodiscard]] double X(double s) const override { return SecondaryX() + Px(s) * s; }
    [[nodiscard]] double Y(double s) const override { return SecondaryY() + Py(s) * s; }
    [[nodiscard]] double Z(double s) const override { return SecondaryZ() + Pz(s) * s; }
    [[nodiscard]] double Px(double s [[maybe_unused]] = 0.) const override { return fMeasuredMom.Px(); }  // PENDING: wtf
    [[nodiscard]] double Py(double s [[maybe_unused]] = 0.) const override { return fMeasuredMom.Py(); }  // PENDING: wtf
    [[nodiscard]] double Pz(double s [[maybe_unused]] = 0.) const override { return fMeasuredMom.Pz(); }  // PENDING: wtf

    [[nodiscard]] double Rapidity() const { return fMeasuredMom.Rapidity(); }
    [[nodiscard]] double Mass() const { return fMeasuredMom.M(); }
    [[nodiscard]] double MassAsDecay() const { return std::sqrt(fNoFermiEnergy * fNoFermiEnergy - fMeasuredMom.P2()); }
    [[nodiscard]] double EnergyAsDecay() const { return fNoFermiEnergy; }

    [[nodiscard]] double CPAwrt(const ROOT::Math::XYZPoint& v) const { return Math::CosinePointingAngle(PxPyPz(), SecondaryVertex(), v); }

   protected:
    size_t fIndex{0};
    double fNoFermiEnergy;
};

}  // namespace Tree2Secondaries

#endif  // T2S_SECONDARY_SEXAQUARK_HXX
