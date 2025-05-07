#ifndef T2S_SECONDARY_CHARGED_HXX
#define T2S_SECONDARY_CHARGED_HXX

#include "Math/Point3D.h"
#include "Math/Vector4D.h"

#include "Math/Constants.hxx"
#include "Math/VtxrResults.hxx"
#include "Secondary/Reconstructed.hxx"

namespace Tree2Secondaries {

//
class Charged final : public Reconstructed {
   public:
    Charged(const Charged &) = delete;
    Charged(Charged &&) noexcept = default;
    Charged &operator=(const Charged &) = delete;
    Charged &operator=(Charged &&) = default;
    ~Charged() final = default;

    Charged(size_t id, int charge, double bz, ROOT::Math::XYZPoint xyz0, ROOT::Math::PxPyPzMVector pxyz0)
        : Reconstructed{charge, std::move(xyz0), std::move(pxyz0)}, fTrackIndex{id}, fBz{bz} {}

    [[nodiscard]] size_t Index() const override { return fTrackIndex; }

    // convenient aliases //
    [[nodiscard]] double X0() const { return fMeasuredVtx.X(); }
    [[nodiscard]] double Y0() const { return fMeasuredVtx.Y(); }
    [[nodiscard]] double Z0() const { return fMeasuredVtx.Z(); }
    [[nodiscard]] double Px0() const { return fMeasuredMom.Px(); }
    [[nodiscard]] double Py0() const { return fMeasuredMom.Py(); }
    [[nodiscard]] double Pz0() const { return fMeasuredMom.Pz(); }
    [[nodiscard]] double Pt0() const { return fMeasuredMom.Pt(); }

    // Parametrization of charged particles: (helix)                    //
    //   X(s)  = X0 + Sin(omega s) Px/omega + (1-Cos(omega s)) Py/omega //
    //   Y(s)  = Y0 - (1-Cos(omega s)) Px/omega + Sin(omega s) Py/omega //
    //   Z(s)  = Z0 + Pz s                                              //
    //   Px(s) = Cos(omega s) Px + Sin(omega s) Py                      //
    //   Py(s) = -Sin(omega s) Px + Cos(omega s) Py                     //
    //   Pz    = (constant)                                             //
    [[nodiscard]] double Omega() const { return fBz * Charge() * Const::Kappa; }
    [[nodiscard]] double X(double s) const override {
        return X0() + std::sin(Omega() * s) * Px0() / Omega() + (1 - std::cos(Omega() * s)) * Py0() / Omega();
    }
    [[nodiscard]] double Y(double s) const override {
        return Y0() - (1 - std::cos(Omega() * s)) * Px0() / Omega() + std::sin(Omega() * s) * Py0() / Omega();
    }
    [[nodiscard]] double Z(double s) const override { return Z0() + Pz0() * s; }

    [[nodiscard]] double Px(double s) const override { return Px0() * std::cos(Omega() * s) + Py0() * std::sin(Omega() * s); }
    [[nodiscard]] double Py(double s) const override { return -Px0() * std::sin(Omega() * s) + Py0() * std::cos(Omega() * s); }
    [[nodiscard]] double Pz(double s [[maybe_unused]]) const override { return Pz0(); }  // PENDING: wtf

    /*
        void SetInitialState(double x0, double y0, double z0, double px0, double py0, double pz0, double mass) {
            fMeasuredVtx.SetCoordinates(x0, y0, z0);
            fMeasuredMom.SetCoordinates(px0, py0, pz0, mass);
        }
        void Propagate(double s) { SetInitialState(X(s), Y(s), Z(s), Px(s), Py(s), Pz0(), Mass()); }
        [[nodiscard]] Charged CopyPropagated(double s) const { return {Index(), Charge(), fBz, XYZ(s), PxPyPzM(s)}; }
     */

   protected:
    size_t fTrackIndex;
    double fBz;  // PENDING: shouldn't it be a singleton or smth?
};

}  // namespace Tree2Secondaries

#endif  // T2S_SECONDARY_CHARGED_HXX
