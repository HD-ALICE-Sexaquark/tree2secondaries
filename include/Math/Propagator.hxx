#ifndef T2S_SECONDARY_PROPAGATOR_HXX
#define T2S_SECONDARY_PROPAGATOR_HXX

#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

#include "Math/Constants.hxx"

namespace Tree2Secondaries::Helper {

class Propagator {
   public:
    // Simple constructor for test purposes.
    explicit Propagator(double bz) : fBz{bz} {}
    ~Propagator() = default;

    // setters //
    void SetBz(double bz) { fBz = bz; }
    void SetPrimaryVertex(double x, double y, double z) { fPV.SetCoordinates(x, y, z); }

    // getters //
    double Bz() const { return fBz; }
    ROOT::Math::XYZPoint PrimaryVertex() const { return fPV; }

    // Parametrization of charged particles: (helix)                    //
    //   X(s)  = X0 + Sin(omega s) Px/omega + (1-Cos(omega s)) Py/omega //
    //   Y(s)  = Y0 - (1-Cos(omega s)) Px/omega + Sin(omega s) Py/omega //
    //   Z(s)  = Z0 + Pz s                                              //
    //   Px(s) = Cos(omega s) Px + Sin(omega s) Py                      //
    //   Py(s) = -Sin(omega s) Px + Cos(omega s) Py                     //
    //   Pz    = (constant)                                             //
    double Omega(double charge) const { return fBz * charge * Const::Kappa; }
    double ChargedX(double s, double charge, double x0, double px0, double py0) const {
        return x0 + (std::sin(Omega(charge) * s) * px0 + (1. - std::cos(Omega(charge) * s)) * py0) / Omega(charge);
    }
    double ChargedY(double s, double charge, double y0, double px0, double py0) const {
        return y0 + (-(1. - std::cos(Omega(charge) * s)) * px0 + std::sin(Omega(charge) * s) * py0) / Omega(charge);
    }
    double ChargedZ(double s, double z0, double pz0) const { return z0 + pz0 * s; }
    double ChargedPx(double s, double charge, double px0, double py0) const {
        return px0 * std::cos(Omega(charge) * s) + py0 * std::sin(Omega(charge) * s);
    }
    double ChargedPy(double s, double charge, double px0, double py0) const {
        return -px0 * std::sin(Omega(charge) * s) + py0 * std::cos(Omega(charge) * s);
    }

    // test //

    double FastChargedX(double sin_theta, double cos_theta, double one_over_omega, double x0, double px0, double py0) const {
        return x0 + (sin_theta * px0 + (1. - cos_theta) * py0) * one_over_omega;
    }
    double FastChargedY(double sin_theta, double cos_theta, double one_over_omega, double y0, double px0, double py0) const {
        return y0 + (-(1. - cos_theta) * px0 + sin_theta * py0) * one_over_omega;
    }
    double FastChargedZ(double theta, double omega, double z0, double pz0) const { return z0 + pz0 * omega * theta; }
    double FastChargedPx(double sin_theta, double cos_theta, double px0, double py0) const { return px0 * cos_theta + py0 * sin_theta; }
    double FastChargedPy(double sin_theta, double cos_theta, double px0, double py0) const { return -px0 * sin_theta + py0 * cos_theta; }

    // Parametrization of neutral particles: (line) //
    //   X(s)  = X0 + Px s                          //
    //   Y(s)  = Y0 + Py s                          //
    //   Z(s)  = Z0 + Pz s                          //
    double NeutralCoord(double s, double xi0, double pi0) const { return xi0 + pi0 * s; }

   protected:
    ROOT::Math::XYZPoint fPV{0., 0., 0.};
    double fBz{0.};
};

}  // namespace Tree2Secondaries::Helper

#endif  // T2S_SECONDARY_PROPAGATOR_HXX
