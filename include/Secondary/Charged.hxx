#ifndef T2S_SECONDARY_CHARGED_HXX
#define T2S_SECONDARY_CHARGED_HXX

#include "Math/Point3D.h"
#include "Math/Vector4D.h"

#include "Math/Propagator.hxx"
#include "Secondary/Particle.hxx"

namespace Tree2Secondaries {

class Charged {
   public:
    Charged(const Charged &) = delete;
    Charged(Charged &&) noexcept = default;
    Charged &operator=(const Charged &) = delete;
    Charged &operator=(Charged &&) = default;
    ~Charged() = default;

    Charged(int entry, int charge, const ROOT::Math::XYZPoint &vertex, const ROOT::Math::PxPyPzMVector &pxpypzm)
        : fState{{pxpypzm.Px(), pxpypzm.Py(), pxpypzm.Pz(), pxpypzm.E()}, vertex},  //
          fCharge{charge},                                                          //
          fEntry{entry} {}
    // Charged(int entry, int charge, const ROOT::Math::XYZPoint &vertex, const ROOT::Math::PxPyPzEVector &momentum)
    //     : fState{momentum, vertex},  //
    //       fCharge{charge},           //
    //       fEntry{entry} {}

    int Entry() const { return fEntry; }
    int Charge() const { return fCharge; }

    double X0() const { return fState.Vertex.X(); }
    double Y0() const { return fState.Vertex.Y(); }
    double Z0() const { return fState.Vertex.Z(); }
    double Px0() const { return fState.Momentum.Px(); }
    double Py0() const { return fState.Momentum.Py(); }
    double Pz0() const { return fState.Momentum.Pz(); }
    double Pt() const { return fState.Momentum.Pt(); }
    double Mass() const { return fState.Momentum.M(); }
    double Energy() const { return fState.Momentum.E(); }

    double Omega(const Helper::Propagator &prop) const { return prop.Omega(fCharge); }
    double X(double s, const Helper::Propagator &prop) const { return prop.ChargedX(s, Charge(), X0(), Px0(), Py0()); }
    double Y(double s, const Helper::Propagator &prop) const { return prop.ChargedY(s, Charge(), Y0(), Px0(), Py0()); }
    double Z(double s, const Helper::Propagator &prop) const { return prop.ChargedZ(s, Z0(), Pz()); }
    ROOT::Math::XYZPoint XYZ(double s, const Helper::Propagator &prop) const { return {X(s, prop), Y(s, prop), Z(s, prop)}; }
    ROOT::Math::XYZPoint FastXYZ(double theta, double sin_theta, double cos_theta, double omega, double one_over_omega,
                                 const Helper::Propagator &prop) const {
        return {prop.FastChargedX(sin_theta, cos_theta, one_over_omega, X0(), Px0(), Py0()),  //
                prop.FastChargedY(sin_theta, cos_theta, one_over_omega, Y0(), Px0(), Py0()),  //
                prop.FastChargedZ(theta, omega, Z0(), Pz0())};
    }

    double Px(double s, const Helper::Propagator &prop) const { return prop.ChargedPx(s, Charge(), Px0(), Py0()); }
    double Py(double s, const Helper::Propagator &prop) const { return prop.ChargedPy(s, Charge(), Px0(), Py0()); }
    double FastPx(double sin_theta, double cos_theta, const Helper::Propagator &prop) const {
        return prop.FastChargedPx(sin_theta, cos_theta, Px0(), Py0());
    }
    double FastPy(double sin_theta, double cos_theta, const Helper::Propagator &prop) const {
        return prop.FastChargedPy(sin_theta, cos_theta, Px0(), Py0());
    }
    double Pz() const { return Pz0(); }
    ROOT::Math::PxPyPzEVector PxPyPzE(double s, const Helper::Propagator &prop) const {
        return {Px(s, prop), Py(s, prop), Pz(), std::sqrt(Px(s, prop) * Px(s, prop) + Py(s, prop) * Py(s, prop) + Pz() * Pz() + Mass() * Mass())};
    }
    ROOT::Math::PxPyPzMVector PxPyPzM(double s, const Helper::Propagator &prop) const { return {Px(s, prop), Py(s, prop), Pz(), Mass()}; }
    Particle::State PropagatedState(double s, const Helper::Propagator &prop) const { return {PxPyPzE(s, prop), XYZ(s, prop)}; }

   protected:
    Particle::State fState;  // initially set arbitrary momentum
    int fCharge;
    int fEntry;
};

}  // namespace Tree2Secondaries

#endif  // T2S_SECONDARY_CHARGED_HXX
