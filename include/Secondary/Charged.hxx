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

    Charged(int entry, int charge, const ROOT::Math::XYZPoint &vertex, const ROOT::Math::PxPyPzMVector &momentum)
        : fState{{momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E()}, vertex},  //
          fCharge{charge},                                                              //
          fEntry{entry} {}
    // Charged(int entry, int charge, const ROOT::Math::XYZPoint &vertex, const ROOT::Math::PxPyPzEVector &momentum)
    //     : fState{momentum, vertex},  //
    //       fCharge{charge},           //
    //       fEntry{entry} {}

    [[nodiscard]] int Entry() const { return fEntry; }
    [[nodiscard]] int Charge() const { return fCharge; }

    [[nodiscard]] double X0() const { return fState.Vertex.X(); }
    [[nodiscard]] double Y0() const { return fState.Vertex.Y(); }
    [[nodiscard]] double Z0() const { return fState.Vertex.Z(); }
    [[nodiscard]] double Px0() const { return fState.Momentum.Px(); }
    [[nodiscard]] double Py0() const { return fState.Momentum.Py(); }
    [[nodiscard]] double Pz0() const { return fState.Momentum.Pz(); }
    [[nodiscard]] double Pt0() const { return fState.Momentum.Pt(); }
    [[nodiscard]] double Mass() const { return fState.Momentum.M(); }

    [[nodiscard]] double Omega(const Helper::Propagator &prop) const { return prop.Omega(fCharge); }
    [[nodiscard]] double X(double s, const Helper::Propagator &prop) const { return prop.ChargedX(s, Charge(), X0(), Px0(), Py0()); }
    [[nodiscard]] double Y(double s, const Helper::Propagator &prop) const { return prop.ChargedY(s, Charge(), Y0(), Px0(), Py0()); }
    [[nodiscard]] double Z(double s, const Helper::Propagator &prop) const { return prop.ChargedZ(s, Z0(), Pz()); }
    [[nodiscard]] ROOT::Math::XYZPoint XYZ(double s, const Helper::Propagator &prop) const { return {X(s, prop), Y(s, prop), Z(s, prop)}; }

    [[nodiscard]] double Px(double s, const Helper::Propagator &prop) const {
        //
        return prop.ChargedPx(s, Charge(), Px0(), Py0());
    }
    [[nodiscard]] double Py(double s, const Helper::Propagator &prop) const {
        //
        return prop.ChargedPy(s, Charge(), Px0(), Py0());
    }
    [[nodiscard]] double Pz() const { return Pz0(); }
    [[nodiscard]] ROOT::Math::PxPyPzEVector PxPyPzE(double s, const Helper::Propagator &prop) const {
        return {Px(s, prop), Py(s, prop), Pz(), std::sqrt(Px(s, prop) * Px(s, prop) + Py(s, prop) * Py(s, prop) + Pz() * Pz() + Mass() * Mass())};
    }
    [[nodiscard]] ROOT::Math::PxPyPzMVector PxPyPzM(double s, const Helper::Propagator &prop) const {
        return {Px(s, prop), Py(s, prop), Pz(), Mass()};
    }
    [[nodiscard]] Particle::State PropagatedState(double s, const Helper::Propagator &prop) const { return {PxPyPzE(s, prop), XYZ(s, prop)}; }

   protected:
    Particle::State fState;  // initially set arbitrary momentum
    int fCharge;
    int fEntry;
};

}  // namespace Tree2Secondaries

#endif  // T2S_SECONDARY_CHARGED_HXX
