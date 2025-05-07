#ifndef T2S_SECONDARY_RECONSTRUCTED_HXX
#define T2S_SECONDARY_RECONSTRUCTED_HXX

#include <cstdlib>
#include <memory>

#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

#include "Secondary/Particle.hxx"
#include "Secondary/True.hxx"

namespace Tree2Secondaries {

//
class Reconstructed : public Particle {
   public:
    Reconstructed(int charge, ROOT::Math::XYZPoint vtx, ROOT::Math::PxPyPzMVector momentum)
        : fMeasuredMom{std::move(momentum)}, fMeasuredVtx{std::move(vtx)}, fCharge{charge} {}
    ~Reconstructed() override = default;

    void SetCharge(int charge) { fCharge = charge; }
    void SetLinkedMc(const std::shared_ptr<True> &mc) { fLinked = mc; }

    [[nodiscard]] int Charge() const { return fCharge; }
    [[nodiscard]] const std::shared_ptr<True> &LinkedMc() const { return fLinked; }

    [[nodiscard]] virtual double X(double s) const = 0;
    [[nodiscard]] virtual double Y(double s) const = 0;
    [[nodiscard]] virtual double Z(double s) const = 0;
    [[nodiscard]] ROOT::Math::XYZPoint XYZ(double s) const { return {X(s), Y(s), Z(s)}; }
    [[nodiscard]] virtual double Px(double s) const = 0;
    [[nodiscard]] virtual double Py(double s) const = 0;
    [[nodiscard]] virtual double Pz(double s [[maybe_unused]]) const = 0;
    [[nodiscard]] ROOT::Math::PxPyPzMVector PxPyPzM(double s) const {  //
        return {Px(s), Py(s), Pz(s), Mass()};
    }
    [[nodiscard]] ROOT::Math::PxPyPzEVector PxPyPzE(double s) const {
        return {Px(s), Py(s), Pz(s), std::sqrt(Px(s) * Px(s) + Py(s) * Py(s) + Pz(s) * Pz(s) + Mass() * Mass())};
    }

    // properties at measured point //
    [[nodiscard]] ROOT::Math::XYZPoint MeasuredPos() const { return fMeasuredVtx; }
    [[nodiscard]] ROOT::Math::XYZVector PxPyPz() const { return fMeasuredMom.Vect(); }
    [[nodiscard]] ROOT::Math::PxPyPzMVector MeasuredMom() const { return fMeasuredMom; }
    [[nodiscard]] double Mass() const { return fMeasuredMom.M(); }
    [[nodiscard]] double Pt() const { return fMeasuredMom.Pt(); }
    [[nodiscard]] double Eta() const { return fMeasuredMom.Eta(); }
    [[nodiscard]] double P() const { return fMeasuredMom.Vect().R(); }
    [[nodiscard]] double P2() const { return fMeasuredMom.Vect().Mag2(); }
    [[nodiscard]] double Energy() const { return fMeasuredMom.E(); }

   protected:
    ROOT::Math::PxPyPzMVector fMeasuredMom{0., 0., 0., 0.};  // initially set arbitrary momentum
    ROOT::Math::XYZPoint fMeasuredVtx{0., 0., 0.};           // initially set arbitrary vertex
    std::shared_ptr<True> fLinked{nullptr};
    int fCharge{0};
};

}  // namespace Tree2Secondaries

#endif  // T2S_SECONDARY_RECONSTRUCTED_HXX
