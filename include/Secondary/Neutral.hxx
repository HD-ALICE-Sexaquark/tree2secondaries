#ifndef T2S_SECONDARY_NEUTRAL_HXX
#define T2S_SECONDARY_NEUTRAL_HXX

#include <utility>

#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

#include "Math/Common.hxx"
#include "Math/Propagator.hxx"
#include "Secondary/Particle.hxx"

namespace Tree2Secondaries {

class alignas(256) Neutral {
   public:
    Neutral(const Neutral&) = delete;
    Neutral(Neutral&&) noexcept = default;
    Neutral& operator=(const Neutral&) = delete;
    Neutral& operator=(Neutral&&) = default;
    ~Neutral() = default;

    Neutral(ROOT::Math::XYZPoint vtx, ROOT::Math::PxPyPzEVector momentum) : fState{std::move(momentum), std::move(vtx)} {}

    Neutral(int neg_idx, int pos_idx, ROOT::Math::XYZPoint vtx, ROOT::Math::PxPyPzEVector momentum)
        : fState{std::move(momentum), std::move(vtx)},  //
          fNegIndex{neg_idx},
          fPosIndex{pos_idx} {}

    Neutral(int neg_idx, int pos_idx, const Particle::Pair& pair)
        : fState{pair.first.Momentum + pair.first.Momentum,                  //
                 Math::MiddlePoint(pair.first.Vertex, pair.second.Vertex)},  //
          fNegIndex{neg_idx},
          fPosIndex{pos_idx} {}

    [[nodiscard]] int NegIndex() const { return fNegIndex; }
    [[nodiscard]] ROOT::Math::XYZPoint NegVertex() const { return fNegative.Vertex; }
    [[nodiscard]] ROOT::Math::XYZPoint PosVertex() const { return fPositive.Vertex; }

    [[nodiscard]] int PosIndex() const { return fPosIndex; }
    [[nodiscard]] ROOT::Math::XYZVector NegMomentum() const { return fNegative.Momentum.Vect(); }
    [[nodiscard]] ROOT::Math::XYZVector PosMomentum() const { return fNegative.Momentum.Vect(); }

    [[nodiscard]] ROOT::Math::XYZPoint DecayVertex() const { return fState.Vertex; }
    [[nodiscard]] double DecayX() const { return fState.Vertex.X(); }
    [[nodiscard]] double DecayY() const { return fState.Vertex.Y(); }
    [[nodiscard]] double DecayZ() const { return fState.Vertex.Z(); }
    [[nodiscard]] double DecayRadius() const { return fState.Vertex.Rho(); }

    [[nodiscard]] double X(double s, const Helper::Propagator& prop) const { return prop.NeutralCoord(s, DecayX(), Px()); }
    [[nodiscard]] double Y(double s, const Helper::Propagator& prop) const { return prop.NeutralCoord(s, DecayY(), Py()); }
    [[nodiscard]] double Z(double s, const Helper::Propagator& prop) const { return prop.NeutralCoord(s, DecayZ(), Pz()); }
    [[nodiscard]] ROOT::Math::XYZPoint XYZ(double s, const Helper::Propagator& prop) const { return {X(s, prop), Y(s, prop), Z(s, prop)}; }

    [[nodiscard]] double Px() const { return fState.Momentum.Px(); }
    [[nodiscard]] double Py() const { return fState.Momentum.Py(); }
    [[nodiscard]] double Pz() const { return fState.Momentum.Pz(); }
    [[nodiscard]] ROOT::Math::XYZVector PxPyPz() const { return fState.Momentum.Vect(); }
    [[nodiscard]] ROOT::Math::PxPyPzEVector PxPyPzE() const { return fState.Momentum; }
    [[nodiscard]] double Eta() const { return fState.Momentum.Eta(); }
    [[nodiscard]] double Pt() const { return fState.Momentum.Pt(); }
    [[nodiscard]] double P2() const { return fState.Momentum.P2(); }
    [[nodiscard]] double Mass() const { return fState.Momentum.M(); }
    [[nodiscard]] double Energy() const { return fState.Momentum.E(); }

    [[nodiscard]] double CPAwrt(const ROOT::Math::XYZPoint& v) const { return Math::CosinePointingAngle(PxPyPz(), fState.Vertex, v); }
    [[nodiscard]] double DCAwrt(const ROOT::Math::XYZPoint& v) const { return Math::FastDCALineVertex(fState.Momentum.Vect(), fState.Vertex, v); }
    [[nodiscard]] double ArmenterosAlpha() const { return Math::ArmenterosAlpha(PxPyPz(), fNegative.Momentum.Vect(), fPositive.Momentum.Vect()); }
    [[nodiscard]] double ArmenterosQt() const { return Math::ArmenterosQt(PxPyPz(), fNegative.Momentum.Vect()); }
    [[nodiscard]] double DCANegWrtV0() const { return (fNegative.Vertex - DecayVertex()).R(); }
    [[nodiscard]] double DCAPosWrtV0() const { return (fPositive.Vertex - DecayVertex()).R(); }
    [[nodiscard]] double DCAbtwDaughters() const { return (fNegative.Vertex - fPositive.Vertex).R(); }
    [[nodiscard]] Particle::State PropagatedState(double s, const Helper::Propagator& prop) const { return {PxPyPzE(), XYZ(s, prop)}; }

   protected:
    Particle::State fState;
    Particle::State fNegative;
    Particle::State fPositive;
    int fNegIndex{0};
    int fPosIndex{0};
};

}  // namespace Tree2Secondaries

#endif  // T2S_SECONDARY_NEUTRAL_HXX
