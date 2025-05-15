#ifndef T2S_SECONDARY_NEUTRAL_HXX
#define T2S_SECONDARY_NEUTRAL_HXX

#include <utility>

#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

#include "Math/Common.hxx"
#include "Math/Propagator.hxx"
#include "Secondary/Particle.hxx"
#include "Utilities/Logger.hxx"

namespace Tree2Secondaries {

class Neutral {
   public:
    Neutral(const Neutral&) = delete;
    Neutral(Neutral&&) noexcept = default;
    Neutral& operator=(const Neutral&) = delete;
    Neutral& operator=(Neutral&&) = default;
    ~Neutral() = default;

    Neutral(ROOT::Math::XYZPoint vtx, ROOT::Math::PxPyPzEVector momentum) : fState{std::move(momentum), std::move(vtx)} {}

    Neutral(int entry, int neg_entry, int pos_entry, int hypothesis_pid, const Particle::Pair& pair)
        : fState{pair.first.Momentum + pair.second.Momentum,                 //
                 Math::MiddlePoint(pair.first.Vertex, pair.second.Vertex)},  //
          fNegative(pair.first),
          fPositive(pair.second),
          fEntry{entry},
          fNegEntry{neg_entry},
          fPosEntry{pos_entry},
          fHypothesisPID{hypothesis_pid} {}

    [[nodiscard]] int Entry() const { return fEntry; }
    [[nodiscard]] int NegEntry() const { return fNegEntry; }
    [[nodiscard]] int PosEntry() const { return fPosEntry; }
    [[nodiscard]] int HypothesisPID() const { return fHypothesisPID; }

    [[nodiscard]] ROOT::Math::XYZPoint NegVertex() const { return fNegative.Vertex; }
    [[nodiscard]] ROOT::Math::XYZPoint PosVertex() const { return fPositive.Vertex; }

    [[nodiscard]] ROOT::Math::XYZVector NegMomentum() const { return fNegative.Momentum.Vect(); }
    [[nodiscard]] ROOT::Math::XYZVector PosMomentum() const { return fPositive.Momentum.Vect(); }

    [[nodiscard]] ROOT::Math::PxPyPzEVector NegPxPyPzE() const { return fNegative.Momentum; }
    [[nodiscard]] ROOT::Math::PxPyPzEVector PosPxPyPzE() const { return fPositive.Momentum; }

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

    void Print(int pdg_code_v0, const Helper::Propagator& prop) {
        INFO("Searching %i {%i,%i},P={%f,%f,%f,m=%f,e=%f},V={%f,%f,%f}", pdg_code_v0, fNegEntry, fPosEntry, fState.Momentum.Px(),
             fState.Momentum.Py(), fState.Momentum.Pz(), fState.Momentum.M(), fState.Momentum.E(), fState.Vertex.X(), fState.Vertex.Y(),
             fState.Vertex.Z());
        fNegative.Print();
        fPositive.Print();
        INFO("m%f dbd%f z%f r%f dn%f dp%f pt%f et%f qt%f a%f cpv%f dpv%f", Mass(), DCAbtwDaughters(), DecayZ(), DecayRadius(), DCANegWrtV0(),
             DCAPosWrtV0(), Pt(), Eta(), ArmenterosQt(), ArmenterosAlpha(), CPAwrt(prop.PrimaryVertex()), DCAwrt(prop.PrimaryVertex()));
    }

   protected:
    Particle::State fState;
    Particle::State fNegative{{0., 0., 0., 0.}, {0., 0., 0.}};
    Particle::State fPositive{{0., 0., 0., 0.}, {0., 0., 0.}};
    int fEntry{0};
    int fNegEntry{0};
    int fPosEntry{0};
    int fHypothesisPID{0};
};

}  // namespace Tree2Secondaries

#endif  // T2S_SECONDARY_NEUTRAL_HXX
