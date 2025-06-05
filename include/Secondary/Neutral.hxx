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

    int Entry() const { return fEntry; }
    int NegEntry() const { return fNegEntry; }
    int PosEntry() const { return fPosEntry; }
    int HypothesisPID() const { return fHypothesisPID; }

    ROOT::Math::XYZPoint NegVertex() const { return fNegative.Vertex; }
    ROOT::Math::XYZPoint PosVertex() const { return fPositive.Vertex; }

    ROOT::Math::XYZVector NegMomentum() const { return fNegative.Momentum.Vect(); }
    ROOT::Math::XYZVector PosMomentum() const { return fPositive.Momentum.Vect(); }

    ROOT::Math::PxPyPzEVector NegPxPyPzE() const { return fNegative.Momentum; }
    ROOT::Math::PxPyPzEVector PosPxPyPzE() const { return fPositive.Momentum; }

    ROOT::Math::XYZPoint DecayVertex() const { return fState.Vertex; }
    double DecayX() const { return fState.Vertex.X(); }
    double DecayY() const { return fState.Vertex.Y(); }
    double DecayZ() const { return fState.Vertex.Z(); }
    double DecayRadius() const { return fState.Vertex.Rho(); }

    double X(double s, const Helper::Propagator& prop) const { return prop.NeutralCoord(s, DecayX(), Px()); }
    double Y(double s, const Helper::Propagator& prop) const { return prop.NeutralCoord(s, DecayY(), Py()); }
    double Z(double s, const Helper::Propagator& prop) const { return prop.NeutralCoord(s, DecayZ(), Pz()); }
    ROOT::Math::XYZPoint XYZ(double s, const Helper::Propagator& prop) const { return {X(s, prop), Y(s, prop), Z(s, prop)}; }

    double Px() const { return fState.Momentum.Px(); }
    double Py() const { return fState.Momentum.Py(); }
    double Pz() const { return fState.Momentum.Pz(); }
    ROOT::Math::XYZVector PxPyPz() const { return fState.Momentum.Vect(); }
    ROOT::Math::PxPyPzEVector PxPyPzE() const { return fState.Momentum; }
    double Eta() const { return fState.Momentum.Eta(); }
    double Pt() const { return fState.Momentum.Pt(); }
    double P2() const { return fState.Momentum.P2(); }
    double Mass() const { return fState.Momentum.M(); }
    double Energy() const { return fState.Momentum.E(); }

    double CPAwrt(const ROOT::Math::XYZPoint& v) const { return Math::CosinePointingAngle(PxPyPz(), fState.Vertex, v); }
    double DCAwrt(const ROOT::Math::XYZPoint& v) const { return Math::FastDCALineVertex(fState.Momentum.Vect(), fState.Vertex, v); }
    double ArmenterosAlpha() const { return Math::ArmenterosAlpha(PxPyPz(), fNegative.Momentum.Vect(), fPositive.Momentum.Vect()); }
    double ArmenterosQt() const { return Math::ArmenterosQt(PxPyPz(), fNegative.Momentum.Vect()); }
    double DCANegWrtV0() const { return (fNegative.Vertex - DecayVertex()).R(); }
    double DCAPosWrtV0() const { return (fPositive.Vertex - DecayVertex()).R(); }
    double DCAbtwDaughters() const { return (fNegative.Vertex - fPositive.Vertex).R(); }
    Particle::State PropagatedState(double s, const Helper::Propagator& prop) const { return {PxPyPzE(), XYZ(s, prop)}; }

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
