#ifndef T2S_VTXR_RESULTS_HXX
#define T2S_VTXR_RESULTS_HXX

#include "Math/Point3D.h"
#include "Math/Vector4D.h"

#include "App/Logger.hxx"

namespace Tree2Secondaries::Particle {

struct State {
    ROOT::Math::PxPyPzEVector Momentum{0., 0., 0., 0.};
    ROOT::Math::XYZPoint Vertex{0., 0., 0.};
    void Print() {
        INFO("P={%f,%f,%f,m=%f,e=%f},V={%f,%f,%f}", Momentum.Px(), Momentum.Py(), Momentum.Pz(), Momentum.M(), Momentum.E(), Vertex.X(), Vertex.Y(),
             Vertex.Z());
    }
};

struct Pair {
    State first{{0., 0., 0., 0.}, {0., 0., 0.}};
    State second{{0., 0., 0., 0.}, {0., 0., 0.}};
};

}  // namespace Tree2Secondaries::Particle

#endif  // T2S_VTXR_RESULTS_HXX
