#ifndef T2S_VTXR_RESULTS_HXX
#define T2S_VTXR_RESULTS_HXX

#include "Math/Point3D.h"
#include "Math/Vector4D.h"

namespace Tree2Secondaries::Particle {

struct State {
    ROOT::Math::PxPyPzEVector Momentum{0., 0., 0., 0.};
    ROOT::Math::XYZPoint Vertex{0., 0., 0.};
};

struct Pair {
    State first{{0., 0., 0., 0.}, {0., 0., 0.}};
    State second{{0., 0., 0., 0.}, {0., 0., 0.}};
};

}  // namespace Tree2Secondaries::Particle

#endif  // T2S_VTXR_RESULTS_HXX
