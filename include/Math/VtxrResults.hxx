#ifndef T2S_VTXR_RESULTS_HXX
#define T2S_VTXR_RESULTS_HXX

#include "Math/Point3D.h"

namespace Tree2Secondaries {

struct VtxrResults {
    double ds{0.};
    ROOT::Math::XYZPoint pca{0., 0., 0.};
};

struct PartPartResults {
    VtxrResults q{0., {0., 0., 0.}};
    VtxrResults t{0., {0., 0., 0.}};
};

}  // namespace Tree2Secondaries

#endif  // T2S_VTXR_RESULTS_HXX
