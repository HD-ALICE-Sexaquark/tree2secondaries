#ifndef T2S_SECONDARY_VERTEXER_HXX
#define T2S_SECONDARY_VERTEXER_HXX

#include "Math/Point3D.h"

#include "Math/VtxrResults.hxx"
#include "Secondary/Charged.hxx"
#include "Secondary/Neutral.hxx"

namespace Tree2Secondaries::Vertexer {

PartPartResults MinimizeDistanceHelixHelix(const Charged& q, const Charged& t);
VtxrResults MinimizeDistanceHelixVertex(const Charged& q, const ROOT::Math::XYZPoint& v);
PartPartResults MinimizeDistanceHelixLine(const Charged& q, const Neutral& n);
PartPartResults MinimizeDistanceLineLine(const Neutral& n, const Neutral& m);
VtxrResults MinimizeDistanceLineVertex(const Neutral& n, const ROOT::Math::XYZPoint& v);

}  // namespace Tree2Secondaries::Vertexer

#endif  // T2S_SECONDARY_VERTEXER_HXX
