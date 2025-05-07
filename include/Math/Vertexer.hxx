#ifndef T2S_SECONDARY_VERTEXER_HXX
#define T2S_SECONDARY_VERTEXER_HXX

#include "Math/Point3D.h"

#include "Math/Propagator.hxx"
#include "Secondary/Charged.hxx"
#include "Secondary/Neutral.hxx"
#include "Secondary/Particle.hxx"

namespace Tree2Secondaries::Vertexer {

Particle::Pair MinimizeDistanceHelixHelix(const Charged& q, const Charged& t, const Helper::Propagator& prop);
Particle::State MinimizeDistanceHelixVertex(const Charged& q, const ROOT::Math::XYZPoint& v, const Helper::Propagator& prop);
Particle::Pair MinimizeDistanceHelixLine(const Charged& q, const Neutral& n, const Helper::Propagator& prop);
Particle::Pair MinimizeDistanceLineLine(const Neutral& n, const Neutral& m, const Helper::Propagator& prop);
Particle::State MinimizeDistanceLineVertex(const Neutral& n, const ROOT::Math::XYZPoint& v, const Helper::Propagator& prop);

}  // namespace Tree2Secondaries::Vertexer

#endif  // T2S_SECONDARY_VERTEXER_HXX
