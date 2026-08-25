#pragma once

#include "KalmanFitter/KalmanFitterParticle.hxx"
#include "KalmanFitter/KalmanFitterVertex.hxx"

// ## Production Vertex Constraint ## //

namespace T2DS::KF::Internal {

[[nodiscard]] Particle AtProductionVertex(const Particle& part, const Vertex& vtx, double bz);

}  // namespace T2DS::KF::Internal
