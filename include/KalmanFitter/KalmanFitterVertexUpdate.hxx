#pragma once

#include "KalmanFitter/KalmanFitterFitTypes.hxx"
#include "KalmanFitter/KalmanFitterTransport.hxx"

// ## Vertex Update ## //

namespace T2DS::KF::Internal {

FitResult VertexUpdate(const TransportedPair& pair);
FitResult VertexUpdateMC(const TransportedPair& pair);

}  // namespace T2DS::KF::Internal
