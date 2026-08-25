#include "KalmanFitter/KalmanFitterFitTypes.hxx"
#include "KalmanFitter/KalmanFitterMassConstraints.hxx"
#include "KalmanFitter/KalmanFitterProdVertexConstraint.hxx"
#include "KalmanFitter/KalmanFitterTransport.hxx"
#include "KalmanFitter/KalmanFitterVertexUpdate.hxx"

#include "KalmanFitter/KalmanFitter.hxx"

namespace T2DS::KF {

FitResult FitVertex(const Component& p1, const Component& p2, double bz, const FitPolicy& policy) {

    // (1) transport both legs onto their PCAs
    const Internal::TransportedPair pair = Internal::TransportToPCA(p1.part, p2.part, p1.seed, p2.seed, bz);

    // (2) vertex update, constraining the daughters' masses or not
    FitResult res = policy.pin_daughters ? Internal::VertexUpdateMC(pair) : Internal::VertexUpdate(pair);

    // (3) constrain mother's mass; after this, the daughters MUST be rescaled to kinematically close the decay tree
    if (policy.mother_mass) {
        if (const auto scale = Internal::PinMotherMass(res.mother, *policy.mother_mass)) res.RescaleDaughters(*scale);
    }

    // (4) constrain production vertex; intended to be executed after step (3)
    if (policy.prod_vertex) res.at_pv = Internal::AtProductionVertex(res.mother, *policy.prod_vertex, bz);

    return res;
}

}  // namespace T2DS::KF
