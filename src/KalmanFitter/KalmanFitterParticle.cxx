#include <cmath>

#include <Eigen/Eigen>

#include "common/POD_PreFoundLambda.hpp"
#include "common/POD_Track.hpp"
#include "common/POD_V0.hpp"

#include "KalmanFitter/BaseKalmanFitter.hxx"
#if T2DS_LEGACY_KF
#include "Legacy/LegacyParticle.hxx"
#endif

#include "KalmanFitter/KalmanFitterParticle.hxx"

namespace T2DS::KF {

// Create a `KF::Particle`, by setting `fP`, `fC` and `fQ` by hard-copying from a track.
Particle Particle::FromTrack(const POD::Track& v, double mass) {

    Particle out;

    out.fP(0) = v.X;
    out.fP(1) = v.Y;
    out.fP(2) = v.Z;
    out.fP(3) = v.Px;
    out.fP(4) = v.Py;
    out.fP(5) = v.Pz;
    out.fP(6) = std::sqrt(mass * mass + out.SquaredMomentum());
    out.fP(7) = 0.;

    for (unsigned int i = 0; i < 6; ++i) {
        for (unsigned int j = 0; j <= i; ++j) {
            out.fC(i, j) = v.CovMatrix[IJ(i, j)];
        }
    }

    // dE/dp_i = p_i/E linear propagation //

    double h0 = out.fP(3) / out.fP(6);
    double h1 = out.fP(4) / out.fP(6);
    double h2 = out.fP(5) / out.fP(6);

    out.fC(6, 0) = h0 * out.fC(3, 0) + h1 * out.fC(4, 0) + h2 * out.fC(5, 0);
    out.fC(6, 1) = h0 * out.fC(3, 1) + h1 * out.fC(4, 1) + h2 * out.fC(5, 1);
    out.fC(6, 2) = h0 * out.fC(3, 2) + h1 * out.fC(4, 2) + h2 * out.fC(5, 2);
    out.fC(6, 3) = h0 * out.fC(3, 3) + h1 * out.fC(4, 3) + h2 * out.fC(5, 3);
    out.fC(6, 4) = h0 * out.fC(4, 3) + h1 * out.fC(4, 4) + h2 * out.fC(5, 4);
    out.fC(6, 5) = h0 * out.fC(5, 3) + h1 * out.fC(5, 4) + h2 * out.fC(5, 5);
    out.fC(6, 6) = (h0 * h0 * out.fC(3, 3) +  //
                    h1 * h1 * out.fC(4, 4) +  //
                    h2 * h2 * out.fC(5, 5) +  //
                    2 * (h0 * h1 * out.fC(4, 3) + h0 * h2 * out.fC(5, 3) + h1 * h2 * out.fC(5, 4)));
    out.fC(7, 7) = Initial_Css;

    out.fQ = v.Charge;

    return out;
}

// Create a `KF::Particle`, by setting `fP`, `fC` and `fQ` from a V0 view.
Particle Particle::FromV0(const POD::V0& v) {

    Particle out;

    out.fP(0) = v.Decay_X;
    out.fP(1) = v.Decay_Y;
    out.fP(2) = v.Decay_Z;
    out.fP(3) = v.Px;
    out.fP(4) = v.Py;
    out.fP(5) = v.Pz;
    out.fP(6) = v.Energy;
    out.fP(7) = 0.;

    for (unsigned int i = 0; i < 7; ++i) {
        for (unsigned int j = 0; j <= i; ++j) {
            out.fC(i, j) = v.CovMatrix[IJ(i, j)];
        }
    }
    out.fC(7, 7) = Initial_Css;

    out.fQ = 0;

    return out;
}

// Create a `KF::Particle`, by setting `fP`, `fC` and `fQ` from a V0 view.
Particle Particle::FromPreFoundLambda(const POD::Extended::PreFoundLambda& l) {

    Particle out;

    out.fP(0) = l.Decay_X;
    out.fP(1) = l.Decay_Y;
    out.fP(2) = l.Decay_Z;
    out.fP(3) = l.Px;
    out.fP(4) = l.Py;
    out.fP(5) = l.Pz;
    out.fP(6) = l.Energy;
    out.fP(7) = 0.;

    for (unsigned int i = 0; i < 7; ++i) {
        for (unsigned int j = 0; j <= i; ++j) {
            out.fC(i, j) = l.CovMatrix[IJ(i, j)];
        }
    }
    out.fC(7, 7) = Initial_Css;

    out.fQ = 0;

    return out;
}

#if T2DS_LEGACY_KF
Particle Particle::FromLegacy(const Legacy::Particle& part) {

    Particle out;

    out.fP = Eigen::Map<const Eigen::Vector<double, 8>>(part.fP);

    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j <= i; ++j) {
            out.fC(i, j) = part.fC[Legacy::IJ(i, j)];
        }
    }

    out.fChi2 = part.fChi2;
    out.fNDF = part.fNDF;
    out.fQ = part.fQ;

    return out;
}
#endif

}  // namespace T2DS::KF
