#include "KalmanFitter/KalmanFitterParticle.hxx"

#include <cmath>

#include <Eigen/Eigen>

#include "Math/Constants.hxx"
#if T2S_LEGACY_KF
#include "Legacy/LegacyParticle.hxx"
#endif

namespace Tree2Secondaries::KalmanFitter {

// Create a `KF::Particle`, by setting `fP`, `fC` and `fQ` by hard-copying from a track.
Particle Particle::FromTrack(const View::VecTracks& v, double mass) {

    Particle out;

    out.fP(0) = v.X();
    out.fP(1) = v.Y();
    out.fP(2) = v.Z();
    out.fP(3) = v.Px();
    out.fP(4) = v.Py();
    out.fP(5) = v.Pz();
    out.fP(6) = std::sqrt(mass * mass + out.SquaredMomentum());
    out.fP(7) = 0.;

    for (unsigned int i = 0; i < 6; ++i) {
        for (unsigned int j = 0; j <= i; ++j) {
            out.fC(i, j) = v.Cov(i, j);
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
    out.fC(7, 7) = Const::Initial_Css;

    out.fQ = v.Charge<int>();

    return out;
}

// Create a `KF::Particle`, by setting `fP`, `fC` and `fQ` from a V0 view.
Particle Particle::FromV0(const View::VecV0s& v) {

    Particle out;

    out.fP(0) = v.X();
    out.fP(1) = v.Y();
    out.fP(2) = v.Z();
    out.fP(3) = v.Px();
    out.fP(4) = v.Py();
    out.fP(5) = v.Pz();
    out.fP(6) = v.Energy();
    out.fP(7) = 0.;

    for (unsigned int i = 0; i < 7; ++i) {
        for (unsigned int j = 0; j <= i; ++j) {
            out.fC(i, j) = v.Cov(i, j);
        }
    }
    out.fC(7, 7) = Const::Initial_Css;

    out.fQ = 0;

    return out;
}

#if T2S_LEGACY_KF
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

}  // namespace Tree2Secondaries::KalmanFitter
