#include <Eigen/Eigen>

#include "KalmanFitter/KalmanFitterParticle.hxx"
#include "Math/Constants.hxx"
#if T2S_LEGACY_KF
#include "Legacy/LegacyParticle.hxx"
#endif

namespace Tree2Secondaries::KalmanFitter {

// # Static Functions # //

// Create a `KF::Particle`, by setting `fP`, `fC` and `fQ` from a track view.
Particle Particle::FromTrack(const View::Rec::Track& track, double mass) {

    Particle out;

    out.fP(0) = track.X();
    out.fP(1) = track.Y();
    out.fP(2) = track.Z();
    out.fP(3) = track.Px();
    out.fP(4) = track.Py();
    out.fP(5) = track.Pz();
    out.fP(6) = std::sqrt(mass * mass + out.SquaredMomentum());
    out.fP(7) = 0.;

    out.fC(0, 0) = track.SigmaX2();
    out.fC(1, 0) = track.SigmaXY();
    out.fC(1, 1) = track.SigmaY2();
    out.fC(2, 0) = track.SigmaXZ();
    out.fC(2, 1) = track.SigmaYZ();
    out.fC(2, 2) = track.SigmaZ2();
    out.fC(3, 0) = track.SigmaXPx();
    out.fC(3, 1) = track.SigmaYPx();
    out.fC(3, 2) = track.SigmaZPx();
    out.fC(3, 3) = track.SigmaPx2();
    out.fC(4, 0) = track.SigmaXPy();
    out.fC(4, 1) = track.SigmaYPy();
    out.fC(4, 2) = track.SigmaZPy();
    out.fC(4, 3) = track.SigmaPxPy();
    out.fC(4, 4) = track.SigmaPy2();
    out.fC(5, 0) = track.SigmaXPz();
    out.fC(5, 1) = track.SigmaYPz();
    out.fC(5, 2) = track.SigmaZPz();
    out.fC(5, 3) = track.SigmaPxPz();
    out.fC(5, 4) = track.SigmaPyPz();
    out.fC(5, 5) = track.SigmaPz2();

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

    out.fQ = track.Charge();

    return out;
}

// Create a `KF::Particle`, by setting `fP`, `fC` and `fQ` from a V0 view.
Particle Particle::FromV0(const View::Rec::V0& v0) {

    Particle out;

    out.fP(0) = v0.X();
    out.fP(1) = v0.Y();
    out.fP(2) = v0.Z();
    out.fP(3) = v0.Px();
    out.fP(4) = v0.Py();
    out.fP(5) = v0.Pz();
    out.fP(6) = v0.Energy();
    out.fP(7) = 0.;

    out.fC(0, 0) = v0.SigmaX2();
    out.fC(1, 0) = v0.SigmaXY();
    out.fC(1, 1) = v0.SigmaY2();
    out.fC(2, 0) = v0.SigmaXZ();
    out.fC(2, 1) = v0.SigmaYZ();
    out.fC(2, 2) = v0.SigmaZ2();
    out.fC(3, 0) = v0.SigmaXPx();
    out.fC(3, 1) = v0.SigmaYPx();
    out.fC(3, 2) = v0.SigmaZPx();
    out.fC(3, 3) = v0.SigmaPx2();
    out.fC(4, 0) = v0.SigmaXPy();
    out.fC(4, 1) = v0.SigmaYPy();
    out.fC(4, 2) = v0.SigmaZPy();
    out.fC(4, 3) = v0.SigmaPxPy();
    out.fC(4, 4) = v0.SigmaPy2();
    out.fC(5, 0) = v0.SigmaXPz();
    out.fC(5, 1) = v0.SigmaYPz();
    out.fC(5, 2) = v0.SigmaZPz();
    out.fC(5, 3) = v0.SigmaPxPz();
    out.fC(5, 4) = v0.SigmaPyPz();
    out.fC(5, 5) = v0.SigmaPz2();

    out.fC(6, 0) = v0.SigmaXE();
    out.fC(6, 1) = v0.SigmaYE();
    out.fC(6, 2) = v0.SigmaZE();
    out.fC(6, 3) = v0.SigmaPxE();
    out.fC(6, 4) = v0.SigmaPyE();
    out.fC(6, 5) = v0.SigmaPzE();
    out.fC(6, 6) = v0.SigmaE2();

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
