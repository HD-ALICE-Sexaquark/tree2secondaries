#include <cmath>

#include <Eigen/Eigen>

#include "common/Constants.hpp"
#include "common/DB_Particles.hpp"
#include "common/POD_Event.hpp"
#include "common/POD_PreFoundLambda.hpp"
#include "common/POD_Track.hpp"
#include "common/POD_V0.hpp"

#include "KalmanFitter/BaseKalmanFitter.hxx"

#include "KalmanFitter/KalmanFitterParticle.hxx"

namespace T2DS::KF {

// Record what the fit is allowed to assume about this particle's mass.
// - An stable particle is always exactly on its own mass shell.
// - An unstable (but reconstructable) particle only carries a lower bound (the sum of its daughters' masses),
//   unless the caller asserts it was already pinned.
void Particle::SetMassBookkeeping(const DB::Particles::Definition& pid, bool on_shell) {

    if (DB::Particles::IsStable(pid)) {
        fMassHypo = pid.mass;
        fSumDaughterMass = pid.mass;
        return;
    }

    if (on_shell) fMassHypo = pid.mass;
    fSumDaughterMass = fMassHypo.value_or(DB::Particles::SumDaughterMass(pid));
}

// Create a `KF::Particle`, by setting `fP`, `fC` and `fQ` from a track.
Particle Particle::FromTrack(const POD::Track& v, const DB::Particles::Definition& pid) {

    const double mass = pid.mass;

    Particle out;

    out.fP(0) = static_cast<double>(v.X);
    out.fP(1) = static_cast<double>(v.Y);
    out.fP(2) = static_cast<double>(v.Z);
    out.fP(3) = static_cast<double>(v.Px);
    out.fP(4) = static_cast<double>(v.Py);
    out.fP(5) = static_cast<double>(v.Pz);
    out.fP(6) = std::sqrt(mass * mass + out.SquaredMomentum());
    out.fP(7) = 0.;

    for (unsigned int i = 0; i < 6; ++i) {
        for (unsigned int j = 0; j <= i; ++j) {
            out.fC(i, j) = static_cast<double>(v.CovMatrix[IJ(i, j)]);
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

    out.SetMassBookkeeping(pid, true);

    return out;
}

// Create a `KF::Particle`, by setting `fP`, `fC` and `fQ` from a V0 view.
Particle Particle::FromV0(const POD::V0& v, const DB::Particles::Definition& pid, bool on_shell) {

    Particle out;

    out.fP(0) = static_cast<double>(v.Decay_X);
    out.fP(1) = static_cast<double>(v.Decay_Y);
    out.fP(2) = static_cast<double>(v.Decay_Z);
    out.fP(3) = static_cast<double>(v.Px);
    out.fP(4) = static_cast<double>(v.Py);
    out.fP(5) = static_cast<double>(v.Pz);
    out.fP(6) = static_cast<double>(v.Energy);
    out.fP(7) = 0.;

    for (unsigned int i = 0; i < 7; ++i) {
        for (unsigned int j = 0; j <= i; ++j) {
            out.fC(i, j) = static_cast<double>(v.CovMatrix[IJ(i, j)]);
        }
    }
    out.fC(7, 7) = Initial_Css;

    out.fQ = 0;

    out.SetMassBookkeeping(pid, on_shell);

    return out;
}

// Create a `KF::Particle`, by setting `fP`, `fC` and `fQ` from a pre-found (anti)lambda view.
Particle Particle::FromPreFoundLambda(const POD::Extended::PreFoundLambda& l, const DB::Particles::Definition& pid, bool on_shell) {

    Particle out;

    out.fP(0) = static_cast<double>(l.Decay_X);
    out.fP(1) = static_cast<double>(l.Decay_Y);
    out.fP(2) = static_cast<double>(l.Decay_Z);
    out.fP(3) = static_cast<double>(l.Px);
    out.fP(4) = static_cast<double>(l.Py);
    out.fP(5) = static_cast<double>(l.Pz);
    out.fP(6) = static_cast<double>(l.Energy);
    out.fP(7) = 0.;

    for (unsigned int i = 0; i < 7; ++i) {
        for (unsigned int j = 0; j <= i; ++j) {
            out.fC(i, j) = static_cast<double>(l.CovMatrix[IJ(i, j)]);
        }
    }
    out.fC(7, 7) = Initial_Css;

    out.fQ = 0;

    out.SetMassBookkeeping(pid, on_shell);

    return out;
}

// Propagate the error of the decay length l = s * p, with s = `fP(7)`.
//     d(l)/d(s) = p
//     d(l)/d(p_i) = s * p_i / p
//  => var(l) = p^2 * var(s) + 2 * s * sum_i p_i * cov(s,p_i) + (s/p)^2 * p^T * cov(p,p) * p
std::optional<double> Particle::DecayLengthErr() const {

    const double p2 = SquaredMomentum();
    if (p2 < Common::AbsAlmostZero) return std::nullopt;  // protection

    const Eigen::Vector<double, 3> mom = fP.segment<3>(3);

    const double variance = p2 * VarS()                                                         //
                            + 2. * S() * (Px() * CovSPx() + Py() * CovSPy() + Pz() * CovSPz())  //
                            + (S() * S() / p2) * mom.dot(fC.block<3, 3>(3, 3).selfadjointView<Eigen::Lower>() * mom);

    if (variance <= 0.) return std::nullopt;  // protection

    return std::sqrt(variance);
}

// Create a `KF::Vertex` from the event's reconstructed primary vertex.
Vertex Vertex::FromEvent(const POD::Event& e) {

    Vertex out;

    out.xyz(0) = static_cast<double>(e.PV_X);
    out.xyz(1) = static_cast<double>(e.PV_Y);
    out.xyz(2) = static_cast<double>(e.PV_Z);

    for (unsigned int i = 0; i < 3; ++i) {
        for (unsigned int j = 0; j <= i; ++j) {
            out.cov(i, j) = static_cast<double>(e.PV_CovMatrix[IJ(i, j)]);
            out.cov(j, i) = out.cov(i, j);
        }
    }

    return out;
}

}  // namespace T2DS::KF
