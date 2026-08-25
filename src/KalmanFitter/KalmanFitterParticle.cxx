#include <cmath>
#include <optional>

#include <Eigen/Core>

#include "common/Constants.hpp"
#include "common/DB_Particles.hpp"
#include "common/POD_PreFoundLambda.hpp"
#include "common/POD_Track.hpp"
#include "common/POD_V0.hpp"

#include "KalmanFitter/KalmanFitterUtils.hxx"

#include "KalmanFitter/KalmanFitterParticle.hxx"

namespace T2DS::KF {

namespace {

template <typename POD_LikeV0>
Particle From(const POD_LikeV0& v, const DB::Particles::Definition& pid, bool on_shell) {
    Particle out;

    out.fP(0) = static_cast<double>(v.Decay_X);
    out.fP(1) = static_cast<double>(v.Decay_Y);
    out.fP(2) = static_cast<double>(v.Decay_Z);
    out.fP(3) = static_cast<double>(v.Px);
    out.fP(4) = static_cast<double>(v.Py);
    out.fP(5) = static_cast<double>(v.Pz);
    out.fP(6) = static_cast<double>(v.Energy);
    out.fP(7) = 0.;

    UnpackSym<7>(out.fC, v.CovMatrix);
    out.fC(7, 7) = Constants::Initial_Css;
    out.fQ = 0;
    out.SetMassBookkeeping(pid, on_shell);

    return out;
}

}  // namespace

// == Constructors == //

// Create a `KF::Particle`, by setting `fP`, `fC` and `fQ` from a track.
Particle Particle::FromTrack(const POD::Track& t, const DB::Particles::Definition& pid) {

    const double mass = pid.mass;

    Particle out;

    out.fP(0) = static_cast<double>(t.X);
    out.fP(1) = static_cast<double>(t.Y);
    out.fP(2) = static_cast<double>(t.Z);
    out.fP(3) = static_cast<double>(t.Px);
    out.fP(4) = static_cast<double>(t.Py);
    out.fP(5) = static_cast<double>(t.Pz);
    out.fP(6) = std::sqrt(mass * mass + out.SquaredMomentum());
    out.fP(7) = 0.;

    UnpackSym<6>(out.fC, t.CovMatrix);

    // dE/dp_i = p_i/E linear propagation //
    // -- h = d(E)/d(px,py,pz), so the energy row is h' * C[3:6, 0:6] and its variance h' * C[3:6, 3:6] * h.
    //    `fC` being full symmetric is what lets both be written as plain blocks.
    const Eigen::Vector<double, 3> h = out.fP.segment<3>(3) / out.fP(6);

    out.fC(6, 6) = h.dot(out.fC.block<3, 3>(3, 3) * h);  // -- first: it reads the momentum block, which row 6 doesn't touch
    out.fC.block<1, 6>(6, 0).noalias() = h.transpose() * out.fC.block<3, 6>(3, 0);
    out.fC.block<6, 1>(0, 6) = out.fC.block<1, 6>(6, 0).transpose();  // -- `fC` is full symmetric, see `KF::Particle`
    out.fC(7, 7) = Constants::Initial_Css;

    out.fQ = t.Charge;

    out.SetMassBookkeeping(pid, true);

    return out;
}
Particle Particle::FromV0(const POD::V0& v, const DB::Particles::Definition& pid, bool on_shell) {  //
    return From<POD::V0>(v, pid, on_shell);
}
Particle Particle::FromPreFoundLambda(const POD::Extended::PreFoundLambda& l, const DB::Particles::Definition& pid, bool on_shell) {
    return From<POD::Extended::PreFoundLambda>(l, pid, on_shell);
}

// == Mass Constraint == //

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

// == Production Vertex Constraint == //

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
                            + (S() * S() / p2) * mom.dot(fC.block<3, 3>(3, 3) * mom);

    if (variance <= 0.) return std::nullopt;  // protection

    return std::sqrt(variance);
}

}  // namespace T2DS::KF
