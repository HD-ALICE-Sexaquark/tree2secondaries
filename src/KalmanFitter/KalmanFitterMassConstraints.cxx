#include <cmath>
#include <optional>

#include <Eigen/Core>

#include "KalmanFitter/KalmanFitterMassConstraints.hxx"

// ## Mass Constraints ## //

namespace T2DS::KF::Internal {

namespace {

constexpr int MaxIter = 10;  // PENDING: maybe it's already done in 2-3 steps
constexpr double Tolerance = 1.E-10;
constexpr double MinDenom = 1.E-10;
constexpr double MinVariance = 1.E-20;

// The outcome of a mass-shell pin: the two rescaling factors, and d(px,py,pz,E)'/d(px,py,pz,E). The position rows
// of the full 7x7 Jacobian are the identity and its off-diagonal blocks are zero, so only this corner is kept.
struct MassPin {
    MassScale scale;
    Eigen::Matrix<double, 4, 4> jacob;
};

// Pin the state vector `p` (and propagate its covariance `c`) onto the mass shell defined by `mass`, by rescaling
//     p -> p/(1 - lambda)
//     E -> E/(1 + lambda)
// where lambda is chosen so that E'^2 - p'^2 = mass^2.
// Substituting gives the quartic:
//     f(lambda) = -m^2 * lambda^4 + a * lambda^2 + b * lambda + c
//     a = E^2 - p^2 + 2m^2
//     b = -2(E^2 + p^2)
//     c = E^2 - p^2 - m^2
// whose root lambda = 0 corresponds to a particle already on shell, since f(0) = c.
// If no usable lambda could be found, returns nullopt; leaving `p` and `c` untouched.
[[nodiscard]] std::optional<MassPin> PinToMassShell(Eigen::Vector<double, 8>& p, Eigen::Matrix<double, 8, 8>& c, double mass) {

    // a negative energy cannot be rescaled onto a physical shell
    if (p(6) < 0.) return std::nullopt;

    const double energy2 = p(6) * p(6);
    const double mom2 = p.segment<3>(3).squaredNorm();
    const double mass2 = mass * mass;

    const double a = energy2 - mom2 + 2. * mass2;
    const double b = -2. * (energy2 + mom2);
    const double c0 = energy2 - mom2 - mass2;

    // solve for lambda //

    // -- seed with the smaller root of the quadratic part, i.e. dropping the -m^2*lambda^4 term.
    //    (energy2 + mom2) == -b/2 and d == (b^2 - 4ac)/4, so the textbook form is (energy2 + mom2 - sqrt(d))/a.
    //    That subtracts two nearly-equal numbers exactly in the near-on-shell case that dominates here, so use
    //    the conjugate instead: the roots multiply to c0/a, hence lambda_- == c0/((energy2 + mom2) + sqrt(d)).
    //    It also degrades gracefully to the linear root -c0/b as a -> 0, so no separate fallback is needed.
    const double d = 4. * energy2 * mom2 - mass2 * (energy2 - mom2 - 2. * mass2);
    const double q = energy2 + mom2 + (d > 0. ? std::sqrt(d) : 0.);

    double lambda = 0.;
    if (q > MinDenom) lambda = c0 / q;

    // -- refine by Newton
    for (int i = 0; i < MaxIter; ++i) {
        const double lambda2 = lambda * lambda;
        const double f = -mass2 * lambda2 * lambda2 + a * lambda2 + b * lambda + c0;
        const double df = -4. * mass2 * lambda2 * lambda + 2. * a * lambda + b;
        if (std::abs(df) < MinDenom) break;
        const double delta = f / df;
        lambda -= delta;
        if (std::abs(delta) < Tolerance) break;
    }

    // -- whatever the loop settled on is applied, as in the original: an iterate that missed the tolerance is still
    //    closer to the shell than not pinning at all, and refusing it would silently drop the
    //    mass(mother) >= sum(mass(daughters)) guarantee that the pin exists to provide

    // -- protection: lambda = +-1 would blow up the rescaling; leave caller's state untouched
    if (!std::isfinite(lambda) || std::abs(1. - lambda) < MinDenom || std::abs(1. + lambda) < MinDenom) return std::nullopt;

    const double lpi = 1. / (1. + lambda);
    const double lmi = 1. / (1. - lambda);
    const double lp2i = lpi * lpi;
    const double lm2i = lmi * lmi;

    // prepare Jacobian Matrix //

    // -- d(lambda)/d(px,py,pz,E), by implicit differentiation of f
    const double lambda2 = lambda * lambda;
    const double dfl = -4. * mass2 * lambda2 * lambda + 2. * a * lambda + b;

    // protection
    if (std::abs(dfl) < MinDenom) return std::nullopt;

    Eigen::Vector<double, 4> dfx;
    dfx(0) = -2. * (1. + lambda) * (1. + lambda) * p(3);
    dfx(1) = -2. * (1. + lambda) * (1. + lambda) * p(4);
    dfx(2) = -2. * (1. + lambda) * (1. + lambda) * p(5);
    dfx(3) = 2. * (1. - lambda) * (1. - lambda) * p(6);

    const Eigen::Vector<double, 4> dlx = -dfx / dfl;

    // -- d(px',py',pz',E')/d(lambda)
    Eigen::Vector<double, 4> dxx;
    dxx(0) = p(3) * lm2i;
    dxx(1) = p(4) * lm2i;
    dxx(2) = p(5) * lm2i;
    dxx(3) = -p(6) * lp2i;

    Eigen::Matrix<double, 4, 4> j;
    j.noalias() = dxx * dlx.transpose();
    j.topLeftCorner<3, 3>().diagonal().array() += lmi;
    j(3, 3) += lpi;

    // apply //

    // -- the position rows of the full 7x7 Jacobian are the identity and its off-diagonal blocks are zero, so
    //    `c`'s position block comes out untouched and only three blocks actually move. That is what keeps the
    //    later cross correction in `VertexUpdateMC` valid.
    // -- row/col 7 (S) is left alone, matching the original
    //    NOTE: that leaves S's cross-covariances with the rescaled momenta stale. Harmless as long as nothing reads
    //          them -- `fP(7)` is never even filled -- but it is a real inconsistency, not just an omission.
    const Eigen::Matrix<double, 4, 3> c_px = c.block<4, 3>(3, 0);
    const Eigen::Matrix<double, 4, 4> c_pp = c.block<4, 4>(3, 3);

    c.block<4, 3>(3, 0).noalias() = j * c_px;
    c.block<3, 4>(0, 3) = c.block<4, 3>(3, 0).transpose();

    // -- the lower triangle is authoritative, exactly as in the original, whose packed covariance has no upper
    //    half that could disagree with it
    const Eigen::Matrix<double, 4, 4> c_pp_new = j * c_pp * j.transpose();
    c.block<4, 4>(3, 3) = c_pp_new.selfadjointView<Eigen::Lower>();

    p.segment<3>(3) *= lmi;
    p(6) *= lpi;

    return MassPin{.scale = MassScale{.p_scale = lmi, .e_scale = lpi}, .jacob = j};
}

}  // namespace

std::optional<Eigen::Matrix<double, 4, 4>> PinDaughterToMassShell(Eigen::Vector<double, 8>& p, Eigen::Matrix<double, 8, 8>& c,
                                                                  const std::optional<double>& mass_hypo, double sum_daughter_mass) {

    if (mass_hypo) {
        if (const auto pin = PinToMassShell(p, c, *mass_hypo)) return pin->jacob;
        return std::nullopt;
    }

    const double m2 = p(6) * p(6) - p.segment<3>(3).squaredNorm();
    if (m2 < sum_daughter_mass * sum_daughter_mass) {
        if (const auto pin = PinToMassShell(p, c, sum_daughter_mass)) return pin->jacob;
    }

    return std::nullopt;
}

std::optional<MassScale> PinMotherMass(Particle& part, double mass) {

    // h = d(m^2)/d(px,py,pz,E)
    Eigen::Vector<double, 4> h;
    h << -2. * part.Px(), -2. * part.Py(), -2. * part.Pz(), 2. * part.E();

    const double var_m2 = h.transpose() * part.fC.block<4, 4>(3, 3) * h;

    // -- the mass error is already 0, so the particle can't be constrained (the original doesn't guard this one,
    //    but its linearised counterpart does)
    if (var_m2 < MinVariance) return std::nullopt;

    const double residual = part.E() * part.E() - part.SquaredMomentum() - mass * mass;

    // -- the pin's own Jacobian is discarded: there's no cross-covariance left to rotate at this point
    const std::optional<MassPin> pin = PinToMassShell(part.fP, part.fC, mass);
    if (!pin) return std::nullopt;

    // -- both bail-outs above happen before anything is mutated, so a failure leaves `part` untouched

    part.fChi2 += residual * residual / var_m2;
    part.fNDF += 1;
    part.fMassHypo = mass;
    part.fSumDaughterMass = mass;

    return pin->scale;
}

}  // namespace T2DS::KF::Internal
