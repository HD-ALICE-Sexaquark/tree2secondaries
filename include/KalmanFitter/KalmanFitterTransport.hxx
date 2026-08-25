#pragma once

#include <Eigen/Core>

#include "KalmanFitter/KalmanFitterParticle.hxx"
#include "Seeder/SeederTypes.hxx"

// ## Transport ## //

namespace T2DS::KF::Internal {

// Move a particle's (x,y,z,px,py,pz) onto a PCA. The (E,S) components are untouched: transport doesn't change
// the energy, and `S` is only ever filled by `AtProductionVertex`.
void SetStateAtPCA(Particle& p, const Seeder::PCA& pca);

// Linearised transport of a state onto a PCA:
// - `jacob` -- d(r')/d(r),       the transport Jacobian
// - `corr`  -- d(r')/d(r_other), the correlation matrix, i.e. how the *other* particle's state moved this one's `ds`
// Both are identity/zero outside their 6x6 (x,y,z,px,py,pz) corner, so only that corner is stored.
// `corr` is never materialised: the *only* way the other particle's state reaches this one is through the scalar
// `ds`, so corr = d(r')/d(ds) * d(ds)/d(r_other) is exactly rank one. Keeping its two factors instead turns every
// `corr * X * corr'` and `corr * X * jacob'` downstream into a scalar (or a row/column) times an outer product.
struct TransportOp {

    // Build the pair from a seed, given
    // - this parametrization:
    //     x' = x0 + sB * px0 + cB * py0
    //     y' = y0 - cB * px0 + sB * py0
    //     z' = z0 + dS * pz0
    //     px' = c * px0 + s * py0
    //     py' = -s * px0 + c * py0
    //     pz' = pz0,
    // - these cached `Seed` variables:
    //     sB = sin(bq * ds)/bq
    //     cB = (1 - cos(bq * ds))/bq
    //     cos = cos(bq * ds)
    //     sin = sin(bq * ds),
    // - and the d(ds)/d(r0), d(ds)/d(r_other) derivatives from `Deriv`, for (x0,y0,z0,px0,py0,pz0) respectively.
    [[nodiscard]] static TransportOp From(const Particle& p, const Seeder::Result& s, double bz);

    // C -> J C J^T. Rows and columns 6-7 (E and S) are transported by the identity, so only three of the four
    // blocks actually move; the full 8x8 triple product would spend most of its work multiplying by zeros.
    void PropagateCov(Eigen::Matrix<double, 8, 8>& c) const;

    Eigen::Matrix<double, 6, 6> jacob;  // NOTE: non-initialized on purpose, `From` fills it whole
    Eigen::Vector<double, 6> df_ds;     // d(r')/d(ds); NOTE: idem
    Eigen::Vector<double, 6> ds_dr1;    // d(ds)/d(r_other), so that corr == df_ds * ds_dr1'; NOTE: idem
};

// Both particles sitting at their PCAs, ready for the vertex update.
struct TransportedPair {
    Particle first;
    Particle second;
    // -- cross-covariance between the two transported positions ("D" in the original). It only depends on
    //    pre-transport quantities, so it is settled here and `jacob`/`corr`/the before-transport covariances
    //    never have to outlive the transport.
    Eigen::Matrix<double, 3, 3> cross;
};

// Transport both legs onto their PCAs, folding each leg's before-transport covariance into the other's.
// NOTE: the pairwise signature is the point -- each leg's `ds` depends on the other leg's state, so neither
//       transport can be applied before both operators exist.
[[nodiscard]] TransportedPair TransportToPCA(const Particle& p_1, const Particle& p_2, const Seeder::Result& s_1, const Seeder::Result& s_2,
                                             double bz);

}  // namespace T2DS::KF::Internal
