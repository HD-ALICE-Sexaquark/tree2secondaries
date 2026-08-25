#include <Eigen/Core>

#include "common/Constants.hpp"

#if T2DS_DEBUG
#include "App/Logger.hxx"
#include "App/Utilities.hxx"
#endif
#include "KalmanFitter/KalmanFitterParticle.hxx"
#include "Seeder/SeederTypes.hxx"

#include "KalmanFitter/KalmanFitterTransport.hxx"

// ## Transport ## //

namespace T2DS::KF::Internal {

// Move a particle's (x,y,z,px,py,pz) onto a PCA.
void SetStateAtPCA(Particle& p, const Seeder::PCA& pca) {
    p.fP.head<3>() = Eigen::Map<const Eigen::Vector<double, 3>>(pca.xyz.data());
    p.fP.segment<3>(3) = Eigen::Map<const Eigen::Vector<double, 3>>(pca.mom.data());
}

// == TransportOp == //

TransportOp TransportOp::From(const Particle& p, const Seeder::Result& s, double bz) {

    const double bq = bz * static_cast<double>(p.Charge()) * Common::Kappa;
    const double px0 = p.Px();
    const double py0 = p.Py();
    const double pz0 = p.Pz();

    TransportOp out;  // non-initialized on purpose

    // prepare Jacobian Matrix //

    // for example, when (i,j)=(0,3)=(x',px0):
    // d(x')/d(px0) = d(sB * px0)/d(px0) + d(cB * py0)/d(px0)
    //              = sB + px0 * d(sB)/d(px0) + py0 * d(cB)/d(px0)
    //              = sB + px0 * d(sB)/d(ds) * d(ds)/d(px0) + py0 * d(cB)/d(ds) * d(ds)/d(px0)
    //              = sB + [ px0 * d(sB)/d(ds) + py0 * d(cB)/d(ds) ] * d(ds)/d(px0)
    //              = sB + d(x')/d(ds) * d(ds)/d(px0)
    //              = FC + DF_DS * DS_DR

    out.jacob.setIdentity();

    // -- first components (non-zero and non-unity)
    out.jacob(0, 3) = s.seed.sB;
    out.jacob(0, 4) = s.seed.cB;
    out.jacob(1, 3) = -s.seed.cB;
    out.jacob(1, 4) = s.seed.sB;
    out.jacob(2, 5) = s.seed.ds;
    out.jacob(3, 3) = s.seed.cos;
    out.jacob(3, 4) = s.seed.sin;
    out.jacob(4, 3) = -s.seed.sin;
    out.jacob(4, 4) = s.seed.cos;

    // -- d(r')/d(ds)
    out.df_ds(0) = s.seed.cos * px0 + s.seed.sin * py0;
    out.df_ds(1) = -s.seed.sin * px0 + s.seed.cos * py0;
    out.df_ds(2) = pz0;
    out.df_ds(3) = -bq * s.seed.sin * px0 + bq * s.seed.cos * py0;
    out.df_ds(4) = -bq * s.seed.cos * px0 - bq * s.seed.sin * py0;
    out.df_ds(5) = 0.;

    const Eigen::Map<const Eigen::RowVector<double, 6>> ds_dr(s.deriv.ds_dr.data());

    out.jacob.noalias() += out.df_ds * ds_dr;

    // prepare Correlation Matrix //

    // -- kept factored; `corr` itself is the outer product `df_ds * ds_dr1'`
    out.ds_dr1 = Eigen::Map<const Eigen::Vector<double, 6>>(s.deriv.ds_dr1.data());

#if T2DS_DEBUG
    Logger::Debug(__FUNCTION__, "jacob  = {}", out.jacob);
    Logger::Debug(__FUNCTION__, "df_ds  = {}", out.df_ds);
    Logger::Debug(__FUNCTION__, "ds_dr1 = {}", out.ds_dr1);
#endif

    return out;
}

void TransportOp::PropagateCov(Eigen::Matrix<double, 8, 8>& c) const {

    const Eigen::Matrix<double, 6, 6> c_rr = c.topLeftCorner<6, 6>();
    const Eigen::Matrix<double, 6, 2> c_re = c.topRightCorner<6, 2>();

    c.topLeftCorner<6, 6>().noalias() = jacob * c_rr * jacob.transpose();
    c.topRightCorner<6, 2>().noalias() = jacob * c_re;
    c.bottomLeftCorner<2, 6>() = c.topRightCorner<6, 2>().transpose();
    // -- `c.bottomRightCorner<2,2>()`, i.e. the (E,S) block, is transported by the identity
}

// Transport both particles onto their PCAs.
// Input arguments: (packed in a `Seeder::Result` per leg)
// - `ds`     -- transport parameter
// - `ds_dr`  -- partial derivatives of ds w.r.t. state parameters = d(ds1)/dr1
// - `ds_dr1` -- partial derivatives of current particle's ds w.r.t. other particle's state parameters = d(ds1)/dr2
// - `pca.xyz`, `pca.mom`, `theta`, `sin`, `cos`, `sB`, `cB`, `ds`
TransportedPair TransportToPCA(const Particle& p_1, const Particle& p_2, const Seeder::Result& s_1, const Seeder::Result& s_2, double bz) {

    const TransportOp op_1 = TransportOp::From(p_1, s_1, bz);
    const TransportOp op_2 = TransportOp::From(p_2, s_2, bz);

    // -- covariances before transport; nothing downstream ever looks past this corner
    const Eigen::Matrix<double, 6, 6> bt_cov_1 = p_1.fC.topLeftCorner<6, 6>();
    const Eigen::Matrix<double, 6, 6> bt_cov_2 = p_2.fC.topLeftCorner<6, 6>();

    TransportedPair out{.first = p_1, .second = p_2, .cross = Eigen::Matrix<double, 3, 3>::Zero()};

    // update states //

    SetStateAtPCA(out.first, s_1.seed.pca);
    SetStateAtPCA(out.second, s_2.seed.pca);

    // update cov matrices //

    // -- with the transport jacobians
    op_1.PropagateCov(out.first.fC);
    op_2.PropagateCov(out.second.fC);

    // -- with the correlation matrices + the other particle's before-transport cov matrix
    //    `corr * V * corr'` = (ds_dr1' * V * ds_dr1) * df_ds * df_ds', i.e. a variance of `ds` times an outer
    //    product: the other particle only ever reaches this one through the scalar `ds`
    const double var_ds_1 = op_1.ds_dr1.dot(bt_cov_2 * op_1.ds_dr1);
    const double var_ds_2 = op_2.ds_dr1.dot(bt_cov_1 * op_2.ds_dr1);

    out.first.fC.topLeftCorner<6, 6>().noalias() += var_ds_1 * op_1.df_ds * op_1.df_ds.transpose();
    out.second.fC.topLeftCorner<6, 6>().noalias() += var_ds_2 * op_2.df_ds * op_2.df_ds.transpose();

    // cross-covariance between the two transported positions //

    // -- both terms carry a rank-1 `corr` factor, so each collapses to an outer product; the parentheses keep the
    //    contraction on the vector side, where it costs a matrix-vector product instead of a matrix-matrix one
    const Eigen::RowVector<double, 3> row = op_2.ds_dr1.transpose() * bt_cov_1 * op_1.jacob.topRows<3>().transpose();
    const Eigen::Vector<double, 3> col = op_2.jacob.topRows<3>() * (bt_cov_2 * op_1.ds_dr1);

    out.cross.noalias() = op_2.df_ds.head<3>() * row;
    out.cross.noalias() += col * op_1.df_ds.head<3>().transpose();

#if T2DS_DEBUG
    Logger::Debug(__FUNCTION__, "first.fP  = {}", out.first.fP);
    Logger::Debug(__FUNCTION__, "first.fC  = {}", out.first.fC);
    Logger::Debug(__FUNCTION__, "second.fP = {}", out.second.fP);
    Logger::Debug(__FUNCTION__, "second.fC = {}", out.second.fC);
    Logger::Debug(__FUNCTION__, "cross     = {}", out.cross);
#endif

    return out;
}

}  // namespace T2DS::KF::Internal
