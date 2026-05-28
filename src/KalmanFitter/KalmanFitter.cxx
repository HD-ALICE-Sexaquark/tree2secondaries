#include "KalmanFitter/BaseKalmanFitter.hxx"

#include <Eigen/Eigen>

#include "common/Constants.hpp"

#include "Seeder/BaseSeeder.hxx"
#if R2DS_DEBUG
#include "App/Logger.hxx"
#include "App/Utilities.hxx"
#endif

namespace R2DS::KalmanFitter {

// Daughter Struct //

// Calculate the transport Jacobian = d(fP new)/d(fP old) and the correlation matrix = d(fP new)/d(r1),
// given
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
// - and the d(ds)/d(r0) derivatives from `Deriv` for (x0,y0,z0,px0,py0,z0), respectively.
void Daughter::PrepareJacobAndCorr(const Seeder::Result& s, double bz) {

    // get initial values //

    double bq = bz * static_cast<double>(Charge()) * Common::Kappa;
    double px0 = Px();
    double py0 = Py();
    double pz0 = Pz();

    // prepare Jacobian Matrix //

    jacob = Eigen::Matrix<double, 8, 8>::Identity();

    // for example, when (i,j)=(0,3)=(x',px0):
    // d(x')/d(px0) = d(sB * px0)/d(px0) + d(cB * py0)/d(px0)
    //              = sB + px0 * d(sB)/d(px0) + py0 * d(cB)/d(px0)
    //              = sB + px0 * d(sB)/d(ds) * d(ds)/d(px0) + py0 * d(cB)/d(ds) * d(ds)/d(px0)
    //              = sB + [ px0 * d(sB)/d(ds) + py0 * d(cB)/d(ds) ] * d(ds)/d(px0)
    //              = sB + d(x')/d(ds) * d(ds)/d(px0)
    //              = FC + DF_DS * DS_DR

    // -- first components (non-zero and non-unity)
    jacob(0, 3) = s.seed.sB;
    jacob(0, 4) = s.seed.cB;
    jacob(1, 3) = -s.seed.cB;
    jacob(1, 4) = s.seed.sB;
    jacob(2, 5) = s.seed.ds;
    jacob(3, 3) = s.seed.cos;
    jacob(3, 4) = s.seed.sin;
    jacob(4, 3) = -s.seed.sin;
    jacob(4, 4) = s.seed.cos;

    // -- d(r')/d(ds)
    Eigen::Vector<double, 6> df_ds;
    df_ds(0) = s.seed.cos * px0 + s.seed.sin * py0;
    df_ds(1) = -s.seed.sin * px0 + s.seed.cos * py0;
    df_ds(2) = pz0;
    df_ds(3) = -bq * s.seed.sin * px0 + bq * s.seed.cos * py0;
    df_ds(4) = -bq * s.seed.cos * px0 - bq * s.seed.sin * py0;
    df_ds(5) = 0.;

    Eigen::RowVector<double, 6> ds_dr = Eigen::Map<const Eigen::RowVector<double, 6>>(s.deriv.ds_dr.data());

    jacob.block<6, 6>(0, 0).noalias() += df_ds * ds_dr;

    // prepare Correlation Matrix //

    Eigen::RowVector<double, 6> ds_dr1 = Eigen::Map<const Eigen::RowVector<double, 6>>(s.deriv.ds_dr1.data());

    corr = df_ds * ds_dr1;

#if R2DS_DEBUG
    Logger::Debug(__FUNCTION__, "jacob  = {}", jacob);
    Logger::Debug(__FUNCTION__, "corr   = {}", corr);
#endif
}

// Transport particle.
// Input arguments: (packed in a single `Seeder::Seed` struct)
// - `ds`     -- transport parameter
// - `ds_dr`  -- partial derivatives of ds w.r.t. state parameters = d(ds1)/dr1
// - `ds_dr1` -- partial derivatives of current particle's ds w.r.t. other particle's state parameters = d(ds1)/dr2
// - `pca.xyz`, `pca.mom`, `theta`, `sin`, `cos`, `sB`, `cB`, `ds`
void Daughter::Transport(const Seeder::PCA& pca, const Eigen::Ref<const Eigen::Matrix<double, 8, 8>>& other_bt_cov) {

    // update state //

    fP(0) = pca.xyz[0];
    fP(1) = pca.xyz[1];
    fP(2) = pca.xyz[2];
    fP(3) = pca.mom[0];
    fP(4) = pca.mom[1];
    fP(5) = pca.mom[2];

    // update cov matrix //

    // -- with jacobian
    fC = jacob * fC.selfadjointView<Eigen::Lower>() * jacob.transpose();

    // -- with corr. matrix + other particle's cov matrix
    fC.block<6, 6>(0, 0).noalias() += corr * other_bt_cov.block<6, 6>(0, 0).selfadjointView<Eigen::Lower>() * corr.transpose();

#if R2DS_DEBUG
    Logger::Debug(__FUNCTION__, "fP = {}", fP);
    Logger::Debug(__FUNCTION__, "fC = {}", fC);
#endif
}

// Main Fitting Method //

Particle GetUpdated(const Daughter& kf_1, const Daughter& kf_2) {

    Particle out = kf_1;  // PENDING: discarded data?

    // sum of position covariances //

    Eigen::Matrix<double, 3, 3> mS = kf_1.fC.block<3, 3>(0, 0) + kf_2.fC.block<3, 3>(0, 0);
    auto ldlt = mS.selfadjointView<Eigen::Lower>().ldlt();

    // residual = measured - estimated //

    Eigen::Vector<double, 3> zeta = kf_2.fP.head<3>() - kf_1.fP.head<3>();

    // update chi2 //

    out.fChi2 = zeta.transpose() * ldlt.solve(zeta);

    // correlation between state and position measurement //

    Eigen::Matrix<double, 7, 3> mCHt = kf_1.fC.block<7, 3>(0, 0);
    mCHt.block<4, 3>(3, 0) -= kf_2.fC.block<4, 3>(3, 0);

    // Kalman gain //

    Eigen::Matrix<double, 7, 3> mK = ldlt.solve(mCHt.transpose()).transpose();

    // add 4-momentum //

    out.fP.segment<4>(3).noalias() += kf_2.fP.segment<4>(3);
    out.fC.block<4, 4>(3, 3).noalias() += kf_2.fC.block<4, 4>(3, 3);

    // state update: P += K * zeta //

    out.fP.head<7>().noalias() += mK * zeta;

    // covariance update: C -= K * CHt' //

    out.fC.block<7, 7>(0, 0).noalias() -= mK * mCHt.transpose();

    // correlation correction //

    Eigen::Matrix<double, 3, 3> F3C1F1T =
        kf_2.corr.block<3, 6>(0, 0) * kf_1.bt_cov.block<6, 6>(0, 0).selfadjointView<Eigen::Lower>() * kf_1.jacob.block<3, 6>(0, 0).transpose();
    Eigen::Matrix<double, 3, 3> F4C2F2T =
        kf_2.jacob.block<3, 6>(0, 0) * kf_2.bt_cov.block<6, 6>(0, 0).selfadjointView<Eigen::Lower>() * kf_1.corr.block<3, 6>(0, 0).transpose();
    Eigen::Matrix<double, 3, 3> D = F3C1F1T + F4C2F2T;

    Eigen::Matrix<double, 3, 3> K = ldlt.solve(kf_1.fC.block<3, 3>(0, 0)).transpose();
    Eigen::Matrix<double, 3, 3> K2 = Eigen::Matrix<double, 3, 3>::Identity() - K.transpose();
    Eigen::Matrix<double, 3, 3> A = D * K2;
    Eigen::Matrix<double, 3, 3> M = K * A;

    out.fC.block<3, 3>(0, 0).noalias() += M + M.transpose();

    // update charge and NDF //

    out.fQ = kf_1.Charge() + kf_2.Charge();
    out.fNDF += 2;

#if R2DS_DEBUG
    Logger::Debug(__FUNCTION__, "mS     = {}", mS);
    Logger::Debug(__FUNCTION__, "zeta   = {}", zeta);
    Logger::Debug(__FUNCTION__, "mCHt   = {}", mCHt);
    Logger::Debug(__FUNCTION__, "mK     = {}", mK);
    Logger::Debug(__FUNCTION__, "D      = {}", D);
    Logger::Debug(__FUNCTION__, "K      = {}", K);
    Logger::Debug(__FUNCTION__, "K2     = {}", K2);
    Logger::Debug(__FUNCTION__, "A      = {}", A);
    Logger::Debug(__FUNCTION__, "M      = {}", M);
    Logger::Debug(__FUNCTION__, "out.fP = {}", out.fP);
    Logger::Debug(__FUNCTION__, "out.fC = {}", out.fC);
#endif

    return out;
}

Particle FitVertex(const KalmanFitter::Particle& part_1, const KalmanFitter::Particle& part_2, const Seeder::Result& s_1, const Seeder::Result& s_2,
                   double bz) {

    Daughter kf_1(part_1);
    Daughter kf_2(part_2);

    kf_1.PrepareJacobAndCorr(s_1, bz);
    kf_2.PrepareJacobAndCorr(s_2, bz);

    kf_1.Transport(s_1.seed.pca, kf_2.bt_cov);
    kf_2.Transport(s_2.seed.pca, kf_1.bt_cov);

    return GetUpdated(kf_1, kf_2);
}

}  // namespace R2DS::KalmanFitter
