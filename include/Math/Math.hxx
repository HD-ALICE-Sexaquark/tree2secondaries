#ifndef T2S_MATH_COMMON_HXX
#define T2S_MATH_COMMON_HXX

#include <tuple>

#include <Math/Point3D.h>
#include <Math/Vector3D.h>
#include <Math/VectorUtil.h>

#include <KFParticle.hxx>
#include <KFParticle_Math.hxx>

#include "Math/Constants.hxx"

namespace Tree2Secondaries::Math {

// Calculate the cosine of the pointing angle of a particle with momentum Px,Py,Pz and vertex X,Y,Z w.r.t. to a reference point.
inline double CosinePointingAngle(const ROOT::Math::XYZVector& mom_v0, const ROOT::Math::XYZPoint& pos_v0, const ROOT::Math::XYZPoint& pos_ref) {
    return std::cos(ROOT::Math::VectorUtil::Angle(mom_v0, pos_v0 - pos_ref));
}

// Overload of `CosinePointingAngle(...)`, using scalars instead of vectors.
inline double CosinePointingAngle(double v0_px, double v0_py, double v0_pz,  //
                                  double v0_x, double v0_y, double v0_z,     //
                                  double ref_x, double ref_y, double ref_z) {
    return CosinePointingAngle({v0_px, v0_py, v0_pz}, {v0_x, v0_y, v0_z}, {ref_x, ref_y, ref_z});
}

// Calculate Armenteros-Podolanski qT.
// Based on https://github.com/alisw/AliRoot (`STEER/ESD/AliESDv0::PtArmV0()`)
inline double ArmenterosQt(const ROOT::Math::XYZVector& mom_v0, const ROOT::Math::XYZVector& mom_dau) {
    return ROOT::Math::VectorUtil::Perp(mom_v0, mom_dau);
}

// Overload of `ArmenterosQt(...)`
inline double ArmenterosQt(double v0_px, double v0_py, double v0_pz,  //
                           double dau_px, double dau_py, double dau_pz) {
    return ArmenterosQt({v0_px, v0_py, v0_pz}, {dau_px, dau_py, dau_pz});
}

// Calculate Armenteros-Podolanski alpha.
// Based on https://github.com/alisw/AliRoot (`STEER/ESD/AliESDv0::AlphaV0()`)
inline double ArmenterosAlpha(const ROOT::Math::XYZVector& mom_v0, const ROOT::Math::XYZVector& mom_neg, const ROOT::Math::XYZVector& mom_pos) {
    double lQlNeg = mom_neg.Dot(mom_v0) / mom_v0.R();
    double lQlPos = mom_pos.Dot(mom_v0) / mom_v0.R();
    if (std::abs(lQlPos + lQlNeg) < 1E-6) return 2.;  // protection
    return (lQlPos - lQlNeg) / (lQlPos + lQlNeg);
}

// Overload of `ArmenterosAlpha()`.
inline double ArmenterosAlpha(double v0_px, double v0_py, double v0_pz,     //
                              double neg_px, double neg_py, double neg_pz,  //
                              double pos_px, double pos_py, double pos_pz) {
    return ArmenterosAlpha({v0_px, v0_py, v0_pz}, {neg_px, neg_py, neg_pz}, {pos_px, pos_py, pos_pz});
}

// Return distance of closest approach (DCA) between a neutral particle (line trajectory) and a point in space.
inline double FastDCALineVertex(const ROOT::Math::XYZVector& mom_v0, const ROOT::Math::XYZPoint& pos_v0, const ROOT::Math::XYZPoint& pos_ref) {
    return std::abs(mom_v0.Dot(pos_v0 - pos_ref) / mom_v0.R());
}

// Calculate the middle-point between two points in 3D space.
inline ROOT::Math::XYZPoint MiddlePoint(const ROOT::Math::XYZPoint& p1, const ROOT::Math::XYZPoint& p2) {
    return {(p1.X() + p2.X()) / 2., (p1.Y() + p2.Y()) / 2., (p1.Z() + p2.Z()) / 2.};
}

// Calculate square of the distance between two points in 3D space.
inline double DistanceSquared(const ROOT::Math::XYZPoint& p1, const ROOT::Math::XYZPoint& p2) { return (p1 - p2).Mag2(); }

// Return distance of closest approach (DCA) between a charged particle (helix trajectory under an homogeneous magnetic field) and a point in space.
inline double FastDCAHelixVertex(const KF::Particle& p, double ref_x, double ref_y, double ref_z, float bz) {

    ROOT::Math::XYZPoint ref{ref_x, ref_y, ref_z};

    double bq{bz * static_cast<double>(p.Charge()) * KF::Const::Kappa};

    double px0{p.Px()};
    double py0{p.Py()};
    double pz0{p.Pz()};
    double pt2{px0 * px0 + py0 * py0};
    double x0{p.X()};
    double y0{p.Y()};
    double z0{p.Z()};

    double dx{ref_x - x0};
    double dy{ref_y - y0};
    double dz{ref_z - z0};
    double a{dx * px0 + dy * py0};

    double abq{bq * a};
    double bbq{bq * (dx * py0 - dy * px0) - pt2};

    // 1.a -- get solution and update cache properties //

    double theta{std::atan2(abq, -bbq)};
    auto [sin, cos] = KF::Math::sincos(theta);
    double sB{sin / bq};
    double cB{(1. - cos) / bq};

    double ds{theta / bq};

    ROOT::Math::XYZPoint pca{x0 + sB * px0 + cB * py0,  //
                             y0 - cB * px0 + sB * py0,  //
                             z0 + ds * pz0};

    // 2 -- add z-component as small correction //

    double cbq{bbq * cos - abq * sin - pz0 * pz0};
    if (std::abs(cbq) < KF::Const::AbsAlmostZero) return (ref - pca).R();  // protection

    double sz{(ds * pz0 - dz) * pz0 / cbq};

    // 2.b -- update ds //

    ds += sz;

    // 2.c -- update rest of cache properties //

    theta = bq * ds;
    std::tie(sin, cos) = KF::Math::sincos(theta);
    sB = sin / bq;
    cB = (1. - cos) / bq;

    pca.SetXYZ(x0 + sB * px0 + cB * py0,  //
               y0 - cB * px0 + sB * py0,  //
               z0 + ds * pz0);

    return (ref - pca).R();
}

// Return distance of closest approach (DCA) between two charged particles (helix trajectories) under an homogeneus magnetic field.
inline double FastDCAHelixHelix(const KF::Particle& p, const KF::Particle& q, float bz) {

    double dca_sq{Const::BigNumber};

    // 1 -- find points of closest approach (PCAs) in XY plane //

    double bq1{bz * static_cast<double>(p.Charge()) * KF::Const::Kappa};
    double bq2{bz * static_cast<double>(q.Charge()) * KF::Const::Kappa};

    double px01{p.Px()};
    double py01{p.Py()};
    double pz01{p.Pz()};
    double px02{q.Px()};
    double py02{q.Py()};
    double pz02{q.Pz()};
    double pt12{px01 * px01 + py01 * py01};
    double pt22{px02 * px02 + py02 * py02};
    double x01{p.X()};
    double y01{p.Y()};
    double z01{p.Z()};
    double x02{q.X()};
    double y02{q.Y()};
    double z02{q.Z()};

    double dx0{x01 - x02};
    double dy0{y01 - y02};
    double dr02{dx0 * dx0 + dy0 * dy0};
    double drp1{dx0 * px01 + dy0 * py01};
    double dxyp1{dx0 * py01 - dy0 * px01};
    double drp2{dx0 * px02 + dy0 * py02};
    double dxyp2{dx0 * py02 - dy0 * px02};
    double p1p2{px01 * px02 + py01 * py02};
    double dp1p2{px01 * py02 - px02 * py01};

    double k11{bq2 * drp1 - dp1p2};
    double k21{bq1 * (bq2 * dxyp1 - p1p2) + bq2 * pt12};
    double k12{bq1 * drp2 - dp1p2};
    double k22{bq2 * (bq1 * dxyp2 + p1p2) - bq1 * pt22};

    double kp{dxyp1 * bq2 - dxyp2 * bq1 - p1p2};
    double kd{dr02 * bq1 * bq2 / 2. + kp};
    double c1{-bq1 * kd - pt12 * bq2};
    double c2{bq2 * kd + pt22 * bq1};

    double d1{std::sqrt(std::max(pt12 * pt22 - kd * kd, 0.))};

    // 1.a -- select solution with minimum distance in 3D  //

    double ds1{0.};
    double ds2{0.};

    ROOT::Math::XYZPoint pca1{};
    ROOT::Math::XYZPoint pca2{};

    double px1{px01};
    double py1{py01};
    double px2{px02};
    double py2{py02};

    for (auto sign : {+1, -1}) {
        // particle 1 //
        double tmp_theta1{std::atan2(bq1 * (k11 * c1 + sign * k21 * d1), sign * bq1 * k11 * d1 * bq1 - k21 * c1)};
        auto [tmp_sin1, tmp_cos1] = KF::Math::sincos(tmp_theta1);
        double tmp_sB1{tmp_sin1 / bq1};
        double tmp_cB1{(1. - tmp_cos1) / bq1};
        double tmp_ds1{tmp_theta1 / bq1};

        ROOT::Math::XYZPoint tmp_pca1{x01 + tmp_sB1 * px01 + tmp_cB1 * py01,  //
                                      y01 - tmp_cB1 * px01 + tmp_sB1 * py01,  //
                                      z01 + tmp_ds1 * pz01};

        // particle 2 //
        double tmp_theta2{std::atan2(bq2 * (k12 * c2 + sign * k22 * d1), sign * bq2 * k12 * d1 * bq2 - k22 * c2)};
        auto [tmp_sin2, tmp_cos2] = KF::Math::sincos(tmp_theta2);
        double tmp_sB2{tmp_sin2 / bq2};
        double tmp_cB2{(1. - tmp_cos2) / bq2};
        double tmp_ds2{tmp_theta2 / bq2};

        ROOT::Math::XYZPoint tmp_pca2{x02 + tmp_sB2 * px02 + tmp_cB2 * py02,  //
                                      y02 - tmp_cB2 * px02 + tmp_sB2 * py02,  //
                                      z02 + tmp_ds2 * pz02};

        // store //
        double tmp_dca_sq{(tmp_pca2 - tmp_pca1).Mag2()};
        if (tmp_dca_sq < dca_sq) {
            ds1 = tmp_ds1;
            ds2 = tmp_ds2;

            pca1 = tmp_pca1;
            pca2 = tmp_pca2;

            dca_sq = tmp_dca_sq;

            px1 = tmp_cos1 * px01 + tmp_sin1 * py01;
            py1 = -tmp_sin1 * px01 + tmp_cos1 * py01;
            px2 = tmp_cos2 * px02 + tmp_sin2 * py02;
            py2 = -tmp_sin2 * px02 + tmp_cos2 * py02;
        }
    }

    // 2 -- add z-component as small correction //

    double p12{px1 * px1 + py1 * py1 + pz01 * pz01};
    double p22{px2 * px2 + py2 * py2 + pz02 * pz02};
    double lp1p2{px1 * px2 + py1 * py2 + pz01 * pz02};

    double detp{lp1p2 * lp1p2 - p12 * p22};  // protection
    if (std::abs(detp) < KF::Const::AbsAlmostZero || detp * detp < KF::Const::AbsAlmostZero) {
        return std::sqrt(dca_sq);
    }

    // 2.b -- update ds //

    double dx{pca2.X() - pca1.X()};
    double dy{pca2.Y() - pca1.Y()};
    double dz{pca2.Z() - pca1.Z()};

    double ldrp1{px1 * dx + py1 * dy + pz01 * dz};
    double ldrp2{px2 * dx + py2 * dy + pz02 * dz};
    double a1{ldrp2 * lp1p2 - ldrp1 * p22};
    double a2{ldrp2 * p12 - ldrp1 * lp1p2};

    ds1 += a1 / detp;
    ds2 += a2 / detp;

    // 2.c -- get DCA //

    double theta1{bq1 * ds1};
    auto [sin1, cos1] = KF::Math::sincos(theta1);
    double sB1{sin1 / bq1};
    double cB1{(1. - cos1) / bq1};

    pca1.SetXYZ(x01 + sB1 * px01 + cB1 * py01,  //
                y01 - cB1 * px01 + sB1 * py01,  //
                z01 + ds1 * pz01);

    double theta2{bq2 * ds2};
    auto [sin2, cos2] = KF::Math::sincos(theta2);
    double sB2{sin2 / bq2};
    double cB2{(1. - cos2) / bq2};

    pca2.SetXYZ(x02 + sB2 * px02 + cB2 * py02,  //
                y02 - cB2 * px02 + sB2 * py02,  //
                z02 + ds2 * pz02);

    dca_sq = (pca2 - pca1).Mag2();
    return std::sqrt(dca_sq);
}

}  // namespace Tree2Secondaries::Math

#endif  // T2S_MATH_COMMON_HXX
