#include <limits>

#include "Math/Point2D.h"
#include "Math/Vector2D.h"
#include "Math/Vector3D.h"

#include "Math/Constants.hxx"
#include "Math/Vertexer.hxx"

using XYPoint = ROOT::Math::XYPoint;
using XYZPoint = ROOT::Math::XYZPoint;
using XYVector = ROOT::Math::XYVector;
using XYZVector = ROOT::Math::XYZVector;

namespace Tree2Secondaries::Vertexer {

// Calculate the DCA between two charged particles (assuming both have helicoidal trajectory).
Particle::Pair MinimizeDistanceHelixHelix(const Charged& q, const Charged& t, const Helper::Propagator& prop) {

    // 1st approximation : DCA in XY //

    double omega1{q.Omega(prop)};
    double omega1_sq{omega1 * omega1};
    double omega2{t.Omega(prop)};
    double omega2_sq{omega2 * omega2};

    double px01{q.Px0()};
    double py01{q.Py0()};
    double pz01{q.Pz0()};

    double px02{t.Px0()};
    double py02{t.Py0()};
    double pz02{t.Pz0()};

    double xc1{q.X0() + py01 / omega1};
    double yc1{q.Y0() - px01 / omega1};
    double z01{q.Z0()};

    double xc2{t.X0() + py02 / omega2};
    double yc2{t.Y0() - px02 / omega2};
    double z02{t.Z0()};

    double rx{xc1 - xc2};
    double ry{yc1 - yc2};
    double r2{rx * rx + ry * ry};  // squared distance between circles' centers

    double pt1_sq{q.Pt() * q.Pt()};
    double pt2_sq{t.Pt() * t.Pt()};
    double Pxy1dotDR{px01 * rx + py01 * ry};
    double Pxy2dotDR{px02 * rx + py02 * ry};
    double Pxy1crossDR{px01 * ry - py01 * rx};
    double Pxy2crossDR{px02 * ry - py02 * rx};

    double radius1{q.Pt() / std::abs(omega1)};
    double radius2{t.Pt() / std::abs(omega2)};
    // double dist_btw_centers{std::sqrt(r2)};
    double radius_sum_sq{(radius1 + radius2) * (radius1 + radius2)};
    double radius_diff_sq{(radius1 - radius2) * (radius1 - radius2)};

    double theta1{0.};
    double sin_theta1{0.};
    double cos_theta1{0.};
    ROOT::Math::XYZPoint pca1{0., 0., 0.};

    double theta2{0.};
    double sin_theta2{0.};
    double cos_theta2{0.};
    ROOT::Math::XYZPoint pca2{0., 0., 0.};

    double dca_3d{std::numeric_limits<double>::max()};

    // opposite of intersection condition, which is: "D < R1+R2" && "D > |R1-R2|"
    if (r2 > radius_sum_sq || r2 < radius_diff_sq) {
        // then, PCA are colinear to both circle centers //
        // 4 possible solutions, choose the one that minimizes distance in 3D //
        for (auto sign1 : {+1, -1}) {
            double tmp_theta1{std::atan2(sign1 * Pxy1dotDR, sign1 * Pxy1crossDR)};
            double tmp_sin_theta1{0.};
            double tmp_cos_theta1{0.};
            sincos(tmp_theta1, &tmp_sin_theta1, &tmp_cos_theta1);
            ROOT::Math::XYZPoint tmp_pca1{{xc1 + (px01 * tmp_sin_theta1 - py01 * tmp_cos_theta1) / omega1},  //
                                          {yc1 + (px01 * tmp_cos_theta1 + py01 * tmp_sin_theta1) / omega1},
                                          {z01 + pz01 * tmp_theta1 / omega1}};
            //
            for (auto sign2 : {+1, -1}) {
                double tmp_theta2{std::atan2(sign2 * Pxy2dotDR, sign2 * Pxy2crossDR)};
                double tmp_sin_theta2{0.};
                double tmp_cos_theta2{0.};
                sincos(tmp_theta2, &tmp_sin_theta2, &tmp_cos_theta2);
                ROOT::Math::XYZPoint tmp_pca2{{xc2 + (px02 * tmp_sin_theta2 - py02 * tmp_cos_theta2) / omega2},  //
                                              {yc2 + (px02 * tmp_cos_theta2 + py02 * tmp_sin_theta2) / omega2},
                                              {z02 + pz02 * tmp_theta2 / omega2}};
                //
                auto tmp_dca{Math::DistanceSquared(tmp_pca1, tmp_pca2)};
                if (tmp_dca < dca_3d) {
                    dca_3d = tmp_dca;

                    theta1 = tmp_theta1;
                    sin_theta1 = tmp_sin_theta1;
                    cos_theta1 = tmp_cos_theta1;
                    pca1 = tmp_pca1;

                    theta2 = tmp_theta2;
                    sin_theta2 = tmp_sin_theta2;
                    cos_theta2 = tmp_cos_theta2;
                    pca2 = tmp_pca2;
                }
            }
        }
    } else {
        // circles intersect : "D < R1+R2" && "D > |R1-R2|"
        double cc{omega1_sq * omega2_sq};
        double kk{omega1_sq * pt2_sq - omega2_sq * pt1_sq};
        double k1{-kk + cc * r2};
        double k2{kk + cc * r2};
        double k3{2 * omega2_sq * pt1_sq + 2 * omega1_sq * pt2_sq - cc * r2};
        double dk{cc * r2 * k3 - kk * kk};
        double sqrt_dk{dk > 0 ? std::sqrt(dk) : 0.};
        double d1{std::abs(omega1) * std::abs(Pxy1dotDR) * sqrt_dk};
        double d2{std::abs(omega2) * std::abs(Pxy2dotDR) * sqrt_dk};
        if (std::abs(Pxy1dotDR) < Const::AlmostZero) Pxy1dotDR = Const::LocalSmall;
        if (std::abs(Pxy2dotDR) < Const::AlmostZero) Pxy2dotDR = Const::LocalSmall;
        // 2 possible solutions, choose the one that minimizes distance in 3D //
        for (auto sign : {+1, -1}) {
            //
            double tmp_theta1{std::atan2(-omega1 * k1 * Pxy1dotDR + sign * d1 * Pxy1crossDR / Pxy1dotDR, -omega1 * k1 * Pxy1crossDR - d1)};
            double tmp_sin_theta1{0.};
            double tmp_cos_theta1{0.};
            sincos(tmp_theta1, &tmp_sin_theta1, &tmp_cos_theta1);
            ROOT::Math::XYZPoint tmp_pca1{{xc1 + (px01 * tmp_sin_theta1 - py01 * tmp_cos_theta1) / omega1},  //
                                          {yc1 + (px01 * tmp_cos_theta1 + py01 * tmp_sin_theta1) / omega1},
                                          {z01 + pz01 * tmp_theta1 / omega1}};
            //
            double tmp_theta2{std::atan2(omega2 * k2 * Pxy2dotDR + sign * d2 * Pxy2crossDR / Pxy2dotDR, omega2 * k2 * Pxy2crossDR - d2)};
            double tmp_sin_theta2{0.};
            double tmp_cos_theta2{0.};
            sincos(tmp_theta2, &tmp_sin_theta2, &tmp_cos_theta2);
            ROOT::Math::XYZPoint tmp_pca2{{xc2 + (px02 * tmp_sin_theta2 - py02 * tmp_cos_theta2) / omega2},  //
                                          {yc2 + (px02 * tmp_cos_theta2 + py02 * tmp_sin_theta2) / omega2},
                                          {z02 + pz02 * tmp_theta2 / omega2}};
            //
            double tmp_dca{Math::DistanceSquared(tmp_pca1, tmp_pca2)};
            if (tmp_dca < dca_3d) {
                dca_3d = tmp_dca;

                theta1 = tmp_theta1;
                sin_theta1 = tmp_sin_theta1;
                cos_theta1 = tmp_cos_theta1;
                pca1 = tmp_pca1;

                theta2 = tmp_theta2;
                sin_theta2 = tmp_sin_theta2;
                cos_theta2 = tmp_cos_theta2;
                pca2 = tmp_pca2;
            }
        }
    }

    // 2nd approximation: add Z-component //

    double px1{px01 * cos_theta1 + py01 * sin_theta1};
    double py1{-px01 * sin_theta1 + py01 * cos_theta1};

    double px2{px02 * cos_theta2 + py02 * sin_theta2};
    double py2{-px02 * sin_theta2 + py02 * cos_theta2};

    double Pxy1dotPxy2{px1 * px2 + py1 * py2};
    double P1dotP2{Pxy1dotPxy2 + pz01 * pz02};
    Pxy1crossDR = px1 * ry - py1 * rx;
    Pxy2crossDR = px2 * ry - py2 * rx;

    double bb{pz01 * pz01 - omega1 * Pxy1crossDR + omega1 / omega2 * Pxy1dotPxy2};
    double dd{-pz02 * pz02 - omega2 * Pxy2crossDR - omega2 / omega1 * Pxy1dotPxy2};
    double den{bb * dd + P1dotP2 * P1dotP2};

    double ene1{q.Energy()};
    double ene2{t.Energy()};

    // protection against zero
    if (std::abs(den) < Const::AlmostZero) {
        return {
            {{px1, py1, pz01, ene1}, pca1},
            {{px2, py2, pz02, ene2}, pca2},
        };
    }

    double rz{z01 - z02};

    double Pxy1crossPxy2{px1 * py2 - py1 * px2};
    double P1dotDR{px1 * rx + py1 * ry + pz01 * rz};
    double P2dotDR{px2 * rx + py2 * ry + pz02 * rz};

    double sr1{theta1 / omega1};
    double sr2{theta2 / omega2};

    double aa{P1dotDR + pz01 * pz01 * sr1 - pz01 * pz02 * sr2 + Pxy1crossPxy2 / omega2};
    double cc{P2dotDR - pz02 * pz02 * sr2 + pz01 * pz02 * sr1 + Pxy1crossPxy2 / omega1};

    double sz1{(-aa * dd - cc * P1dotP2) / den};
    double sz2{(-bb * cc + aa * P1dotP2) / den};

    // update final properties //

    double fnl_theta1{omega1 * (sr1 + sz1)};
    double fnl_sin_theta1{std::sin(fnl_theta1)};
    double fnl_cos_theta1{std::cos(fnl_theta1)};

    double fnl_theta2{omega2 * (sr2 + sz2)};
    double fnl_sin_theta2{std::sin(fnl_theta2)};
    double fnl_cos_theta2{std::cos(fnl_theta2)};

    return {{{px01 * fnl_cos_theta1 + py01 * fnl_sin_theta1,   //
              -px01 * fnl_sin_theta1 + py01 * fnl_cos_theta1,  //
              pz01, ene1},
             {xc1 + (px01 * fnl_sin_theta1 - py01 * fnl_cos_theta1) / omega1,  //
              yc1 + (px01 * fnl_cos_theta1 + py01 * fnl_sin_theta1) / omega1,  //
              z01 + pz01 * (sr1 + sz1)}},
            {{px02 * fnl_cos_theta2 + py02 * fnl_sin_theta2,   //
              -px02 * fnl_sin_theta2 + py02 * fnl_cos_theta2,  //
              pz02, ene2},
             {xc2 + (px02 * fnl_sin_theta2 - py02 * fnl_cos_theta2) / omega2,  //
              yc2 + (px02 * fnl_cos_theta2 + py02 * fnl_sin_theta2) / omega2,  //
              z02 + pz02 * (sr2 + sz2)}}};
}

// Calculate the DCA between a charged particle (assuming a helicoidal trajectory) and a point in space.
Particle::State MinimizeDistanceHelixVertex(const Charged& q, const ROOT::Math::XYZPoint& v, const Helper::Propagator& prop) {

    // 1st approximation : DCA in XY                //
    // helix projection on xy-plane : circumference //

    double omega{q.Omega(prop)};
    double xc{q.X0() + q.Py0() / omega};
    double yc{q.Y0() - q.Px0() / omega};
    double rx{xc - v.X()};
    double ry{yc - v.Y()};
    double PxydotDR{q.Px0() * rx + q.Py0() * ry};
    double PxycrossDR{q.Px0() * ry - q.Py0() * rx};

    double dca_3d{std::numeric_limits<double>::max()};
    double sr{0.};

    // PCA is colinear to circle center and vertex //
    // 2 possible solutions, choose the one that minimizes distance in 3D //
    for (auto sign : {+1, -1}) {
        double tmp_sr{std::atan2(sign * PxydotDR, sign * PxycrossDR) / omega};
        double tmp_dca{(q.XYZ(tmp_sr, prop) - v).Mag2()};
        if (tmp_dca < dca_3d) {
            dca_3d = tmp_dca;
            sr = tmp_sr;
        }
    }

    // 2nd approximation : add z-component //

    double rz{q.Z0() - v.Z()};

    double px{q.Px(sr, prop)};
    double py{q.Py(sr, prop)};
    PxycrossDR = px * ry - py * rx;
    double P1dotDR = px * rx + py * ry + q.Pz0() * rz;

    double sz{0.};
    double det{q.Pz0() * q.Pz0() - omega * PxycrossDR};
    if (std::abs(det) > Const::AlmostZero) {
        sz = (-q.Pz0() * q.Pz0() * sr - P1dotDR) / det;
    }

    return q.PropagatedState(sr + sz, prop);
}

// Calculate the DCA between a charged particle (assuming a helicoidal trajectory) and a neutral particule (assuming a straight line trajectory).
Particle::Pair MinimizeDistanceHelixLine(const Charged& q, const Neutral& n, const Helper::Propagator& prop) {

    // 1st approximation : DCA in XY //

    double omega1{q.Omega(prop)};
    double omega1_sq{omega1 * omega1};

    double rx{q.X0() + q.Py0() / omega1 - n.DecayX()};
    double ry{q.Y0() - q.Px0() / omega1 - n.DecayY()};

    double pt1sq{q.Pt() * q.Pt()};
    double pt2sq{n.Pt() * n.Pt()};

    double Pxy1dotPxy2{q.Px0() * n.Px() + q.Py0() * n.Py()};
    double Pxy1crossPxy2{q.Px0() * n.Py() - q.Py0() * n.Px()};
    double Pxy2dotDR{n.Px() * rx + n.Py() * ry};
    double Pxy2crossDR{n.Px() * ry - n.Py() * rx};

    double dca_3d{std::numeric_limits<double>::max()};
    double sr1{0.}, sr2{0.};

    double discrim{omega1_sq * (pt1sq * pt2sq - omega1_sq * Pxy2crossDR * Pxy2crossDR)};
    if (discrim < 0.) {
        // no intersection in the XY plane //
        sr2 = Pxy2dotDR / pt2sq;
        // 2 possible solutions, choose the one that minimizes distance in 3D //
        for (auto sign : {+1, -1}) {
            double tmp_sr1{std::atan2(sign * Pxy1crossPxy2, -sign * Pxy1dotPxy2) / omega1};
            double tmp_dca{(q.XYZ(tmp_sr1, prop) - n.XYZ(sr2, prop)).Mag2()};
            if (tmp_dca < dca_3d) {
                dca_3d = tmp_dca;
                sr1 = tmp_sr1;
            }
        }
    } else {
        // circle and line intersect //
        // 2 possible solutions, choose the one that minimizes distance in 3D //
        for (auto sign : {+1, -1}) {
            discrim = std::sqrt(discrim);
            double tmp_sr1{std::atan2(omega1 * Pxy2crossDR * Pxy1crossPxy2 - sign * discrim * Pxy1dotPxy2 / omega1,   //
                                      -omega1 * Pxy2crossDR * Pxy1dotPxy2 - sign * discrim * Pxy1crossPxy2 / omega1)  //
                           / omega1};
            double tmp_sr2{(Pxy2dotDR * omega1_sq - sign * discrim) / (omega1_sq * pt2sq)};
            double tmp_dca{(q.XYZ(tmp_sr1, prop) - n.XYZ(tmp_sr2, prop)).Mag2()};
            if (tmp_dca < dca_3d) {
                dca_3d = tmp_dca;
                sr1 = tmp_sr1;
                sr2 = tmp_sr2;
            }
        }
    }

    // 2nd approximation: add Z-component //

    double rz{q.Z0() - n.DecayZ()};

    double px1{q.Px(sr1, prop)};
    double py1{q.Py(sr1, prop)};

    Pxy1dotPxy2 = px1 * n.Px() + py1 * n.Py();
    Pxy1crossPxy2 = px1 * n.Py() - py1 * n.Px();
    double P1dotP2{Pxy1dotPxy2 + q.Pz0() * n.Pz()};
    double P1dotDR{px1 * rx + py1 * ry + q.Pz0() * rz};
    double P2dotDR{n.Px() * rx + n.Py() * ry + n.Pz() * rz};
    double Pxy1crossDR{px1 * ry - py1 * rx};

    double aa{q.Pz0() * q.Pz0() * sr1 + P1dotDR - sr2 * P1dotP2};
    double bb{q.Pz0() * q.Pz0() - omega1 * Pxy1crossDR + omega1 * sr2 * Pxy1crossPxy2};
    double cc{P2dotDR + q.Pz0() * n.Pz() * sr1 - n.P2() * sr2 + Pxy1crossPxy2 / omega1};

    double sz1{0.};
    double sz2{0.};
    double den{P1dotP2 * P1dotP2 - bb * n.P2()};
    if (std::abs(den) > Const::AlmostZero) {
        sz1 = (aa * n.P2() - cc * P1dotP2) / den;
        sz2 = (-bb * cc + aa * P1dotP2) / den;
    }

    return {q.PropagatedState(sr1 + sz1, prop), n.PropagatedState(sr2 + sz2, prop)};
}

// Calculate the DCA between two neutral particles (assuming both have a straight line trajectory).
Particle::Pair MinimizeDistanceLineLine(const Neutral& n, const Neutral& m, const Helper::Propagator& prop) {

    double p1p1{n.P2()};
    double p1p2{n.PxPyPz().Dot(m.PxPyPz())};
    double p2p2{m.P2()};

    double det{p1p2 * p1p2 - p1p1 * p2p2};
    if (std::abs(det) < Const::AlmostZero) {
        // rare case : particles are parallel, minimize both to a same arbitrary vertex (origin) //
        return {MinimizeDistanceLineVertex(n, {0., 0., 0.}, prop),  //
                MinimizeDistanceLineVertex(m, {0., 0., 0.}, prop)};
    }

    XYZVector dr{n.DecayVertex() - m.DecayVertex()};
    double p1dr{n.PxPyPz().Dot(dr)};
    double p2dr{m.PxPyPz().Dot(dr)};

    return {n.PropagatedState((p2p2 * p1dr - p1p2 * p2dr) / det, prop), m.PropagatedState((p1p1 * p2dr - p1p2 * p1dr) / det, prop)};
}

// Calculate the DCA between a neutral particle (assuming a straight line trajectory) and a point in space.
Particle::State MinimizeDistanceLineVertex(const Neutral& n, const ROOT::Math::XYZPoint& v, const Helper::Propagator& prop) {
    return n.PropagatedState(-n.PxPyPz().Dot(n.DecayVertex() - v) / n.P2(), prop);
}

}  // namespace Tree2Secondaries::Vertexer
