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

    double omega1 = q.Omega(prop);
    double omega2 = t.Omega(prop);

    double px01{q.Px0()};
    double py01{q.Py0()};
    double pz01{q.Pz0()};

    double px02{t.Px0()};
    double py02{t.Py0()};
    double pz02{t.Pz0()};

    double xc1{q.X0() + py01 / omega1};
    double yc1{q.Y0() - px01 / omega1};
    double xc2{t.X0() + py02 / omega2};
    double yc2{t.Y0() - px02 / omega2};
    double rx{xc1 - xc2};
    double ry{yc1 - yc2};

    double omega1sq{omega1 * omega1};
    double omega2sq{omega2 * omega2};
    double r2{rx * rx + ry * ry};
    double pt1sq{q.Pt0() * q.Pt0()};
    double pt2sq{t.Pt0() * t.Pt0()};
    double Pxy1dotDR{px01 * rx + py01 * ry};
    double Pxy2dotDR{px02 * rx + py02 * ry};
    double Pxy1crossDR{px01 * ry - py01 * rx};
    double Pxy2crossDR{px02 * ry - py02 * rx};

    double radius1{q.Pt0() / std::abs(omega1)};
    double radius2{t.Pt0() / std::abs(omega2)};
    double dist_btw_centers = std::sqrt(r2);
    bool intersection = dist_btw_centers < radius1 + radius2 && dist_btw_centers > std::abs(radius1 - radius2);

    double dca_3d{std::numeric_limits<double>::max()};
    double sr1{0.}, sr2{0.};

    if (!intersection) {
        // PCA are colinear to both circle centers //
        // 4 possible solutions, choose the one that minimizes distance in 3D //
        for (auto sign1 : {+1, -1}) {
            for (auto sign2 : {+1, -1}) {
                double tmp_sr1{std::atan2(sign1 * Pxy1dotDR, sign1 * Pxy1crossDR) / omega1};
                double tmp_sr2{std::atan2(sign2 * Pxy2dotDR, sign2 * Pxy2crossDR) / omega2};
                double tmp_dca{(q.XYZ(tmp_sr1, prop) - t.XYZ(tmp_sr2, prop)).Mag2()};
                if (tmp_dca < dca_3d) {
                    dca_3d = tmp_dca;
                    sr1 = tmp_sr1;
                    sr2 = tmp_sr2;
                }
            }
        }
    } else {
        // circles intersect //
        double kk{omega1sq * pt2sq - omega2sq * pt1sq};
        double k1{-kk + omega1sq * omega2sq * r2};
        double k2{kk + omega1sq * omega2sq * r2};
        double k3{2 * omega2sq * pt1sq + 2 * omega1sq * pt2sq - omega1sq * omega2sq * r2};
        double dd{omega1sq * omega2sq * r2 * k3 - kk * kk};
        double d1{dd > 0 ? std::abs(omega1) * std::abs(Pxy1dotDR) * std::sqrt(dd) : 0.};
        double d2{dd > 0 ? std::abs(omega2) * std::abs(Pxy2dotDR) * std::sqrt(dd) : 0.};
        if (std::abs(Pxy1dotDR) < Const::AlmostZero) Pxy1dotDR = Const::LocalSmall;
        if (std::abs(Pxy2dotDR) < Const::AlmostZero) Pxy2dotDR = Const::LocalSmall;
        // 2 possible solutions, choose the one that minimizes distance in 3D //
        for (auto sign : {+1, -1}) {
            double tmp_sr1{std::atan2(-omega1 * k1 * Pxy1dotDR + sign * d1 * Pxy1crossDR / Pxy1dotDR,  //
                                      -omega1 * k1 * Pxy1crossDR - d1) /
                           omega1};
            double tmp_sr2{std::atan2(omega2 * k2 * Pxy2dotDR + sign * d2 * Pxy2crossDR / Pxy2dotDR,  //
                                      omega2 * k2 * Pxy2crossDR - d2) /
                           omega2};
            double tmp_dca{(q.XYZ(tmp_sr1, prop) - t.XYZ(tmp_sr2, prop)).Mag2()};
            if (tmp_dca < dca_3d) {
                dca_3d = tmp_dca;
                sr1 = tmp_sr1;
                sr2 = tmp_sr2;
            }
        }
    }

    // 2nd approximation: add Z-component //

    double rz{q.Z0() - t.Z0()};

    double px1{q.Px(sr1, prop)};
    double py1{q.Py(sr1, prop)};

    double px2{t.Px(sr2, prop)};
    double py2{t.Py(sr2, prop)};

    double Pxy1crossPxy2{px1 * py2 - py1 * px2};
    double Pxy1dotPxy2{px1 * px2 + py1 * py2};
    double P1dotDR{px1 * rx + py1 * ry + pz01 * rz};
    double P2dotDR{px2 * rx + py2 * ry + pz02 * rz};
    double P1dotP2{Pxy1dotPxy2 + pz01 * pz02};
    Pxy1crossDR = px1 * ry - py1 * rx;
    Pxy2crossDR = px2 * ry - py2 * rx;

    double aa{P1dotDR + pz01 * pz01 * sr1 - pz01 * pz02 * sr2 + Pxy1crossPxy2 / omega2};
    double bb{pz01 * pz01 - omega1 * Pxy1crossDR + omega1 * Pxy1dotPxy2 / omega2};
    double cc{P2dotDR - pz02 * pz02 * sr2 + pz01 * pz02 * sr1 + Pxy1crossPxy2 / omega1};
    double dd{-pz02 * pz02 - omega2 * Pxy2crossDR - omega2 * Pxy1dotPxy2 / omega1};
    double den{bb * dd + P1dotP2 * P1dotP2};

    double sz1{0.};
    double sz2{0.};
    if (std::abs(den) > Const::AlmostZero) {
        sz1 = (-aa * dd - cc * P1dotP2) / den;
        sz2 = (-bb * cc + aa * P1dotP2) / den;
    }

    return {q.PropagatedState(sr1 + sz1, prop), t.PropagatedState(sr2 + sz2, prop)};
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
    double omega1sq{omega1 * omega1};

    double rx{q.X0() + q.Py0() / omega1 - n.DecayX()};
    double ry{q.Y0() - q.Px0() / omega1 - n.DecayY()};

    double pt1sq{q.Pt0() * q.Pt0()};
    double pt2sq{n.Pt() * n.Pt()};

    double Pxy1dotPxy2{q.Px0() * n.Px() + q.Py0() * n.Py()};
    double Pxy1crossPxy2{q.Px0() * n.Py() - q.Py0() * n.Px()};
    double Pxy2dotDR{n.Px() * rx + n.Py() * ry};
    double Pxy2crossDR{n.Px() * ry - n.Py() * rx};

    double dca_3d{std::numeric_limits<double>::max()};
    double sr1{0.}, sr2{0.};

    double discrim{omega1sq * (pt1sq * pt2sq - omega1sq * Pxy2crossDR * Pxy2crossDR)};
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
            double tmp_sr2{(Pxy2dotDR * omega1sq - sign * discrim) / (omega1sq * pt2sq)};
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
