#ifndef ALICE_VERTEXER_HXX
#define ALICE_VERTEXER_HXX

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>

#include "ALICE/ESD.hxx"
#include "ALICE/Math.hxx"

// Based on `AliRoot/STEER/ESD/AliV0vertexer`
namespace ALICE::Vertexer {

// Rewrite of `AliExternalTrackParam::Evaluate()`.
// Calculate position of a point on a track and some derivatives
// Input arguments:
// - `h`
// - `t`
// Output arguments:
// - `r` : position?
// - `g` : first derivatives
// - `gg` : second derivatives
static void Evaluate(const std::array<double, 8> &h, double t, std::array<double, 3> &r, std::array<double, 3> &g, std::array<double, 3> &gg) {
    double phase{h[4] * t + h[2]};
    auto [sn, cs] = Math::sincos(phase);

    r[0] = h[5];
    r[1] = h[0];
    if (std::abs(h[4]) > Const::AlmostZero) {
        r[0] += (sn - h[6]) / h[4];
        r[1] -= (cs - h[7]) / h[4];
    } else {
        r[0] += t * cs;
        r[1] -= -t * sn;
    }
    r[2] = h[1] + h[3] * t;

    g[0] = cs;
    g[1] = sn;
    g[2] = h[3];

    gg[0] = -h[4] * sn;
    gg[1] = h[4] * cs;
    gg[2] = 0.;
}

// Rewrite of `Preoptimize()`
// This function pre-optimizes a two-track pair in the XY plane
// and provides two X values for the tracks if successful
static bool Preoptimize(const ALICE::Track &nt, const ALICE::Track &pt, double &lPreprocessxn, double &lPreprocessxp, double bz) {

    double lMinimumX{-3.0};
    double lMaximumX{300.0};

    std::array<double, 6> nhelix{nt.GetHelixParameters(bz)};
    std::array<double, 6> phelix{pt.GetHelixParameters(bz)};

    // Negative track parameters in XY
    auto [xNegCenter, yNegCenter] = nt.GetHelixCenter(bz);
    double NegRadius{std::abs(1. / nhelix[4])};

    // Positive track parameters in XY
    auto [xPosCenter, yPosCenter] = pt.GetHelixCenter(bz);
    double PosRadius{std::abs(1. / phelix[4])};

    // Define convenient coordinate system
    // Logical zero: position of negative center
    double ux{xPosCenter - xNegCenter};
    double uy{yPosCenter - yNegCenter};

    // Check center-to-center distance
    double lDist{std::sqrt(std::pow(xNegCenter - xPosCenter, 2) + std::pow(yNegCenter - yPosCenter, 2))};
    // Normalize ux, uz to unit vector
    ux /= lDist;
    uy /= lDist;

    // Calculate perpendicular vector (normalized)
    double vx{-uy};
    double vy{+ux};

    double lPreprocessDCAxy{1E3};  // define outside scope
    lPreprocessxp = pt.GetX();     // start at current location
    lPreprocessxn = nt.GetX();     // start at current location

    // Pre-optimization in the XY plane: cases considered here
    //  Case 1: Circles do not touch, centers far away
    //          (D > R1 + R2)
    //  Case 2: Circles touch, centers at reasonable distance wrt D
    //          (D < R1 + R2) && (D > |R1-R2|)
    //  Case 3: Circles do not touch, one inside the other
    //          (D < |R1-R2|)
    //
    //  Cases 1 and 2 are treated. Case 3 is not treated (unlikely
    //  to be a problem with unlike-sign charged tracks): brute
    //  force minimization takes place in any case

    // Case 1: distance bigger than sum of radii ("gamma-like")
    if (lDist > NegRadius + PosRadius) {
        // re-position tracks along the center-to-center axis
        // -- negative track
        double xNegOptPosition{xNegCenter + NegRadius * ux};
        double yNegOptPosition{yNegCenter + NegRadius * uy};
        double csNeg{std::cos(nt.GetAlpha())};
        double snNeg{std::sin(nt.GetAlpha())};
        double xThisNeg{xNegOptPosition * csNeg + yNegOptPosition * snNeg};

        // -- positive track
        double xPosOptPosition{xPosCenter - PosRadius * ux};
        double yPosOptPosition{yPosCenter - PosRadius * uy};
        double csPos{std::cos(pt.GetAlpha())};
        double snPos{std::sin(pt.GetAlpha())};
        double xThisPos{xPosOptPosition * csPos + yPosOptPosition * snPos};

        if (xThisNeg < lMaximumX && xThisPos < lMaximumX && xThisNeg > lMinimumX && xThisPos > lMinimumX) {
            std::array<double, 3> lCase1NegR{};
            std::array<double, 3> lCase1PosR{};
            if (nt.GetXYZAt(xThisNeg, bz, lCase1NegR) && pt.GetXYZAt(xThisPos, bz, lCase1PosR)) {
                lPreprocessDCAxy = std::sqrt(std::pow(lCase1NegR[0] - lCase1PosR[0], 2) + std::pow(lCase1NegR[1] - lCase1PosR[1], 2) +
                                             std::pow(lCase1NegR[2] - lCase1PosR[2], 2));
                // Pass coordinates
                if (lPreprocessDCAxy < 999) {
                    lPreprocessxp = xThisPos;
                    lPreprocessxn = xThisNeg;
                }
            }
        }
    }

    // Case 2: distance smaller than sum of radii (cowboy/sailor configs)
    if (lDist > std::abs(NegRadius - PosRadius) && lDist < NegRadius + PosRadius) {

        // calculate coordinate for radical line
        double lRadical{(lDist * lDist - PosRadius * PosRadius + NegRadius * NegRadius) / (2 * lDist)};

        // calculate absolute displacement from center-to-center axis
        double lDisplace{(0.5 / lDist) * std::sqrt((-lDist + PosRadius - NegRadius) * (-lDist - PosRadius + NegRadius) *
                                                   (-lDist + PosRadius + NegRadius) * (lDist + PosRadius + NegRadius))};

        // 3D distances in the two cases studied (prefer smallest)
        double lCase2aDCA{1E3};
        double lCase2bDCA{1E3};

        // 2 cases: positive and negative displacement
        std::array<double, 2> xNegOptPosition{};
        std::array<double, 2> yNegOptPosition{};
        std::array<double, 2> xPosOptPosition{};
        std::array<double, 2> yPosOptPosition{};
        std::array<double, 2> xThisNeg{};
        std::array<double, 2> xThisPos{};

        auto [snNeg, csNeg] = Math::sincos(nt.GetAlpha());
        auto [snPos, csPos] = Math::sincos(pt.GetAlpha());

        // Case 2a: Positive displacement along v vector
        // Re-position negative track
        xNegOptPosition[0] = xNegCenter + lRadical * ux + lDisplace * vx;
        yNegOptPosition[0] = yNegCenter + lRadical * uy + lDisplace * vy;
        xThisNeg[0] = xNegOptPosition[0] * csNeg + yNegOptPosition[0] * snNeg;
        // Re-position positive track
        xPosOptPosition[0] = xNegCenter + lRadical * ux + lDisplace * vx;
        yPosOptPosition[0] = yNegCenter + lRadical * uy + lDisplace * vy;
        xThisPos[0] = xPosOptPosition[0] * csPos + yPosOptPosition[0] * snPos;

        // Case 2b: Negative displacement along v vector
        // Re-position negative track
        xNegOptPosition[1] = xNegCenter + lRadical * ux - lDisplace * vx;
        yNegOptPosition[1] = yNegCenter + lRadical * uy - lDisplace * vy;
        xThisNeg[1] = xNegOptPosition[1] * csNeg + yNegOptPosition[1] * snNeg;
        // Re-position positive track
        xPosOptPosition[1] = xNegCenter + lRadical * ux - lDisplace * vx;
        yPosOptPosition[1] = yNegCenter + lRadical * uy - lDisplace * vy;
        xThisPos[1] = xPosOptPosition[1] * csPos + yPosOptPosition[1] * snPos;

        // Case 2a
        if (xThisNeg[0] < lMaximumX && xThisPos[0] < lMaximumX && xThisNeg[0] > lMinimumX && xThisPos[0] > lMinimumX) {
            std::array<double, 3> lCase2aNegR{};
            std::array<double, 3> lCase2aPosR{};
            if (nt.GetXYZAt(xThisNeg[0], bz, lCase2aNegR) && pt.GetXYZAt(xThisPos[0], bz, lCase2aPosR)) {
                lCase2aDCA = std::sqrt(std::pow(lCase2aNegR[0] - lCase2aPosR[0], 2) + std::pow(lCase2aNegR[1] - lCase2aPosR[1], 2) +
                                       std::pow(lCase2aNegR[2] - lCase2aPosR[2], 2));
            }
        }

        // Case 2b
        if (xThisNeg[1] < lMaximumX && xThisPos[1] < lMaximumX && xThisNeg[1] > lMinimumX && xThisPos[1] > lMinimumX) {
            std::array<double, 3> lCase2bNegR{};
            std::array<double, 3> lCase2bPosR{};
            if (nt.GetXYZAt(xThisNeg[1], bz, lCase2bNegR) && pt.GetXYZAt(xThisPos[1], bz, lCase2bPosR)) {
                lCase2bDCA = std::sqrt(std::pow(lCase2bNegR[0] - lCase2bPosR[0], 2) + std::pow(lCase2bNegR[1] - lCase2bPosR[1], 2) +
                                       std::pow(lCase2bNegR[2] - lCase2bPosR[2], 2));
            }
        }

        // Minor detail: all things being equal, prefer closest X
        double lCase2aSumX{xThisPos[0] + xThisNeg[0]};
        double lCase2bSumX{xThisPos[1] + xThisNeg[1]};

        double lDCAxySmallestR{lCase2aDCA};
        double lxpSmallestR{xThisPos[0]};
        double lxnSmallestR{xThisNeg[0]};

        double lDCAxyLargestR{lCase2bDCA};
        double lxpLargestR{xThisPos[1]};
        double lxnLargestR{xThisNeg[1]};

        if (lCase2bSumX + 1E-6 < lCase2aSumX) {
            lDCAxySmallestR = lCase2bDCA;
            lxpSmallestR = xThisPos[1];
            lxnSmallestR = xThisNeg[1];
            lDCAxyLargestR = lCase2aDCA;
            lxpLargestR = xThisPos[0];
            lxnLargestR = xThisNeg[0];
        }

        // Pass conclusion to lPreprocess variables
        lPreprocessDCAxy = lDCAxySmallestR;
        lPreprocessxp = lxpSmallestR;
        lPreprocessxn = lxnSmallestR;
        if (lDCAxyLargestR + 1E-6 < lDCAxySmallestR) {  // beware epsilon: numerical calculations are unstable here
            lPreprocessDCAxy = lDCAxyLargestR;
            lPreprocessxp = lxpLargestR;
            lPreprocessxn = lxnLargestR;
        }

        // Protection against something too crazy
        if (lPreprocessDCAxy > 999) {
            lPreprocessxp = pt.GetX();  // start at current location
            lPreprocessxn = nt.GetX();  // start at current location
        }
    }

    // End of preprocessing stage
    // at this point lPreprocessxp, lPreprocessxn are already good starting points: update helixparams
    if (lPreprocessDCAxy < 999) return true;
    return false;
}

// Rewrite of `AliExternalTrackParam::GetDCA()`.
// Return the (weighted!) distance of closest approach between this track and the track `p`.
// Other returned values:
// `x_this`, `x_t` - coordinates of tracks' reference planes at the DCA
static double Preoptimize_Numerically(const ALICE::Track &p1, const ALICE::Track &p2, double &x_p1, double &x_p2, double bz) {

    double dy2{p1.GetSigmaY2() + p2.GetSigmaY2()};
    double dz2{p1.GetSigmaZ2() + p2.GetSigmaZ2()};
    double dx2{dy2};

    std::array<double, 6> h1{p1.GetHelixParameters(bz)};
    std::array<double, 8> params1{};
    std::copy(h1.begin(), h1.end(), params1.begin());
    std::tie(params1[6], params1[7]) = Math::sincos(params1[2]);

    std::array<double, 6> h2{p2.GetHelixParameters(bz)};
    std::array<double, 8> params2{};
    std::copy(h2.begin(), h2.end(), params2.begin());
    std::tie(params2[6], params2[7]) = Math::sincos(params2[2]);

    std::array<double, 3> r1{};
    std::array<double, 3> g1{};
    std::array<double, 3> gg1{};
    double t1{0.};
    Evaluate(params1, t1, r1, g1, gg1);

    std::array<double, 3> r2{};
    std::array<double, 3> g2{};
    std::array<double, 3> gg2{};
    double t2{0.};
    Evaluate(params2, t2, r2, g2, gg2);

    double dx{r2[0] - r1[0]};
    double dy{r2[1] - r1[1]};
    double dz{r2[2] - r1[2]};
    double dm{dx * dx / dx2 + dy * dy / dy2 + dz * dz / dz2};

    int max{27};
    while (max--) {
        double gt1{-(dx * g1[0] / dx2 + dy * g1[1] / dy2 + dz * g1[2] / dz2)};
        double gt2{+(dx * g2[0] / dx2 + dy * g2[1] / dy2 + dz * g2[2] / dz2)};
        double h11{(g1[0] * g1[0] - dx * gg1[0]) / dx2 + (g1[1] * g1[1] - dy * gg1[1]) / dy2 + (g1[2] * g1[2] - dz * gg1[2]) / dz2};
        double h22{(g2[0] * g2[0] + dx * gg2[0]) / dx2 + (g2[1] * g2[1] + dy * gg2[1]) / dy2 + (g2[2] * g2[2] + dz * gg2[2]) / dz2};
        double h12{-(g1[0] * g2[0] / dx2 + g1[1] * g2[1] / dy2 + g1[2] * g2[2] / dz2)};

        double det{h11 * h22 - h12 * h12};

        double dt1;
        double dt2;
        if (std::abs(det) < 1.E-33) {
            //(quasi)singular Hessian
            dt1 = -gt1;
            dt2 = -gt2;
        } else {
            dt1 = -(gt1 * h22 - gt2 * h12) / det;
            dt2 = -(h11 * gt2 - h12 * gt1) / det;
        }

        if (dt1 * gt1 + dt2 * gt2 > 0.) {
            dt1 = -dt1;
            dt2 = -dt2;
        }

        // check delta(phase1) ?
        // check delta(phase2) ?

        if (std::abs(dt1) / (std::abs(t1) + 1.E-3) < 1.E-4)
            if (std::abs(dt2) / (std::abs(t2) + 1.E-3) < 1.E-4) {
#ifdef T2S_L2_DEBUG
                if ((gt1 * gt1 + gt2 * gt2) > 1.E-4 / dy2 / dy2) {
                    std::cout << __FUNCTION__ << " :: stopped at not a stationary point !" << '\n';
                }
#endif
                double lmb{h11 + h22};
                lmb = lmb - std::sqrt(lmb * lmb - 4 * det);
#ifdef T2S_L2_DEBUG
                if (lmb < 0.) {
                    std::cout << __FUNCTION__ << " :: stopped at not a minimum !" << '\n';
                }
#endif
                break;
            }

        double dd{dm};
        for (int div{1};; div *= 2) {
            Evaluate(params1, t1 + dt1, r1, g1, gg1);
            Evaluate(params2, t2 + dt2, r2, g2, gg2);
            dx = r2[0] - r1[0];
            dy = r2[1] - r1[1];
            dz = r2[2] - r1[2];
            dd = dx * dx / dx2 + dy * dy / dy2 + dz * dz / dz2;
            if (dd < dm) break;
            dt1 *= 0.5;
            dt2 *= 0.5;
            if (div > 512) {
#ifdef T2S_L2_DEBUG
                std::cout << __FUNCTION__ << " :: overshoot !" << '\n';
#endif
                break;
            }
        }
        dm = dd;

        t1 += dt1;
        t2 += dt2;
    }

#ifdef T2S_L2_DEBUG
    if (max <= 0) {
        std::cout << __FUNCTION__ << " :: too many iterations !" << '\n';
    }
#endif

    auto [sn_this, cs_this] = Math::sincos(p1.GetAlpha());
    x_p1 = r1[0] * cs_this + r1[1] * sn_this;

    auto [sn_p, cs_p] = Math::sincos(p2.GetAlpha());
    x_p2 = r2[0] * cs_p + r2[1] * sn_p;

    return std::sqrt(dm * std::sqrt(dy2 * dz2));
}

}  // namespace ALICE::Vertexer

#endif  // ALICE_VERTEXER_HXX
