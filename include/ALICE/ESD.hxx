#ifndef ALICE_ESD_HXX
#define ALICE_ESD_HXX

#include <array>
#include <cmath>
#include <tuple>

#include <TMatrixD.h>

#include "ALICE/Constants.hxx"
#include "ALICE/Math.hxx"

namespace ALICE {

// `ALICE::Track`
// (based on `AliRoot/STEER/STEERBase/AliExternalTrackParam`)
class Track {

   public:
    Track() = delete;
    Track(double x, const std::array<double, 5> &params, const std::array<double, 15> &cov, double alpha, int charge)
        : fC{cov}, fP{params}, fX{x}, fAlpha{alpha}, fCharge{charge} {};

    double GetAlpha() const { return fAlpha; };
    int Charge() const { return fCharge; };
    double GetX() const { return fX; }
    double GetY() const { return fP[0]; }
    double GetZ() const { return fP[1]; }
    double GetSigmaY2() const { return fC[0]; };
    double GetSigmaZY() const { return fC[1]; };
    double GetSigmaZ2() const { return fC[2]; };
    double GetSnp() const { return fP[2]; };
    double GetTgl() const { return fP[3]; };

    // Rewrite of `AliVParticle::Local2GlobalPosition()`.
    // Perform local->global transformation of the track position.
    // When called, the arguments are:
    //    r[0] = local x
    //    r[1] = local y
    //    r[2] = local z
    //   alpha - rotation angle
    // The result is returned as:
    //    r[0] = global x
    //    r[1] = global y
    //    r[2] = global z
    static bool Local2GlobalPosition(std::array<double, 3> &r, double alpha) {
        auto [sn, cs] = Math::sincos(alpha);
        double x{r[0]};
        r[0] = x * cs - r[1] * sn;
        r[1] = x * sn + r[1] * cs;
        return true;
    }

    // Rewrite of `AliVParticle::Local2GlobalMomentum()`.
    // This function performs local->global transformation of the track momentum.
    // When called, the arguments are:
    //    p[0] = 1/pt * charge of the track;
    //    p[1] = sine of local azim. angle of the track momentum;
    //    p[2] = tangent of the track momentum dip angle;
    //   alpha - rotation angle.
    // The result is returned as:
    //    p[0] = px
    //    p[1] = py
    //    p[2] = pz
    // Results for (nearly) straight tracks are meaningless !
    static bool Local2GlobalMomentum(std::array<double, 3> &p, double alpha) {
        if (std::abs(p[0]) <= Const::AlmostZero) return false;
        if (std::abs(p[1]) > 1.) return false;

        double pt{1. / std::abs(p[0])};
        auto [sn, cs] = Math::sincos(alpha);
        double r{std::sqrt((1. - p[1]) * (1. + p[1]))};
        p[0] = pt * (r * cs - p[1] * sn);
        p[1] = pt * (p[1] * cs + r * sn);
        p[2] = pt * p[2];

        return true;
    }

    // Rewrite of `AliExternalTrackParam::GetPxPyPz()`
    // Return the global track momentum components
    // Results for (nearly) straight tracks are meaningless !
    bool GetPxPyPz(std::array<double, 3> &p) const {
        p[0] = fP[4];
        p[1] = fP[2];
        p[2] = fP[3];
        return Local2GlobalMomentum(p, fAlpha);
    }

    // Rewrite of `AliV0ReaderV1::GetHelixCenter()`.
    // Get center of the helix track parametrization.
    // (based on AliRoot/STEER/ESD/AliV0vertexer.cxx)
    std::pair<double, double> GetHelixCenter(double bz) const {

        double xpos{fP[5]};
        double ypos{fP[0]};
        double radius{std::abs(1. / fP[4])};
        double phi{fP[2]};
        if (phi < 0.) phi = phi + 2 * Const::Pi;

        phi -= Const::Pi / 2.;
        auto [ypoint, xpoint] = Math::sincos(phi);
        xpoint *= radius;
        ypoint *= radius;

        if (bz < 0 && fCharge > 0) {
            xpoint = -xpoint;
            ypoint = -ypoint;
        }

        if (bz > 0 && fCharge < 0) {
            xpoint = -xpoint;
            ypoint = -ypoint;
        }

        return {xpos + xpoint, ypos + ypoint};
    };

    // Rewrite of `AliExternalTrackParam::GetHelixParameters()`.
    // External track parameters -> helix parameters
    // `bz` - magnetic field (kG)
    std::array<double, 6> GetHelixParameters(double bz) const {
        auto [sn, cs] = Math::sincos(fAlpha);
        return {fX * sn + fP[0] * cs,       // y0
                fP[1],                      // z0
                std::asin(fP[2]) + fAlpha,  // phi0
                fP[3],                      // tgl
                GetC(bz),                   // C
                fX * cs - fP[0] * sn};      // x0
    };

    // Rewrite of `AliExternalTrackParam::GetXYZ()`.
    // Return the global track position.
    // Output arguments:
    // - `r`
    bool GetXYZ(std::array<double, 3> &r) const {
        r[0] = fX;
        r[1] = fP[0];
        r[2] = fP[1];
        return Local2GlobalPosition(r, fAlpha);
    }

    // Rewrite of `AliExternalTrackParam::GetXYZAt()`.
    // Return the global track position extrapolated to the radial position `x` (cm) under the magnetic field `bz` (kG)
    bool GetXYZAt(double x, double bz, std::array<double, 3> &r) const {
        double dx{x - fX};
        if (std::abs(dx) <= Const::AlmostZero) return GetXYZ(r);

        double crv{GetC(bz)};
        double x2r{crv * dx};
        double f1{fP[2]};
        double f2{f1 + dx * crv};

        if (std::abs(f1) >= 1.) return false;
        if (std::abs(f2) >= 1.) return false;

        double r1{std::sqrt((1. - f1) * (1. + f1))};
        double r2{std::sqrt((1. - f2) * (1. + f2))};
        double dy2dx{(f1 + f2) / (r1 + r2)};
        r[0] = x;
        r[1] = fP[0] + dx * dy2dx;
        if (std::abs(x2r) < 0.05) {
            r[2] = fP[1] + dx * (r2 + f2 * dy2dx) * fP[3];  // Thanks to Andrea & Peter
        } else {
            // for small dx/R the linear apporximation of the arc by the segment is OK,
            // but at large dx/R the error is very large and leads to incorrect Z propagation
            // angle traversed delta = 2*asin(dist_start_end / R / 2), hence the arc is: R*deltaPhi
            // The dist_start_end is obtained from sqrt(dx^2+dy^2) = x/(r1+r2)*sqrt(2+f1*f2+r1*r2)
            // Similarly, the rotation angle in linear in dx only for dx<<R
            double chord{dx * std::sqrt(1 + dy2dx * dy2dx)};  // distance from old position to new one
            double rot{2 * std::asin(0.5 * chord * crv)};     // angular difference seen from the circle center
            r[2] = fP[1] + rot / crv * fP[3];
        }

        return Local2GlobalPosition(r, fAlpha);
    };

    // Rewrite of `AliExternalTrackParam::GetC()`.
    double GetC(double bz) const { return fP[4] * bz * Const::B2C; }

    // Rewrite of `AliExternalTrackParam::GetLinearD()`.
    // Calculate the transverse impact parameter with respect to a point with global coordinates (`x`, `y`) neglecting the track curvature.
    double GetDCAxy_Linear(double x, double y) const {
        auto [sn, cs] = Math::sincos(fAlpha);
        double x_new{x * cs + y * sn};
        double y_new{-x * sn + y * cs};
        double d{(fX - x_new) * fP[2] - (fP[0] - y_new) * std::sqrt((1. - fP[2]) * (1. + fP[2]))};
        return -d;
    }

    // Rewrite of `AliExternalTrackParam::GetD()`.
    // Calculate the transverse impact parameter with respect to a point with global coordinates (`x`, `y`) in the magnetic field `bz` (kG)
    double GetDCAxy(double x, double y, double bz) const {
        if (std::abs(bz) < Const::AlmostZero) return GetDCAxy_Linear(x, y);
        double rp4{GetC(bz)};

        double xt{fX};
        double yt{fP[0]};

        auto [sn, cs] = Math::sincos(fAlpha);
        double a{x * cs + y * sn};
        double y_new{-x * sn + y * cs};
        double x_new{a};
        xt -= x_new;
        yt -= y_new;

        sn = rp4 * xt - fP[2];
        cs = rp4 * yt + std::sqrt((1. - fP[2]) * (1. + fP[2]));
        a = 2 * (xt * fP[2] - yt * std::sqrt((1. - fP[2]) * (1. + fP[2]))) - rp4 * (xt * xt + yt * yt);
        return -a / (1 + std::sqrt(sn * sn + cs * cs));
    }

    // Rewrite of `AliExternalTrackParam::PropagateTo()`.
    // Propagate this track to the plane X=xk (cm) in the field `bz` (kG)
    // NOTE: it modifies the state of `ALICE::Track`
    bool PropagateTo(double xk, double bz) {
        double dx{xk - fX};
        if (std::abs(dx) <= Const::AlmostZero) return true;

        double crv{GetC(bz)};
        if (std::abs(bz) < Const::AlmostZero) crv = 0.;

        double x2r{crv * dx};
        double f1{fP[2]};
        double f2{f1 + x2r};
        if (std::abs(f1) >= 1.) return false;
        if (std::abs(f2) >= 1.) return false;
        if (std::abs(fP[4]) < Const::AlmostZero) return false;

        double &fP0{fP[0]};
        double &fP1{fP[1]};
        double &fP2{fP[2]};
        double &fP3{fP[3]};
        double &fP4{fP[4]};
        double &fC00{fC[0]};
        double &fC10{fC[1]};
        double &fC11{fC[2]};
        double &fC20{fC[3]};
        double &fC21{fC[4]};
        double &fC22{fC[5]};
        double &fC30{fC[6]};
        double &fC31{fC[7]};
        double &fC32{fC[8]};
        double &fC33{fC[9]};
        double &fC40{fC[10]};
        double &fC41{fC[11]};
        double &fC42{fC[12]};
        double &fC43{fC[13]};
        double &fC44{fC[14]};

        double r1{std::sqrt((1. - f1) * (1. + f1))};
        double r2{std::sqrt((1. - f2) * (1. + f2))};
        if (std::abs(r1) < Const::AlmostZero) return false;
        if (std::abs(r2) < Const::AlmostZero) return false;

        fX = xk;
        double dy2dx{(f1 + f2) / (r1 + r2)};
        fP0 += dx * dy2dx;
        fP2 += x2r;
        if (std::abs(x2r) < 0.05) {
            fP1 += dx * (r2 + f2 * dy2dx) * fP3;  // Many thanks to P.Hristov !
        } else {
            // for small dx/R the linear apporximation of the arc by the segment is OK,
            // but at large dx/R the error is very large and leads to incorrect Z propagation
            // angle traversed delta = 2*asin(dist_start_end / R / 2), hence the arc is: R*deltaPhi
            // The dist_start_end is obtained from sqrt(dx^2+dy^2) = x/(r1+r2)*sqrt(2+f1*f2+r1*r2)
            //    double chord = dx*std::sqrt(1+dy2dx*dy2dx);   // distance from old position to new one
            //    double rot = 2*TMath::ASin(0.5*chord*crv); // angular difference seen from the circle center
            //    fP1 += rot/crv*fP3;
            double rot{std::asin(r1 * f2 - r2 * f1)};    // more economic version from Yura.
            if (f1 * f1 + f2 * f2 > 1 && f1 * f2 < 0) {  // special cases of large rotations or large abs angles
                if (f2 > 0)
                    rot = Const::Pi - rot;  //
                else
                    rot = -Const::Pi - rot;
            }
            fP1 += fP3 / crv * rot;
        }

        // f = F - 1
        double rinv{1. / r1};
        double r3inv{rinv * rinv * rinv};
        double f24{x2r / fP4};
        double f02{dx * r3inv};
        double f04{0.5 * f24 * f02};
        double f12{f02 * fP3 * f1};
        double f14{0.5 * f24 * f02 * fP3 * f1};
        double f13{dx * rinv};

        // b = C*ft
        double b00{f02 * fC20 + f04 * fC40};
        double b01{f12 * fC20 + f14 * fC40 + f13 * fC30};
        double b02{f24 * fC40};
        double b10{f02 * fC21 + f04 * fC41};
        double b11{f12 * fC21 + f14 * fC41 + f13 * fC31};
        double b12{f24 * fC41};
        double b20{f02 * fC22 + f04 * fC42};
        double b21{f12 * fC22 + f14 * fC42 + f13 * fC32};
        double b22{f24 * fC42};
        double b40{f02 * fC42 + f04 * fC44};
        double b41{f12 * fC42 + f14 * fC44 + f13 * fC43};
        double b42{f24 * fC44};
        double b30{f02 * fC32 + f04 * fC43};
        double b31{f12 * fC32 + f14 * fC43 + f13 * fC33};
        double b32{f24 * fC43};

        // a = f*b = f*C*ft
        double a00{f02 * b20 + f04 * b40};
        double a01{f02 * b21 + f04 * b41};
        double a02{f02 * b22 + f04 * b42};
        double a11{f12 * b21 + f14 * b41 + f13 * b31};
        double a12{f12 * b22 + f14 * b42 + f13 * b32};
        double a22{f24 * b42};

        // F*C*Ft = C + (b + bt + a)
        fC00 += b00 + b00 + a00;
        fC10 += b10 + b01 + a01;
        fC20 += b20 + b02 + a02;
        fC30 += b30;
        fC40 += b40;
        fC11 += b11 + b11 + a11;
        fC21 += b21 + b12 + a12;
        fC31 += b31;
        fC41 += b41;
        fC22 += b22 + b22 + a22;
        fC32 += b32;
        fC42 += b42;

        // CheckCovariance(); // PENDING

        return true;
    };

    // Rewrite of `AliExternalTrackParam::CheckCovariance()`.
    // Force the diagonal elements of the covariance matrix to be positive.
    // In case the diagonal element is bigger than the maximal allowed value, it is set to
    // the limit and the off-diagonal elements that correspond to it are set to zero.
    // NOTE: it modifies the state of the `ALICE::Track`.
    void CheckCovariance() {

        fC[0] = std::abs(fC[0]);
        if (fC[0] > Const::Max_C0) {
            double scl{std::sqrt(Const::Max_C0 / fC[0])};
            fC[0] = Const::Max_C0;
            fC[1] *= scl;
            fC[3] *= scl;
            fC[6] *= scl;
            fC[10] *= scl;
        }

        fC[2] = std::abs(fC[2]);
        if (fC[2] > Const::Max_C2) {
            double scl{std::sqrt(Const::Max_C2 / fC[2])};
            fC[2] = Const::Max_C2;
            fC[1] *= scl;
            fC[4] *= scl;
            fC[7] *= scl;
            fC[11] *= scl;
        }

        fC[5] = std::abs(fC[5]);
        if (fC[5] > Const::Max_C5) {
            double scl{std::sqrt(Const::Max_C5 / fC[5])};
            fC[5] = Const::Max_C5;
            fC[3] *= scl;
            fC[4] *= scl;
            fC[8] *= scl;
            fC[12] *= scl;
        }

        fC[9] = std::abs(fC[9]);
        if (fC[9] > Const::Max_C9) {
            double scl{std::sqrt(Const::Max_C9 / fC[9])};
            fC[9] = Const::Max_C9;
            fC[6] *= scl;
            fC[7] *= scl;
            fC[8] *= scl;
            fC[13] *= scl;
        }

        fC[14] = std::abs(fC[14]);
        if (fC[14] > Const::Max_C14) {
            double scl{std::sqrt(Const::Max_C14 / fC[14])};
            fC[14] = Const::Max_C14;
            fC[10] *= scl;
            fC[11] *= scl;
            fC[12] *= scl;
            fC[13] *= scl;
        }
    }

   private:
    // Covariance matrix of track parameters.
    std::array<double, 15> fC;

    // Track parameters.
    // `[0]` : local y-coordinate of a track (cm).
    // `[1]` : local z-coordinate of a track (cm).
    // `[2]` : local sine of the track momentum azimuthal angle.
    // `[3]` : tangent of the track momentum dip angle.
    // `[4]` : 1/pt (1/(GeV/c)).
    std::array<double, 5> fP;

    double fX;      // x-coordinate for the point of parametrization
    double fAlpha;  // local-to-global coord. system rotation angle
    int fCharge;
};

// `ALICE::V0`
// (based on `AliRoot/STEER/ESD/AliESDv0`)
class V0 {

   public:
    V0() = delete;

    // Rewrite of `AliESDv0::AliESDv0()`.
    // Main constructor.
    V0(const ALICE::Track &neg, const ALICE::Track &pos) : fParamN{neg}, fParamP{pos} {

        // get neg particle params
        double alpha{fParamN.GetAlpha()};
        auto [sn, cs] = Math::sincos(alpha);

        std::array<double, 3> tmp{};
        fParamN.GetPxPyPz(tmp);
        double px1{tmp[0]};
        double py1{tmp[1]};
        double pz1{tmp[2]};
        fParamN.GetXYZ(tmp);
        double x1{tmp[0]};
        double y1{tmp[1]};
        double z1{tmp[2]};
        double sx1{sn * sn * fParamN.GetSigmaY2() + Const::ResMisAlignPrec};
        double sy1{cs * cs * fParamN.GetSigmaY2() + Const::ResMisAlignPrec};

        // get pos particle params
        alpha = fParamP.GetAlpha();
        std::tie(sn, cs) = Math::sincos(alpha);
        fParamP.GetPxPyPz(tmp);
        double px2{tmp[0]};
        double py2{tmp[1]};
        double pz2{tmp[2]};
        fParamP.GetXYZ(tmp);
        double x2{tmp[0]};
        double y2{tmp[1]};
        double z2{tmp[2]};
        double sx2{sn * sn * fParamP.GetSigmaY2() + Const::ResMisAlignPrec};
        double sy2{cs * cs * fParamP.GetSigmaY2() + Const::ResMisAlignPrec};

        double sz1{fParamN.GetSigmaZ2()};
        double sz2{fParamP.GetSigmaZ2()};
        double wx1{sx2 / (sx1 + sx2)};
        double wx2{1. - wx1};
        double wy1{sy2 / (sy1 + sy2)};
        double wy2{1. - wy1};
        double wz1{sz2 / (sz1 + sz2)};
        double wz2{1. - wz1};

        // update v0 position
        fPos[0] = wx1 * x1 + wx2 * x2;
        fPos[1] = wy1 * y1 + wy2 * y2;
        fPos[2] = wz1 * z1 + wz2 * z2;

        fNmom[0] = px1;
        fNmom[1] = py1;
        fNmom[2] = pz1;
        fPmom[0] = px2;
        fPmom[1] = py2;
        fPmom[2] = pz2;

        /*
        // PENDING: will this become useful?
            for (int i{0}; i < 6; ++i) {
                fClusters[0][i] = 0;
                fClusters[1][i] = 0;
            }
            fNormDCAPrim[0] = fNormDCAPrim[1] = 0;
            for (int i{0}; i < 3; ++i) {
                fAngle[i] = 0;
            }
            for (int i{0}; i < 4; ++i) {
                fCausality[i] = 0;
            }
        */
    }

    double X() const { return fPos[0]; };
    double Y() const { return fPos[1]; };
    double Z() const { return fPos[2]; };
    // double Neg_PCA_XYZ() const { return; };
    // double Neg_PCA_XYZ() const { return; };
    // double Neg_PCA_XYZ() const { return; };
    // double Pos_PCA_XYZ() const { return; };
    // double Pos_PCA_XYZ() const { return; };
    // double Pos_PCA_XYZ() const { return; };
    // double Mass() const { return; };
    // double DCA_Daughters() const { return; };
    double Radius2D() const { return std::sqrt(fPos[0] * fPos[0] + fPos[1] * fPos[1]); };
    // double DCA_Neg_V0() const { return; };
    // double DCA_Pos_V0() const { return; };
    // double Pt() const { return; };
    // double Eta() const { return; };
    // double ArmenterosQt() const { return; };
    // double ArmenterosAlpha() const { return; };

    // Rewrite of `static AliESDv0/GetWeight()`.
    // Return the global weight matrix w = Transpose[G2P]*Inverse[Cpar]*G2P ,
    // where the matrix Cpar is the transverse part of the t covariance
    // in "parallel" system (i.e. the system with X axis parallel to momentum).
    // The matrix G2P performs the transformation global -> "parallel".
    static bool GetWeight(TMatrixD &w, const ALICE::Track &t) {

        double phi{t.GetAlpha() + std::asin(t.GetSnp())};
        auto [sp, cp] = Math::sincos(phi);

        double tgl{t.GetTgl()};
        double cl{1 / std::sqrt(1. + tgl * tgl)};
        double sl{tgl * cl};

        TMatrixD g2p(3, 3);  // global -> parallel
        g2p(0, 0) = cp * cl;
        g2p(0, 1) = sp * cl;
        g2p(0, 2) = sl;
        g2p(1, 0) = -sp;
        g2p(1, 1) = cp;
        g2p(1, 2) = 0.;
        g2p(2, 0) = -sl * cp;
        g2p(2, 1) = -sl * sp;
        g2p(2, 2) = cl;

        double alpha{t.GetAlpha()};
        auto [s, c] = Math::sincos(alpha);
        TMatrixD l2g(3, 3);  // local -> global
        l2g(0, 0) = c;
        l2g(0, 1) = -s;
        l2g(0, 2) = 0;
        l2g(1, 0) = s;
        l2g(1, 1) = c;
        l2g(1, 2) = 0;
        l2g(2, 0) = 0;
        l2g(2, 1) = 0;
        l2g(2, 2) = 1;

        double sy2{t.GetSigmaY2()};
        double syz{t.GetSigmaZY()};
        double sz2{t.GetSigmaZ2()};
        TMatrixD cvl(3, 3);  // local covariance
        cvl(0, 0) = 0;
        cvl(0, 1) = 0;
        cvl(0, 2) = 0;
        cvl(1, 0) = 0;
        cvl(1, 1) = sy2;
        cvl(1, 2) = syz;
        cvl(2, 0) = 0;
        cvl(2, 1) = syz;
        cvl(2, 2) = sz2;

        TMatrixD l2p(g2p, TMatrixD::kMult, l2g);
        TMatrixD cvp(3, 3);  // parallel covariance
        cvp = l2p * cvl * TMatrixD(TMatrixD::kTransposed, l2p);

        double det{cvp(1, 1) * cvp(2, 2) - cvp(1, 2) * cvp(2, 1)};
        if (std::abs(det) < Const::AlmostZero) return false;

        double m{100. * 100.};  // A large uncertainty in the momentum direction
        double eps{1 / m};
        TMatrixD u(3, 3);  // Inverse of the transverse part of the parallel covariance
        u(0, 0) = eps;
        u(0, 1) = 0;
        u(0, 2) = 0;
        u(1, 0) = 0;
        u(1, 1) = cvp(2, 2) / det;
        u(1, 2) = -cvp(2, 1) / det;
        u(2, 0) = 0;
        u(2, 1) = -cvp(1, 2) / det;
        u(2, 2) = cvp(1, 1) / det;

        w = TMatrixD(TMatrixD::kTransposed, g2p) * u * g2p;

        return true;
    }

    // Rewrite of `AliESDv0::Refit()`.
    // Refit this vertex
    int Refit() {

        fStatus = 0;

        // save the daughters' momenta
        fParamN.GetPxPyPz(fNmom);
        fParamP.GetPxPyPz(fPmom);

        // trivial estimation of the V0 vertex parameters
        std::array<double, 3> r1{};
        fParamN.GetXYZ(r1);
        std::array<double, 3> r2{};
        fParamP.GetXYZ(r2);

        for (int i{0}; i < 3; ++i) fPos[i] = 0.5 * (r1[i] + r2[i]);
        fPosCov[1] = fPosCov[3] = fPosCov[4] = 0.;
        fPosCov[0] = (fPos[0] - r1[0]) * (fPos[0] - r1[0]) / 12.;
        fPosCov[2] = (fPos[1] - r1[1]) * (fPos[1] - r1[1]) / 12.;
        fPosCov[5] = (fPos[2] - r1[2]) * (fPos[2] - r1[2]) / 12.;
        fChi2V0 = 12.;

        // try to improve the V0 vertex parameters
        TMatrixD w1(3, 3);
        if (!GetWeight(w1, fParamN)) return fStatus;
        TMatrixD w2(3, 3);
        if (!GetWeight(w2, fParamP)) return fStatus;
        TMatrixD cv(w1);
        cv += w2;
        cv.Invert();
        if (!cv.IsValid()) return fStatus;

        // covariance of the V0 vertex
        fPosCov[0] = cv(0, 0);
        fPosCov[1] = cv(1, 0);
        fPosCov[2] = cv(1, 1);
        fPosCov[3] = cv(2, 0);
        fPosCov[4] = cv(2, 1);
        fPosCov[5] = cv(2, 2);

        // position of the V0 vertex
        TMatrixD cw1(cv, TMatrixD::kMult, w1);
        for (int i{0}; i < 3; ++i) {
            fPos[i] = r2[i];
            for (int j{0}; j < 3; ++j) {
                fPos[i] += cw1(i, j) * (r1[j] - r2[j]);
            }
        }

        // chi2 of the v0 vertex
        fChi2V0 = 0.;
        std::array<double, 3> res1{r1[0] - fPos[0], r1[1] - fPos[1], r1[2] - fPos[2]};
        std::array<double, 3> res2{r2[0] - fPos[0], r2[1] - fPos[1], r2[2] - fPos[2]};
        for (int i{0}; i < 3; ++i) {
            for (int j{0}; j < 3; ++j) {
                fChi2V0 += res1[i] * res1[j] * w1(i, j) + res2[i] * res2[j] * w2(i, j);
            }
        }

        // successful fit
        fStatus = 1;

        return fStatus;
    }

    ALICE::Track fParamN;
    ALICE::Track fParamP;
    int fStatus{};

   private:
    std::array<double, 6> fPosCov{};  // covariance matrix of the vertex position
    std::array<double, 3> fPos{};     // V0's position (global)
    std::array<double, 3> fNmom{};    // momentum of the negative daughter (global)
    std::array<double, 3> fPmom{};    // momentum of the positive daughter (global)
    double fChi2V0{};                 // V0's chi2 value
};

}  // namespace ALICE

#endif  // ALICE_ESD_HXX
