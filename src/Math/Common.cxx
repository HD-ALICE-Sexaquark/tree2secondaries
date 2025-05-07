#include "Math/VectorUtil.h"

#include "Math/Common.hxx"

using namespace ROOT::Math;

namespace Tree2Secondaries::Math {

// Calculate the cosine of the pointing angle of a particle with momentum Px,Py,Pz and vertex X,Y,Z w.r.t. to a reference point.
double CosinePointingAngle(const XYZVector& mom_v0, const XYZPoint& pos_v0, const XYZPoint& pos_ref) {
    return std::cos(VectorUtil::Angle(mom_v0, pos_v0 - pos_ref));
}

// Overload of `CosinePointingAngle(...)`, using scalars instead of vectors.
double CosinePointingAngle(double v0_px, double v0_py, double v0_pz,  //
                           double v0_x, double v0_y, double v0_z,     //
                           double ref_x, double ref_y, double ref_z) {
    XYZVector MomV0(v0_px, v0_py, v0_pz);
    XYZPoint PosV0(v0_x, v0_y, v0_z);
    XYZPoint PosRef(ref_x, ref_y, ref_z);
    return CosinePointingAngle(MomV0, PosV0, PosRef);
}

// Calculate Armenteros-Podolanski qT.
// Based on https://github.com/alisw/AliRoot (STEER/ESD/AliESDv0::PtArmV0())
double ArmenterosQt(const XYZVector& mom_v0, const XYZVector& mom_dau) { return VectorUtil::Perp(mom_v0, mom_dau); }

// Overload of `ArmenterosQt(...)`
double ArmenterosQt(double v0_px, double v0_py, double v0_pz,  //
                    double dau_px, double dau_py, double dau_pz) {
    XYZVector MomV0(v0_px, v0_py, v0_pz);
    XYZVector MomDau(dau_px, dau_py, dau_pz);
    return ArmenterosQt(MomV0, MomDau);
}

// Calculate Armenteros-Podolanski alpha.
// Based on https://github.com/alisw/AliRoot (STEER/ESD/AliESDv0::AlphaV0())
double ArmenterosAlpha(const XYZVector& mom_v0, const XYZVector& mom_neg, const XYZVector& mom_pos) {
    double lQlNeg = mom_neg.Dot(mom_v0) / mom_v0.R();
    double lQlPos = mom_pos.Dot(mom_v0) / mom_v0.R();
    if (std::abs(lQlPos + lQlNeg) < 1E-6) return 2.;  // protection
    return (lQlPos - lQlNeg) / (lQlPos + lQlNeg);
}

// Overload of `ArmenterosAlpha()`.
double ArmenterosAlpha(double v0_px, double v0_py, double v0_pz,     //
                       double neg_px, double neg_py, double neg_pz,  //
                       double pos_px, double pos_py, double pos_pz) {
    XYZVector MomV0(v0_px, v0_py, v0_pz);
    XYZVector MomNeg(neg_px, neg_py, neg_pz);
    XYZVector MomPos(pos_px, pos_py, pos_pz);
    return ArmenterosAlpha(MomV0, MomNeg, MomPos);
}

}  // namespace Tree2Secondaries::Math
