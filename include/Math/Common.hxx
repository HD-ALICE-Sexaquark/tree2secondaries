#ifndef T2S_MATH_COMMON_HXX
#define T2S_MATH_COMMON_HXX

#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/VectorUtil.h"

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
// Based on https://github.com/alisw/AliRoot (STEER/ESD/AliESDv0::PtArmV0())
inline double ArmenterosQt(const ROOT::Math::XYZVector& mom_v0, const ROOT::Math::XYZVector& mom_dau) {
    return ROOT::Math::VectorUtil::Perp(mom_v0, mom_dau);
}

// Overload of `ArmenterosQt(...)`
inline double ArmenterosQt(double v0_px, double v0_py, double v0_pz,  //
                           double dau_px, double dau_py, double dau_pz) {
    return ArmenterosQt({v0_px, v0_py, v0_pz}, {dau_px, dau_py, dau_pz});
}

// Calculate Armenteros-Podolanski alpha.
// Based on https://github.com/alisw/AliRoot (STEER/ESD/AliESDv0::AlphaV0())
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

// Calculate the half-point between two points in 3D space.
inline ROOT::Math::XYZPoint MiddlePoint(const ROOT::Math::XYZPoint& p1, const ROOT::Math::XYZPoint& p2) {
    return {0.5 * (p1.X() + p2.X()), 0.5 * (p1.Y() + p2.Y()), 0.5 * (p1.Z() + p2.Z())};
}

}  // namespace Tree2Secondaries::Math

#endif  // T2S_MATH_COMMON_HXX
