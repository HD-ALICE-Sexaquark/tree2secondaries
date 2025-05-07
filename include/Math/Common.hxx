#ifndef T2S_MATH_COMMON_HXX
#define T2S_MATH_COMMON_HXX

#include "Math/Point3D.h"
#include "Math/Vector3D.h"

namespace Tree2Secondaries::Math {

double CosinePointingAngle(const ROOT::Math::XYZVector& mom_v0, const ROOT::Math::XYZPoint& pos_v0, const ROOT::Math::XYZPoint& pos_ref);
double CosinePointingAngle(double v0_px, double v0_py, double v0_pz, double v0_x, double v0_y, double v0_z, double ref_x, double ref_y, double ref_z);

double ArmenterosQt(const ROOT::Math::XYZVector& mom_v0, const ROOT::Math::XYZVector& mom_dau);
double ArmenterosQt(double v0_px, double v0_py, double v0_pz, double dau_px, double dau_py, double dau_pz);

double ArmenterosAlpha(const ROOT::Math::XYZVector& mom_v0, const ROOT::Math::XYZVector& mom_neg, const ROOT::Math::XYZVector& mom_pos);
double ArmenterosAlpha(double v0_px, double v0_py, double v0_pz, double neg_px, double neg_py, double neg_pz, double pos_px, double pos_py,
                       double pos_pz);

}  // namespace Tree2Secondaries::Math

#endif  // T2S_MATH_COMMON_HXX
