#pragma once

namespace Tree2Secondaries::Cuts {

namespace Track {
constexpr double AbsMax_PID_NSigma{3.};
}  // namespace Track

namespace Lambda {
// kinematics //
constexpr double Min_Pt{0.5};
constexpr double Min_Mass{1.0};
constexpr double Max_Mass{1.2};
constexpr double AbsMax_Eta{1.};
constexpr double Min_CPAwrtPV{0.35};
constexpr double Max_CPAwrtPV{0.99};
constexpr double Min_DCAwrtPV{10.};
constexpr double AbsMax_ArmQtOverAlpha{0.2};
// geometric //
constexpr double AbsMax_Zv{100.};
constexpr double Min_Radius2D{20.};
constexpr double Max_Radius2D{180.};
constexpr double Max_DCAnegV0{10.};
constexpr double Max_DCAposV0{10.};
constexpr double Max_DCAbtwDau{10.};
}  // namespace Lambda

namespace KaonZeroShort {
// kinematics //
constexpr double Min_Pt{0.5};
constexpr double Min_Mass{0.4};
constexpr double Max_Mass{0.6};
constexpr double AbsMax_Eta{1.};
constexpr double Min_CPAwrtPV{0.35};
constexpr double Max_CPAwrtPV{0.99};
constexpr double Min_DCAwrtPV{10.};
// geometric //
constexpr double AbsMax_Zv{50.};
constexpr double Min_Radius2D{20.};
constexpr double Max_Radius2D{180.};
constexpr double Max_DCAnegV0{10.};
constexpr double Max_DCAposV0{10.};
constexpr double Max_DCAbtwDau{10.};
}  // namespace KaonZeroShort

}  // namespace Tree2Secondaries::Cuts
