#ifndef T2S_PACKAGER_CUTS_HXX
#define T2S_PACKAGER_CUTS_HXX

namespace Tree2Secondaries::Cuts {

namespace Track {
constexpr double AbsMax_PID_NSigma{3.};
}  // namespace Track

namespace Lambda {
// kinematics //
constexpr double Min_Pt{1.};
constexpr double Min_Mass{1.08};
constexpr double Max_Mass{1.16};
constexpr double AbsMax_Eta{0.9};
constexpr double Min_CPAwrtPV{0.35};
constexpr double Max_CPAwrtPV{0.9};
constexpr double Min_DCAwrtPV{4.};
constexpr double AbsMax_ArmQtOverAlpha{0.2};
// geometric //
constexpr double AbsMax_Zv{50.};
constexpr double Min_Radius{45.};
constexpr double Max_Radius{140.};
constexpr double Max_DCAnegV0{2.};
constexpr double Max_DCAposV0{2.};
constexpr double Max_DCAbtwDau{2.};
}  // namespace Lambda

namespace KaonZeroShort {
// kinematics //
constexpr double Min_Pt{1.};
constexpr double Min_Mass{0.475};
constexpr double Max_Mass{0.525};
constexpr double AbsMax_Eta{0.8};
constexpr double Min_CPAwrtPV{0.35};
constexpr double Max_CPAwrtPV{0.9};
constexpr double Min_DCAwrtPV{4.};
// geometric //
constexpr double AbsMax_Zv{50.};
constexpr double Min_Radius{20.};
constexpr double Max_Radius{180.};
constexpr double Max_DCAnegV0{2.};
constexpr double Max_DCAposV0{2.};
constexpr double Max_DCAbtwDau{2.};
}  // namespace KaonZeroShort

}  // namespace Tree2Secondaries::Cuts

#endif  // T2S_PACKAGER_CUTS_HXX
