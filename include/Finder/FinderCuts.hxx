#pragma once

namespace Tree2Secondaries::Cuts {

namespace ChannelA {
// kinematics //
// constexpr double Min_CPAwrtPV{0.9};
// constexpr double Max_DCAwrtPV = 15.;  // (cm)
// geometry //
constexpr double Min_Radius2D = 40.;   // (cm)
constexpr double Max_DCALaSV = 10.;    // (cm)
constexpr double Max_DCAK0SV = 10.;    // (cm)
constexpr double Max_DCAbtwV0s = 10.;  // (cm)
// constexpr double Min_La_CPAwrtSV = 0.5;   // PENDING
// constexpr double Min_K0S_CPAwrtSV = 0.5;  // PENDING
// constexpr double Max_DecayLengthLa = 100.; // PENDING
// constexpr double Max_DecayLengthK0 = 100.; // PENDING
}  // namespace ChannelA

namespace ChannelD {
// kinematics //
constexpr double AbsMax_Rapidity = 0.7;
constexpr double Min_CPAwrtPV = 0.9;
constexpr double Max_CPAwrtPV = 1.;
// geometry //
constexpr double Min_Radius2D = 20.;    // (cm)
constexpr double Max_Radius2D = 180.;   // (cm)
constexpr double Max_DCALaSV = 10.;     // (cm)
constexpr double Max_DCALaNegSV = 10.;  // (cm)
constexpr double Max_DCALaPosSV = 10.;  // (cm)
constexpr double Max_DCAKaSV = 10.;     // (cm)
constexpr double Max_DCAKaLa = 10.;     // (cm)
constexpr double Max_DCALaNegKa = 10.;  // (cm)
constexpr double Max_DCALaPosKa = 10.;  // (cm)
}  // namespace ChannelD

}  // namespace Tree2Secondaries::Cuts
