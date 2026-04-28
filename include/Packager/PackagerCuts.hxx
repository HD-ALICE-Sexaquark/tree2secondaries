#pragma once

namespace Tree2Secondaries::Cuts {

namespace Proton {
constexpr double AbsMax_NSigmaProton = 3.;
// constexpr double AbsMin_DCAxy = 30.;  // PENDING
}  // namespace Proton

namespace Kaon {
constexpr double AbsMax_NSigmaKaon = 3.;
}  // namespace Kaon

namespace Pion {
constexpr double AbsMax_NSigmaPion = 3.;
// constexpr double AbsMin_DCAxy = 30.;  // PENDING
}  // namespace Pion

namespace Lambda {
// kinematics //
// constexpr double Min_Pt = 1.25;    // (GeV/c)
constexpr double Min_Mass = 1.09;  // (GeV/c^2)
constexpr double Max_Mass = 1.14;  // (GeV/c^2)
constexpr double AbsMax_Rapidity = 1.;
constexpr double Min_CPAwrtPV = 0.3;
constexpr double Max_CPAwrtPV = 0.95;  // PENDING: I want 0.95!
constexpr double Min_DCAwrtPV = 40.;   // (cm) // PENDING: I want 45!
constexpr double AbsMax_ArmQtOverAlpha = 0.2;
// geometric //
constexpr double Min_Radius2D = 65.;   // (cm)
constexpr double Max_DCAnegV0 = 0.2;   // (cm)
constexpr double Max_DCAposV0 = 0.2;   // (cm)
constexpr double Max_DCAbtwDau = 0.4;  // (cm)
}  // namespace Lambda

namespace KaonZeroShort {
// kinematics //
// constexpr double Min_Pt = 0.75;    // (GeV/c)
constexpr double Min_Mass = 0.47;  // (GeV/c^2)
constexpr double Max_Mass = 0.53;  // (GeV/c^2)
constexpr double AbsMax_Rapidity = 1.;
constexpr double Min_CPAwrtPV = 0.1;
constexpr double Max_CPAwrtPV = 0.95;  // PENDING: I want 0.95!
constexpr double Min_DCAwrtPV = 40.;   // (cm) // PENDING: I want 45!
// geometric //
constexpr double Min_Radius2D = 50.;    // (cm)
constexpr double Max_DCAnegV0 = 0.2;    // (cm)
constexpr double Max_DCAposV0 = 0.2;    // (cm)
constexpr double Max_DCAbtwDau = 0.25;  // (cm)
}  // namespace KaonZeroShort

}  // namespace Tree2Secondaries::Cuts
