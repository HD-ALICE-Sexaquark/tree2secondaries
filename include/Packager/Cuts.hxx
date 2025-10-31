#pragma once

namespace Tree2Secondaries::Cuts {

namespace Proton {
constexpr double AbsMax_NSigmaProton{2.};
constexpr double Min_TPCSignal{65.};
constexpr double Max_P{1.25};
constexpr double Max_DCAxy{5.};
}  // namespace Proton

namespace Kaon {
constexpr double AbsMax_NSigmaKaon{2.};
constexpr double Min_TPCSignal{60.};
constexpr double Max_DCAxy{5.};
}  // namespace Kaon

namespace Pion {
constexpr double AbsMax_NSigmaPion{2.5};
constexpr double Max_P{1.};
constexpr double Max_DCAxy{5.};
}  // namespace Pion

namespace Lambda {
// kinematics //
constexpr double Min_Pt{0.5};     // GeV/c
constexpr double Min_Mass{1.1};   // GeV/c^2
constexpr double Max_Mass{1.13};  // GeV/c^2
constexpr double AbsMax_Eta{1.};
constexpr double Min_CPAwrtPV{0.35};
constexpr double Max_CPAwrtPV{0.8};
constexpr double Min_DCAwrtPV{35.};  // cm
constexpr double AbsMax_ArmQtOverAlpha{0.2};
// geometric //
constexpr double AbsMax_Zv{40.};      // cm
constexpr double Min_Radius2D{45.};   // cm
constexpr double Max_Radius2D{180.};  // cm
constexpr double Max_DCAnegV0{0.5};   // cm
constexpr double Max_DCAposV0{0.5};   // cm
constexpr double Max_DCAbtwDau{0.5};  // cm
}  // namespace Lambda

namespace KaonZeroShort {
// kinematics //
constexpr double Min_Pt{1.};      // GeV/c
constexpr double Min_Mass{0.46};  // GeV/c^2
constexpr double Max_Mass{0.54};  // GeV/c^2
constexpr double AbsMax_Eta{1.};
constexpr double Min_CPAwrtPV{0.35};
constexpr double Max_CPAwrtPV{0.85};
constexpr double Min_DCAwrtPV{60.};  // cm
// geometric //
constexpr double AbsMax_Zv{50.};       // cm
constexpr double Min_Radius2D{50.};    // cm
constexpr double Max_Radius2D{180.};   // cm
constexpr double Max_DCAnegV0{0.25};   // cm
constexpr double Max_DCAposV0{0.25};   // cm
constexpr double Max_DCAbtwDau{0.25};  // cm
}  // namespace KaonZeroShort

}  // namespace Tree2Secondaries::Cuts
