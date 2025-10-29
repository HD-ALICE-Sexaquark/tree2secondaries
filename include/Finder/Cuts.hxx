#pragma once

namespace Tree2Secondaries::Cuts {

namespace ChannelA {
// kinematics //
constexpr double AbsMax_Rapidity{0.7};
constexpr double Min_CPAwrtPV{0.9};
constexpr double Max_CPAwrtPV{1.};
constexpr double Min_MassMinusNucleon{0.};
constexpr double Max_MassMinusNucleon{6.};
// geometry //
constexpr double Min_Radius2D{20.};
constexpr double Max_Radius2D{150.};
constexpr double Max_DCALaSV{10.};
constexpr double Max_DCALaNegSV{10.};
constexpr double Max_DCALaPosSV{10.};
constexpr double Max_DCAK0SV{10.};
constexpr double Max_DCAK0NegSV{10.};
constexpr double Max_DCAK0PosSV{10.};
constexpr double Max_DCAbtwV0s{10.};
constexpr double Max_DecayLengthLa{100.};
constexpr double Max_DecayLengthK0{100.};
}  // namespace ChannelA

namespace ChannelD {
// kinematics //
constexpr double AbsMax_Rapidity{0.7};
constexpr double Min_CPAwrtPV{0.9};
constexpr double Max_CPAwrtPV{1.};
// geometry //
constexpr double Min_Radius2D{20.};
constexpr double Max_Radius2D{180.};
constexpr double Max_DCALaSV{10.};
constexpr double Max_DCALaNegSV{10.};
constexpr double Max_DCALaPosSV{10.};
constexpr double Max_DCAKaSV{10.};
constexpr double Max_DCAKaLa{10.};
constexpr double Max_DCALaNegKa{10.};
constexpr double Max_DCALaPosKa{10.};
}  // namespace ChannelD

}  // namespace Tree2Secondaries::Cuts
