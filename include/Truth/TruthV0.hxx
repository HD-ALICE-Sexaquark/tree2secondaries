#pragma once

#include "Math/Constants.hxx"
#include "Truth/TruthTrack.hxx"
#include "View/MC/ViewMcParticle.hxx"

namespace Tree2Secondaries::Truth::V0 {

[[nodiscard]] inline bool IsTrue(const View::MC::V0& v0, PID_V0 pid) {
    return Const::V0_PdgCode[pid] == v0.PdgCode() &&             //
           Track::IsTrue(v0.Neg, Const::V0_NegativePID[pid]) &&  //
           Track::IsTrue(v0.Pos, Const::V0_PositivePID[pid]);
}
[[nodiscard]] inline bool IsSignal(const View::MC::V0& v0, PID_V0 pid) {  //
    return IsTrue(v0, pid) && v0.Generator() == 2 &&
           Track::ReactionID(v0.Neg, v0, Const::V0_NegativePID[pid]) == Track::ReactionID(v0.Pos, v0, Const::V0_PositivePID[pid]);
}
[[nodiscard]] inline bool IsSecondary(const View::MC::V0& v0, PID_V0 pid) {  //
    return v0.IsSecFromMat() || v0.IsSecFromWeak() || IsSignal(v0, pid);
}
[[nodiscard]] inline int ReactionID(const View::MC::V0& v0, PID_V0 pid) {  //
    return IsSignal(v0, pid) ? v0.Status() : Const::DummyInt;
}
[[nodiscard]] inline bool IsHybrid(const View::MC::V0& v0, PID_V0 pid) {
    return !IsSignal(v0, pid) && ((Track::IsSignal(v0.Neg, Const::V0_NegativePID[pid]) && !Track::IsSignal(v0.Pos, Const::V0_PositivePID[pid])) ||
                                  (!Track::IsSignal(v0.Neg, Const::V0_NegativePID[pid]) && Track::IsSignal(v0.Pos, Const::V0_PositivePID[pid])));
}

[[nodiscard]] inline float DecayX(const View::MC::V0& v0) { return v0.Neg.X(); }
[[nodiscard]] inline float DecayY(const View::MC::V0& v0) { return v0.Neg.Y(); }
[[nodiscard]] inline float DecayZ(const View::MC::V0& v0) { return v0.Neg.Z(); }

}  // namespace Tree2Secondaries::Truth::V0
