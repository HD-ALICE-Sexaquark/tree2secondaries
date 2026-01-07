#pragma once

#include <cmath>

#include "Math/Constants.hxx"
#include "View/MC/ViewMcParticle.hxx"

namespace Tree2Secondaries::Truth::Particle {

// As Track //

[[nodiscard]] inline bool AsTrack_IsTrue(const View::MC::Particle& view, EParticle pid) {  //
    return Const::Particle_PdgCode[pid] == view.PdgCode();
}
[[nodiscard]] inline bool AsTrack_IsSignal(const View::MC::Particle& view, EParticle pid) {  //
    return AsTrack_IsTrue(view, pid) && view.Generator() == 2;
}
[[nodiscard]] inline bool AsTrack_IsSecondary(const View::MC::Particle& view, EParticle pid) {  //
    return view.IsSecFromMat() || view.IsSecFromWeak() || AsTrack_IsSignal(view, pid);
}
[[nodiscard]] inline int AsTrack_ReactionID(const View::MC::Particle& view, const View::MC::Particle& mother, EParticle pid) {  //
    return AsTrack_IsSignal(view, pid) ? mother.Status() : Const::DummyInt;
}

// As V0 //

[[nodiscard]] inline bool AsV0_IsTrue(const View::MC::V0& v0_view, EParticle v0_pid) {
    return Const::Particle_PdgCode[v0_pid] == v0_view.PdgCode() &&        //
           AsTrack_IsTrue(v0_view.Neg, Const::V0_NegativePID[v0_pid]) &&  //
           AsTrack_IsTrue(v0_view.Pos, Const::V0_PositivePID[v0_pid]);
}
[[nodiscard]] inline bool AsV0_IsSignal(const View::MC::V0& v0_view, EParticle v0_pid) {  //
    return AsV0_IsTrue(v0_view, v0_pid) && v0_view.Generator() == 2 &&
           AsTrack_ReactionID(v0_view.Neg, v0_view, Const::V0_NegativePID[v0_pid]) ==
               AsTrack_ReactionID(v0_view.Pos, v0_view, Const::V0_PositivePID[v0_pid]);
}
[[nodiscard]] inline bool AsV0_IsSecondary(const View::MC::V0& v0_view, EParticle v0_pid) {  //
    return v0_view.IsSecFromMat() || v0_view.IsSecFromWeak() || AsV0_IsSignal(v0_view, v0_pid);
}
[[nodiscard]] inline int AsV0_ReactionID(const View::MC::V0& v0_view, EParticle v0_pid) {  //
    return AsV0_IsSignal(v0_view, v0_pid) ? v0_view.Status() : Const::DummyInt;
}
[[nodiscard]] inline bool AsV0_IsHybrid(const View::MC::V0& v0_view, EParticle v0_pid) {
    return !AsV0_IsSignal(v0_view, v0_pid) &&
           ((AsTrack_IsSignal(v0_view.Neg, Const::V0_NegativePID[v0_pid]) && !AsTrack_IsSignal(v0_view.Pos, Const::V0_PositivePID[v0_pid])) ||
            (!AsTrack_IsSignal(v0_view.Neg, Const::V0_NegativePID[v0_pid]) && AsTrack_IsSignal(v0_view.Pos, Const::V0_PositivePID[v0_pid])));
}

[[nodiscard]] inline float AsV0_DecayX(const View::MC::V0& v0_view) { return v0_view.IsValid() ? v0_view.Neg.X() : Const::DummyFloat; }
[[nodiscard]] inline float AsV0_DecayY(const View::MC::V0& v0_view) { return v0_view.IsValid() ? v0_view.Neg.Y() : Const::DummyFloat; }
[[nodiscard]] inline float AsV0_DecayZ(const View::MC::V0& v0_view) { return v0_view.IsValid() ? v0_view.Neg.Z() : Const::DummyFloat; }

}  // namespace Tree2Secondaries::Truth::Particle
