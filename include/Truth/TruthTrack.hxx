#pragma once

#include "Math/Constants.hxx"
#include "View/BaseView.hxx"
#include "View/MC/ViewMcParticle.hxx"

namespace Tree2Secondaries::Truth::Track {

[[nodiscard]] inline bool IsTrue(const View::MC::Particle& mc, PID_StableParticle pid) { return Const::Particle_PdgCode[pid] == mc.PdgCode(); }
[[nodiscard]] inline bool IsSignal(const View::MC::Particle& mc, PID_StableParticle pid) { return IsTrue(mc, pid) && mc.Generator() == 2; }
[[nodiscard]] inline bool IsSecondary(const View::MC::Particle& mc, PID_StableParticle pid) {
    return mc.IsSecFromMat() || mc.IsSecFromWeak() || IsSignal(mc, pid);
}

[[nodiscard]] inline int ReactionID(const View::MC::Particle& mc, const View::MC::Particle& mother, PID_StableParticle pid) {
    if (IsSignal(mc, pid)) {
        if (View::IsValid(mother)) {
            return mother.Status();
        } else {
            return mc.Status();
        }
    }
    return Const::DummyInt;
}

}  // namespace Tree2Secondaries::Truth::Track
