#pragma once

#include <cmath>

#include "Math/Constants.hxx"
#include "View/MC/ViewMcInjected.hxx"
#include "View/MC/ViewMcPackedTrack.hxx"
#include "View/MC/ViewMcPackedV0.hxx"

namespace Tree2Secondaries::Truth::Sexaquark {

[[nodiscard]] inline float AsInjected_Energy(const View::MC::Injected& view, double mass) {  //
    return static_cast<float>(std::sqrt(view.P2() + mass * mass));
}

[[nodiscard]] inline float AsInjected_NucleonEnergy(const View::MC::Injected& view, EReactionChannel channel) {  //
    double mass{Const::Particle_Mass[Const::ReactionNucleonPID[channel]]};
    return static_cast<float>(std::sqrt(view.Nucleon_P2() + mass * mass));
}

// Channel A //

[[nodiscard]] inline double AsChannelA_AfterPx(const View::MC::ChannelA& a_view) {
    return a_view.IsValid() ? a_view.V0A.Px() + a_view.V0B.Px() : Const::DummyDouble;
}
[[nodiscard]] inline double AsChannelA_AfterPy(const View::MC::ChannelA& a_view) {
    return a_view.IsValid() ? a_view.V0A.Py() + a_view.V0B.Py() : Const::DummyDouble;
}
[[nodiscard]] inline double AsChannelA_AfterPz(const View::MC::ChannelA& a_view) {
    return a_view.IsValid() ? a_view.V0A.Pz() + a_view.V0B.Pz() : Const::DummyDouble;
}
[[nodiscard]] inline double AsChannelA_AfterE(const View::MC::ChannelA& a_view) {
    return a_view.IsValid() ? a_view.V0A.Energy() + a_view.V0B.Energy() : Const::DummyDouble;
}

[[nodiscard]] inline bool AsChannelA_IsSignal(const View::MC::ChannelA& a_view) {
    return a_view.IsValid() ? a_view.V0A.IsSignal() && a_view.V0B.IsSignal() : false;
}
[[nodiscard]] inline bool AsChannelA_IsHybrid(const View::MC::ChannelA& a_view) {
    // NOTE: don't check IsValid(), on purpose
    return !AsChannelA_IsSignal(a_view) && (a_view.V0A.IsHybrid() || a_view.V0B.IsHybrid() || (a_view.V0A.IsSignal() && !a_view.V0B.IsSignal()) ||
                                            (!a_view.V0A.IsSignal() && a_view.V0B.IsSignal()) ||
                                            (a_view.V0A.IsSignal() && a_view.V0B.IsSignal() && a_view.V0A.ReactionID() != a_view.V0B.ReactionID()));
}

// Channel D //

[[nodiscard]] inline double AsChannelD_AfterPx(const View::MC::ChannelD& d_view) {
    return d_view.IsValid() ? d_view.V0.Px() + d_view.Kaon.Px() : Const::DummyDouble;
}
[[nodiscard]] inline double AsChannelD_AfterPy(const View::MC::ChannelD& d_view) {
    return d_view.IsValid() ? d_view.V0.Py() + d_view.Kaon.Py() : Const::DummyDouble;
}
[[nodiscard]] inline double AsChannelD_AfterPz(const View::MC::ChannelD& d_view) {
    return d_view.IsValid() ? d_view.V0.Pz() + d_view.Kaon.Pz() : Const::DummyDouble;
}
[[nodiscard]] inline double AsChannelD_AfterE(const View::MC::ChannelD& d_view) {
    return d_view.IsValid() ? d_view.V0.Energy() + d_view.Kaon.Energy() : Const::DummyDouble;
}

[[nodiscard]] inline bool AsChannelD_IsSignal(const View::MC::ChannelD& d_view) {
    return d_view.IsValid() ? d_view.V0.IsSignal() && d_view.Kaon.IsSignal() : false;
}
[[nodiscard]] inline bool AsChannelD_IsHybrid(const View::MC::ChannelD& d_view) {
    // NOTE: don't check IsValid(), on purpose
    return !AsChannelD_IsSignal(d_view) &&
           (d_view.V0.IsHybrid() || (d_view.V0.IsSignal() && !d_view.Kaon.IsSignal()) || (!d_view.V0.IsSignal() && d_view.Kaon.IsSignal()) ||
            (d_view.V0.IsSignal() && d_view.Kaon.IsSignal() && d_view.V0.ReactionID() != d_view.Kaon.ReactionID()));
}

}  // namespace Tree2Secondaries::Truth::Sexaquark
