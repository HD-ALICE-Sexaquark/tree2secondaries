#pragma once

#include <cmath>

#include "View/MC/ViewMcInjected.hxx"
#include "View/MC/ViewMcPackedTrack.hxx"
#include "View/MC/ViewMcPackedV0.hxx"

namespace Tree2Secondaries::Truth::Sexaquark {

[[nodiscard]] inline int ReactionID(const View::MC::Injected& sexa) { return sexa.ReactionID(); }

[[nodiscard]] inline double Energy(const View::MC::Injected& sexa, double mass) {
    auto px = static_cast<double>(sexa.Px());
    auto py = static_cast<double>(sexa.Py());
    auto pz = static_cast<double>(sexa.Pz());
    double squared_energy = px * px + py * py + pz * pz + mass * mass;
    return std::sqrt(squared_energy);
}

[[nodiscard]] inline double NucleonEnergy(const View::MC::Injected& sexa, double mass) {
    auto px = static_cast<double>(sexa.Nucleon_Px());
    auto py = static_cast<double>(sexa.Nucleon_Py());
    auto pz = static_cast<double>(sexa.Nucleon_Pz());
    double squared_energy = px * px + py * py + pz * pz + mass * mass;
    return std::sqrt(squared_energy);
}

// Channel A //

[[nodiscard]] inline double AfterPx(const View::MC::ChannelA& sexa) {
    return static_cast<double>(sexa.V0A.Px()) + static_cast<double>(sexa.V0B.Px());
}
[[nodiscard]] inline double AfterPy(const View::MC::ChannelA& sexa) {
    return static_cast<double>(sexa.V0A.Py()) + static_cast<double>(sexa.V0B.Py());
}
[[nodiscard]] inline double AfterPz(const View::MC::ChannelA& sexa) {
    return static_cast<double>(sexa.V0A.Pz()) + static_cast<double>(sexa.V0B.Pz());
}
[[nodiscard]] inline double AfterE(const View::MC::ChannelA& sexa) {
    return static_cast<double>(sexa.V0A.Energy()) + static_cast<double>(sexa.V0B.Energy());
}

[[nodiscard]] inline bool IsSignal(const View::MC::ChannelA& sexa) {
    return sexa.V0A.IsSignal() && sexa.V0B.IsSignal() && sexa.V0A.ReactionID() == sexa.V0B.ReactionID();
}
[[nodiscard]] inline bool IsHybrid(const View::MC::ChannelA& sexa) {
    return !IsSignal(sexa) && (sexa.V0A.IsHybrid() || sexa.V0B.IsHybrid() ||     //
                               (sexa.V0A.IsSignal() && !sexa.V0B.IsSignal()) ||  //
                               (!sexa.V0A.IsSignal() && sexa.V0B.IsSignal()) ||  //
                               (sexa.V0A.IsSignal() && sexa.V0B.IsSignal() && sexa.V0A.ReactionID() != sexa.V0B.ReactionID()));
}

// Channel D //

[[nodiscard]] inline double AfterPx(const View::MC::ChannelD& sexa) {
    return static_cast<double>(sexa.V0.Px()) + static_cast<double>(sexa.Kaon.Px());
}
[[nodiscard]] inline double AfterPy(const View::MC::ChannelD& sexa) {
    return static_cast<double>(sexa.V0.Py()) + static_cast<double>(sexa.Kaon.Py());
}
[[nodiscard]] inline double AfterPz(const View::MC::ChannelD& sexa) {
    return static_cast<double>(sexa.V0.Pz()) + static_cast<double>(sexa.Kaon.Pz());
}
[[nodiscard]] inline double AfterE(const View::MC::ChannelD& sexa) {
    return static_cast<double>(sexa.V0.Energy()) + static_cast<double>(sexa.Kaon.Energy());
}

[[nodiscard]] inline bool IsSignal(const View::MC::ChannelD& sexa) {
    return sexa.V0.IsSignal() && sexa.Kaon.IsSignal() && sexa.V0.ReactionID() == sexa.Kaon.ReactionID();
}
[[nodiscard]] inline bool IsHybrid(const View::MC::ChannelD& sexa) {
    return !IsSignal(sexa) && (sexa.V0.IsHybrid() ||                             //
                               (sexa.V0.IsSignal() && !sexa.Kaon.IsSignal()) ||  //
                               (!sexa.V0.IsSignal() && sexa.Kaon.IsSignal()) ||  //
                               (sexa.V0.IsSignal() && sexa.Kaon.IsSignal() && sexa.V0.ReactionID() != sexa.Kaon.ReactionID()));
}

}  // namespace Tree2Secondaries::Truth::Sexaquark
