#pragma once

#include <cmath>
#include <cstddef>

#include "Math/Constants.hxx"
#include "View/ViewVectorInjected.hxx"
#include "View/ViewVectorMcV0s.hxx"

namespace Tree2Secondaries::MonteCarlo {

struct ChannelA : View::VecInjected {
   public:
    // ## Constructor ## //

    ChannelA() = delete;
    static std::optional<ChannelA> Create(const Schema::VecInjected* df, const View::VecMcV0s& mc_v0a, const View::VecMcV0s& mc_v0b) {
        if (mc_v0a.ReactionID() < 0 || mc_v0b.ReactionID() < 0) return std::nullopt;
        if (mc_v0a.ReactionID() != mc_v0b.ReactionID()) return std::nullopt;
        auto injected_idx = static_cast<std::size_t>(mc_v0a.ReactionID() - Const::ReactionID_Offset);
        return ChannelA(df, injected_idx, mc_v0a, mc_v0b);
    }

    // ## Getters ## //

    [[nodiscard]] float After_Px() const { return after_x; }
    [[nodiscard]] float After_Py() const { return after_y; }
    [[nodiscard]] float After_Pz() const { return after_z; }
    [[nodiscard]] float After_Energy() const { return after_energy; }
    [[nodiscard]] bool IsSignal() const { return is_signal; }
    [[nodiscard]] bool IsHybrid() const { return is_hybrid; }

    // ## Member Variables ## //

    float after_x{Const::DummyFloat};
    float after_y{Const::DummyFloat};
    float after_z{Const::DummyFloat};
    float after_energy{Const::DummyFloat};
    bool is_signal{false};
    bool is_hybrid{false};

   protected:
    // ## Constructor ## //

    ChannelA(const Schema::VecInjected* df, std::size_t protected_idx, const View::VecMcV0s& mc_v0a, const View::VecMcV0s& mc_v0b)
        : View::VecInjected{df, protected_idx},
          after_x{mc_v0a.Px() + mc_v0b.Px()},
          after_y{mc_v0a.Py() + mc_v0b.Py()},
          after_z{mc_v0a.Pz() + mc_v0b.Pz()},
          after_energy{mc_v0a.Energy() + mc_v0b.Energy()},
          is_signal{true},
          is_hybrid{(mc_v0a.IsSignal() != mc_v0b.IsSignal()) || mc_v0a.IsHybrid() || mc_v0b.IsHybrid()} {}
};

}  // namespace Tree2Secondaries::MonteCarlo
