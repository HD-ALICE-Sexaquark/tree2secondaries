#pragma once

#include <cmath>
#include <cstddef>

#include "Math/Constants.hxx"
#include "View/ViewVectorInjected.hxx"
#include "View/ViewVectorMcTracks.hxx"
#include "View/ViewVectorMcV0s.hxx"

namespace Tree2Secondaries::MonteCarlo {

struct ChannelD : View::VecInjected {
   public:
    // ## Constructor ## //

    ChannelD() = delete;
    static std::optional<ChannelD> Create(const Schema::VecInjected* df, const View::VecMcV0s& mc_v0, const View::VecMcTracks& mc_kaon) {
        if (mc_kaon.ReactionID() < 0 || mc_v0.ReactionID() < 0) return std::nullopt;
        if (mc_kaon.ReactionID() != mc_v0.ReactionID()) return std::nullopt;
        auto injected_idx = static_cast<std::size_t>(mc_kaon.ReactionID() - Const::ReactionID_Offset);
        return ChannelD(df, injected_idx, mc_v0, mc_kaon);
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

    ChannelD(const Schema::VecInjected* df, std::size_t protected_idx, const View::VecMcV0s& mc_v0, const View::VecMcTracks& mc_kaon)
        : View::VecInjected{df, protected_idx},
          after_x{mc_kaon.Px() + mc_v0.Px()},
          after_y{mc_kaon.Py() + mc_v0.Py()},
          after_z{mc_kaon.Pz() + mc_v0.Pz()},
          after_energy{mc_kaon.Energy() + mc_v0.Energy()},
          is_signal{true},
          is_hybrid{(mc_kaon.IsSignal() != mc_v0.IsSignal()) || mc_v0.IsHybrid()} {}
};

}  // namespace Tree2Secondaries::MonteCarlo
