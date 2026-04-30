#pragma once

#include <optional>

#include "Math/Constants.hxx"
#include "MonteCarlo/MonteCarloParticle.hxx"

namespace Tree2Secondaries::MonteCarlo {

struct V0 : Particle {
   public:
    // ## Constructor ## //

    V0() = delete;
    static std::optional<V0> Create(const std::optional<Particle>& mc_neg, const std::optional<Particle>& mc_pos,
                                    const ReactionChannels::Definition& r_channel, int hypothesis_pdg = Const::DummyInt) {
        if (!mc_neg || !mc_pos) return std::nullopt;
        if (mc_neg->Mother_McEntry() < 0 || mc_pos->Mother_McEntry() < 0) return std::nullopt;
        if (mc_neg->Mother_McEntry() != mc_pos->Mother_McEntry()) return std::nullopt;
        return V0(mc_neg.value(), mc_pos.value(), r_channel, hypothesis_pdg);
    }

    // ## Getters ## //

    [[nodiscard]] bool IsHybrid() const { return is_hybrid; }

    // ## Member Variables ## //

    bool is_hybrid{false};

   protected:
    // ## Constructor ## //

    V0(const Particle& mc_neg, const Particle& mc_pos, const ReactionChannels::Definition& r_channel, int hypothesis_pdg)
        : Particle{mc_neg.source,
                   static_cast<std::size_t>(mc_neg.Mother_McEntry()),  //
                   r_channel,                                          //
                   hypothesis_pdg},
          is_hybrid{mc_neg.IsSignal() != mc_pos.IsSignal()} {}
};

}  // namespace Tree2Secondaries::MonteCarlo
