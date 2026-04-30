#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>

#include "App/DB_Particles.hxx"
#include "App/DB_ReactionChannels.hxx"
#include "Math/Constants.hxx"
#include "View/ViewVectorMcParticles.hxx"

namespace Tree2Secondaries::MonteCarlo {

// Extended version of `View::VecMcParticles` with three main features:
// (1) caches some ascendants and descendants info;
// (2) caches the classification info of the particle, whether it's true, signal, etc.;
// (3) provides a factory method `Create` which returns `std::nullopt` when given an invalid index.
struct Particle : View::VecMcParticles {
   public:
    // ## Constructors ## //

    Particle() = delete;

    template <typename T>
    static std::optional<Particle> Create(const Schema::VecMcParticles* df, T unprotected_idx, const ReactionChannels::Definition& r_channel,
                                          int hypothesis_pdg) {
        if (unprotected_idx < 0) return std::nullopt;
        return Particle(df, static_cast<std::size_t>(unprotected_idx), r_channel, hypothesis_pdg);
    }

    // ## Getters ## //

    [[nodiscard]] float Decay_X() const { return decay_x; }
    [[nodiscard]] float Decay_Y() const { return decay_y; }
    [[nodiscard]] float Decay_Z() const { return decay_z; }
    [[nodiscard]] int McEntry() const { return mc_entry; }
    [[nodiscard]] int ReactionID() const { return reaction_id; }
    [[nodiscard]] int Mother_PdgCode() const { return mother_pdg_code; }
    [[nodiscard]] int GrandMother_McEntry() const { return gm_mc_entry; }
    [[nodiscard]] int GrandMother_PdgCode() const { return gm_pdg_code; }
    [[nodiscard]] bool IsTrue() const { return is_true; }
    [[nodiscard]] bool IsSecondary() const { return IsSecFromMat() || IsSecFromWeak() || IsSignal(); }
    [[nodiscard]] bool IsFirstGenSignal() const { return is_gen1_signal; }
    [[nodiscard]] bool IsSecondGenSignal() const { return is_gen2_signal; }
    [[nodiscard]] bool IsSignal() const { return IsFirstGenSignal() || IsSecondGenSignal(); }

    // ## Member Variables ## //

    float decay_x{Const::DummyFloat};
    float decay_y{Const::DummyFloat};
    float decay_z{Const::DummyFloat};
    int mother_pdg_code{Const::DummyInt};
    unsigned int mother_status{0};
    int gm_mc_entry{Const::DummyInt};
    int gm_pdg_code{Const::DummyInt};
    int mc_entry{Const::DummyInt};
    int reaction_id{Const::DummyInt};
    bool mother_is_primary{false};
    bool is_true{false};
    bool is_gen1_signal{false};
    bool is_gen2_signal{false};

   protected:
    // ## Constructor ## //

    Particle(const Schema::VecMcParticles* df, std::size_t protected_idx, const ReactionChannels::Definition& r_channel, int hypothesis_pdg)
        : View::VecMcParticles(df, protected_idx),  //
          mc_entry{static_cast<int>(protected_idx)} {
        CacheAscendantsInfo();
        CacheDescendantsInfo();
        Classify(r_channel, hypothesis_pdg);
    }

    // ## Operations ## //

    // Fill `mother_{pdg_code/status/is_primary}` and `gm_{mc_entry/pgd_code}`.
    void CacheAscendantsInfo() {
        // -- mother
        if (Mother_McEntry() < 0) return;
        View::VecMcParticles mother(source, static_cast<std::size_t>(Mother_McEntry()));
        mother_pdg_code = mother.PdgCode();
        mother_status = mother.Status();
        mother_is_primary = mother.IsLogicalPrimary();
        // -- grandmother
        gm_mc_entry = mother.Mother_McEntry();
        if (gm_mc_entry < 0) return;
        View::VecMcParticles gm(source, static_cast<std::size_t>(gm_mc_entry));
        gm_pdg_code = gm.PdgCode();
    }

    // Fill `decay_{x/y/z}`.
    void CacheDescendantsInfo() {
        if (FirstDau_McEntry() < 0) return;
        View::VecMcParticles dau(source, static_cast<std::size_t>(FirstDau_McEntry()));
        decay_x = dau.Origin_X();
        decay_y = dau.Origin_Y();
        decay_z = dau.Origin_Z();
    }

    // Fill `is_true`, `is_signal`, `reaction_id`.
    // NOTE: only Gen1 and Gen2 products are considered.
    void Classify(const ReactionChannels::Definition& r_channel, int hypothesis_pdg) {
        // (1) particle should be true
        // NOTE: it's always true when there's no hypothesis
        int this_pdg_code = PdgCode();
        is_true = hypothesis_pdg == Const::DummyInt || hypothesis_pdg == this_pdg_code;
        if (!is_true) return;
        // (2) should come from the anti-sexaquark reaction generator
        if (Generator() != Const::EGenerator::kInjectedAntiSexaquarkReaction) return;
        unsigned int mc_status = Status();
        // -- in case it could be first-gen product
        //    (a.1) pdg is found in reaction products?
        if (std::ranges::find(r_channel.products_pdg, this_pdg_code) != r_channel.products_pdg.end()) {
            //    (a.2) has to be marked as primary
            //    (a.3) mc status has to be [600,620[
            is_gen1_signal = IsLogicalPrimary() &&  //
                             Const::ReactionID_Offset <= mc_status && mc_status < Const::ReactionID_Offset + Const::NInjectedPerEvent;
        }
        if (is_gen1_signal) {
            reaction_id = static_cast<int>(mc_status);
            return;
        }
        // -- in case it could be a second-gen product
        //    (b.1) needs to have a mother (implicit)
        //    (b.2) NOTE: don't check on mother's Generator, because it should have the same as daughter's, by construction
        //    (b.3) mother should be marked as primary
        if (!mother_is_primary) return;
        //    (b.4) mother should have a valid reaction id
        if (mother_status < Const::ReactionID_Offset || mother_status >= Const::ReactionID_Offset + Const::NInjectedPerEvent) return;
        //    (b.5) mother should be in first-gen products
        if (std::ranges::find(r_channel.products_pdg, mother_pdg_code) == r_channel.products_pdg.end()) return;
        //    (b.6) mother should be relevant and have detectable-charged decay modes
        auto mother = Particles::FindParticle(mother_pdg_code);
        if (!mother || mother->daughters_pdg.empty()) return;
        //    (b.7) current particle should appear in mother's decay mode
        if (std::ranges::find(mother->daughters_pdg, PdgCode()) == mother->daughters_pdg.end()) return;
        is_gen2_signal = true;
        reaction_id = static_cast<int>(mother_status);
    }
};

}  // namespace Tree2Secondaries::MonteCarlo
