#ifndef T2S_SECONDARY_TRUE_HXX
#define T2S_SECONDARY_TRUE_HXX

#include <cstdlib>

#include "Secondary/Particle.hxx"

namespace Tree2Secondaries {

// True information from MC.
class True final : public Particle {
   public:
    True(const True &) = delete;
    True(True &&) = default;
    True &operator=(const True &) = delete;
    True &operator=(True &&) = default;
    ~True() final = default;

    True() : fMC_Index{}, fMother_McEntry{-1}, fPdgCode{0}, fReactionID{0}, fIsSignal{false} {}
    True(size_t mc_idx, long long mother_mc_idx, int pdg_code, int reaction_id, bool is_signal)
        : fMC_Index{mc_idx},  //
          fMother_McEntry{mother_mc_idx},
          fPdgCode{pdg_code},
          fReactionID{reaction_id},
          fIsSignal{is_signal} {}

    [[nodiscard]] size_t Index() const override { return fMC_Index; }
    [[nodiscard]] long long MotherMcEntry() const { return fMother_McEntry; }
    [[nodiscard]] int PdgCode() const { return fPdgCode; }
    [[nodiscard]] int ReactionID() const { return fReactionID; }
    [[nodiscard]] bool IsSignal() const { return fIsSignal; }
    // [[nodiscard]] bool IsPrimary() const { return fIsPrimary; } // PENDING

   protected:
    size_t fMC_Index;
    long long fMother_McEntry;
    int fPdgCode;
    // bool fIsPrimary; // PENDING
    int fReactionID;
    bool fIsSignal;
};

}  // namespace Tree2Secondaries

#endif  // T2S_SECONDARY_TRUE_HXX
