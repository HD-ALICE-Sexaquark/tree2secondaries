#ifndef T2S_TRUTH_HANDLER_HXX
#define T2S_TRUTH_HANDLER_HXX

#include <vector>

#include "Analysis/InputFormat.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::Helper {

// Contains all logical operations that involve true information.
class TruthHandler {
   public:
    ~TruthHandler() = default;

    void InitMap(int n) { fMcEntry.resize(n, -1); };
    void Link(int track_entry, int mc_entry) { fMcEntry[track_entry] = mc_entry; }
    void Clear() { fMcEntry.clear(); };

    [[nodiscard]] int McEntry(int track_entry) const { return fMcEntry[track_entry]; }

    [[nodiscard]] float Xv(int mc_entry) const { return fInput_MC.Xv->at(mc_entry); }
    [[nodiscard]] float Yv(int mc_entry) const { return fInput_MC.Yv->at(mc_entry); }
    [[nodiscard]] float Zv(int mc_entry) const { return fInput_MC.Zv->at(mc_entry); }
    [[nodiscard]] float Px(int mc_entry) const { return fInput_MC.Px->at(mc_entry); }
    [[nodiscard]] float Py(int mc_entry) const { return fInput_MC.Py->at(mc_entry); }
    [[nodiscard]] float Pz(int mc_entry) const { return fInput_MC.Pz->at(mc_entry); }
    [[nodiscard]] float E(int mc_entry) const { return fInput_MC.E->at(mc_entry); }

    [[nodiscard]] int PdgCode(int mc_entry) const { return fInput_MC.PdgCode->at(mc_entry); }
    [[nodiscard]] int MotherEntry(int mc_entry) const { return fInput_MC.MotherEntry->at(mc_entry); }
    [[nodiscard]] int Status(int mc_entry) const { return fInput_MC.Status->at(mc_entry); }
    [[nodiscard]] int Generator(int mc_entry) const { return fInput_MC.Generator->at(mc_entry); }
    [[nodiscard]] bool IsPrimary(int mc_entry) const { return fInput_MC.IsPrimary->at(mc_entry); }
    [[nodiscard]] bool IsSecFromMat(int mc_entry) const { return fInput_MC.IsSecFromMat->at(mc_entry); }
    [[nodiscard]] bool IsSecFromWeak(int mc_entry) const { return fInput_MC.IsSecFromWeak->at(mc_entry); }

    // derived operations //

    [[nodiscard]] bool IsSignal(int mc_entry) const { return Generator(mc_entry) == 2; }
    [[nodiscard]] bool IsSignal(int v0_mc_entry, int hypothesis_pid) const {
        return Generator(v0_mc_entry) == 2 && PdgCode(v0_mc_entry) == hypothesis_pid;
    }
    [[nodiscard]] int ReactionID(int mc_entry) const {
        if (IsSignal(mc_entry)) {
            int mother_entry{MotherEntry(mc_entry)};
            return mother_entry > Const::DummyInt ? Status(mother_entry) : Status(mc_entry);
        }
        return Const::DummyInt;
    }
    [[nodiscard]] int SameMother(int neg_mc_entry, int pos_mc_entry) const {
        const auto mother_neg_mc_entry{MotherEntry(neg_mc_entry)};
        if (mother_neg_mc_entry == MotherEntry(pos_mc_entry)) return mother_neg_mc_entry;
        return Const::DummyInt;
    }

    InputSOA::MC fInput_MC;

   private:
    std::vector<int> fMcEntry;
};

}  // namespace Tree2Secondaries::Helper

#endif
