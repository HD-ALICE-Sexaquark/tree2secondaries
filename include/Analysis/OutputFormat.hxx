#ifndef T2S_OUTPUT_FORMAT_HXX
#define T2S_OUTPUT_FORMAT_HXX

#include <vector>

namespace Tree2Secondaries::OutputSOA {

struct Particle {
    std::vector<int>* Entry{nullptr};

    std::vector<float>* Xv{nullptr};
    std::vector<float>* Yv{nullptr};
    std::vector<float>* Zv{nullptr};
    std::vector<float>* Px{nullptr};
    std::vector<float>* Py{nullptr};
    std::vector<float>* Pz{nullptr};
    std::vector<float>* E{nullptr};

    void ClearParticle() {
        Entry->clear();
        Xv->clear();
        Yv->clear();
        Zv->clear();
        Px->clear();
        Py->clear();
        Pz->clear();
        E->clear();
    }
};

struct MC : public Particle {
    std::vector<int>* PdgCode{nullptr};
    std::vector<int>* MotherEntry{nullptr};
    std::vector<bool>* IsSignal{nullptr};
    std::vector<bool>* IsPrimary{nullptr};
    std::vector<bool>* IsSecFromMat{nullptr};
    std::vector<bool>* IsSecFromWeak{nullptr};
    std::vector<int>* ReactionID{nullptr};

    void ClearMC() {
        PdgCode->clear();
        MotherEntry->clear();
        IsSignal->clear();
        IsPrimary->clear();
        IsSecFromMat->clear();
        IsSecFromWeak->clear();
        ReactionID->clear();
    }
};

struct Tracks : public Particle {
    MC True;
    void Clear(bool is_mc = false) {
        ClearParticle();
        if (is_mc) True.ClearMC();
    }
};

struct V0s : public Particle {
    Tracks Neg;
    Tracks Pos;
    MC True;
    void Clear(bool is_mc = false) {
        ClearParticle();
        Neg.ClearParticle();
        Pos.ClearParticle();
        if (is_mc) True.ClearMC();
    }
};

}  // namespace Tree2Secondaries::OutputSOA

#endif  // T2S_OUTPUT_FORMAT_HXX
