#ifndef T2S_STRUCTS_PACKED_HXX
#define T2S_STRUCTS_PACKED_HXX

#include <cstddef>
#include <vector>

namespace Tree2Secondaries::PackedEvents {

struct alignas(32) Cov {
    std::vector<float>* X2{nullptr};
    std::vector<float>* XY{nullptr};
    std::vector<float>* Y2{nullptr};
    std::vector<float>* XZ{nullptr};
    std::vector<float>* YZ{nullptr};
    std::vector<float>* Z2{nullptr};
    std::vector<float>* XPx{nullptr};
    std::vector<float>* YPx{nullptr};
    std::vector<float>* ZPx{nullptr};
    std::vector<float>* Px2{nullptr};
    std::vector<float>* XPy{nullptr};
    std::vector<float>* YPy{nullptr};
    std::vector<float>* ZPy{nullptr};
    std::vector<float>* PxPy{nullptr};
    std::vector<float>* Py2{nullptr};
    std::vector<float>* XPz{nullptr};
    std::vector<float>* YPz{nullptr};
    std::vector<float>* ZPz{nullptr};
    std::vector<float>* PxPz{nullptr};
    std::vector<float>* PyPz{nullptr};
    std::vector<float>* Pz2{nullptr};
    std::vector<float>* XE{nullptr};
    std::vector<float>* YE{nullptr};
    std::vector<float>* ZE{nullptr};
    std::vector<float>* PxE{nullptr};
    std::vector<float>* PyE{nullptr};
    std::vector<float>* PzE{nullptr};
    std::vector<float>* E2{nullptr};
    void Clear() {
        X2->clear();
        XY->clear();
        Y2->clear();
        XZ->clear();
        YZ->clear();
        Z2->clear();
        XPx->clear();
        YPx->clear();
        ZPx->clear();
        Px2->clear();
        XPy->clear();
        YPy->clear();
        ZPy->clear();
        PxPy->clear();
        Py2->clear();
        XPz->clear();
        YPz->clear();
        ZPz->clear();
        PxPz->clear();
        PyPz->clear();
        Pz2->clear();
        XE->clear();
        YE->clear();
        ZE->clear();
        PxE->clear();
        PyE->clear();
        PzE->clear();
        E2->clear();
    }
};

struct alignas(32) State {
    std::vector<size_t>* Entry{nullptr};
    std::vector<float>* X{nullptr};
    std::vector<float>* Y{nullptr};
    std::vector<float>* Z{nullptr};
    std::vector<float>* Px{nullptr};
    std::vector<float>* Py{nullptr};
    std::vector<float>* Pz{nullptr};
    std::vector<float>* E{nullptr};
    void ClearState() {
        Entry->clear();
        X->clear();
        Y->clear();
        Z->clear();
        Px->clear();
        Py->clear();
        Pz->clear();
        E->clear();
    }
};

struct alignas(32) Particle : State {
    Cov Sigma;
    void ClearParticle() {
        ClearState();
        Sigma.Clear();
    }
};

struct alignas(32) MC : public State {
    std::vector<int>* PdgCode{nullptr};
    std::vector<int>* MotherEntry{nullptr};
    std::vector<bool>* IsSignal{nullptr};
    std::vector<bool>* IsPrimary{nullptr};
    std::vector<bool>* IsSecFromMat{nullptr};
    std::vector<bool>* IsSecFromWeak{nullptr};
    std::vector<int>* ReactionID{nullptr};
    void ClearMC() {
        ClearState();
        PdgCode->clear();
        MotherEntry->clear();
        IsSignal->clear();
        IsPrimary->clear();
        IsSecFromMat->clear();
        IsSecFromWeak->clear();
        ReactionID->clear();
    }
};

struct alignas(32) Tracks : public Particle {
    MC True;
    void Clear(bool is_mc = false) {
        ClearParticle();
        Sigma.Clear();
        if (is_mc) True.ClearMC();
    }
};

struct alignas(32) V0s : public Particle {
    MC True;
    State Neg;
    State Pos;
    void Clear(bool is_mc = false) {
        ClearParticle();
        Neg.ClearState();
        Pos.ClearState();
        if (is_mc) True.ClearMC();
    }
};

}  // namespace Tree2Secondaries::PackedEvents

#endif  // T2S_STRUCTalignas(32)S_PACKED_HXX
