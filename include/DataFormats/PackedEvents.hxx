#ifndef T2S_STRUCTS_PACKED_HXX
#define T2S_STRUCTS_PACKED_HXX

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
    void ClearCov() {
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
    std::vector<int>* Entry{nullptr};
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
        Sigma.ClearCov();
    }
};

struct alignas(32) Tracks : Particle {
    void Clear() { ClearParticle(); }
};

struct alignas(32) V0s : Particle {
    State Neg;
    State Pos;
    std::vector<float>* Neg_X_AtPCA{nullptr};
    std::vector<float>* Neg_Y_AtPCA{nullptr};
    std::vector<float>* Neg_Z_AtPCA{nullptr};
    std::vector<float>* Neg_Px_AtPCA{nullptr};
    std::vector<float>* Neg_Py_AtPCA{nullptr};
    std::vector<float>* Neg_Pz_AtPCA{nullptr};
    std::vector<float>* Pos_X_AtPCA{nullptr};
    std::vector<float>* Pos_Y_AtPCA{nullptr};
    std::vector<float>* Pos_Z_AtPCA{nullptr};
    std::vector<float>* Pos_Px_AtPCA{nullptr};
    std::vector<float>* Pos_Py_AtPCA{nullptr};
    std::vector<float>* Pos_Pz_AtPCA{nullptr};
    void Clear() {
        ClearParticle();
        Neg.ClearState();
        Pos.ClearState();
        Neg_X_AtPCA->clear();
        Neg_Y_AtPCA->clear();
        Neg_Z_AtPCA->clear();
        Neg_Px_AtPCA->clear();
        Neg_Py_AtPCA->clear();
        Neg_Pz_AtPCA->clear();
        Pos_X_AtPCA->clear();
        Pos_Y_AtPCA->clear();
        Pos_Z_AtPCA->clear();
        Pos_Px_AtPCA->clear();
        Pos_Py_AtPCA->clear();
        Pos_Pz_AtPCA->clear();
    }
};

struct alignas(32) MC_Tracks : State {
    std::vector<int>* Mother_Entry{nullptr};
    std::vector<int>* GrandMother_Entry{nullptr};
    std::vector<int>* PdgCode{nullptr};
    std::vector<int>* Mother_PdgCode{nullptr};
    std::vector<int>* GrandMother_PdgCode{nullptr};
    std::vector<int>* ReactionID{nullptr};
    std::vector<bool>* IsTrue{nullptr};
    std::vector<bool>* IsSignal{nullptr};
    std::vector<bool>* IsSecondary{nullptr};
    void Clear() {
        ClearState();
        Mother_Entry->clear();
        GrandMother_Entry->clear();
        PdgCode->clear();
        Mother_PdgCode->clear();
        GrandMother_PdgCode->clear();
        ReactionID->clear();
        IsTrue->clear();
        IsSignal->clear();
        IsSecondary->clear();
    }
};

struct alignas(32) MC_V0s : State {
    std::vector<float>* DecayX{nullptr};
    std::vector<float>* DecayY{nullptr};
    std::vector<float>* DecayZ{nullptr};

    std::vector<int>* PdgCode{nullptr};
    std::vector<int>* Mother_Entry{nullptr};
    std::vector<int>* Mother_PdgCode{nullptr};
    std::vector<bool>* IsTrue{nullptr};
    std::vector<bool>* IsSignal{nullptr};
    std::vector<bool>* IsSecondary{nullptr};
    std::vector<int>* ReactionID{nullptr};
    std::vector<bool>* IsHybrid{nullptr};

    // neg //
    std::vector<int>* Neg_Entry{nullptr};
    std::vector<float>* Neg_Px{nullptr};
    std::vector<float>* Neg_Py{nullptr};
    std::vector<float>* Neg_Pz{nullptr};
    std::vector<int>* Neg_PdgCode{nullptr};
    std::vector<bool>* Neg_IsTrue{nullptr};
    std::vector<bool>* Neg_IsSignal{nullptr};
    std::vector<bool>* Neg_IsSecondary{nullptr};
    std::vector<int>* Neg_ReactionID{nullptr};

    // pos //
    std::vector<int>* Pos_Entry{nullptr};
    std::vector<float>* Pos_Px{nullptr};
    std::vector<float>* Pos_Py{nullptr};
    std::vector<float>* Pos_Pz{nullptr};
    std::vector<int>* Pos_PdgCode{nullptr};
    std::vector<bool>* Pos_IsTrue{nullptr};
    std::vector<bool>* Pos_IsSignal{nullptr};
    std::vector<bool>* Pos_IsSecondary{nullptr};
    std::vector<int>* Pos_ReactionID{nullptr};

    void Clear() {
        ClearState();
        DecayX->clear();
        DecayY->clear();
        DecayZ->clear();

        PdgCode->clear();
        Mother_Entry->clear();
        Mother_PdgCode->clear();
        IsTrue->clear();
        IsSignal->clear();
        IsSecondary->clear();
        ReactionID->clear();
        IsHybrid->clear();

        // neg //
        Neg_Entry->clear();
        Neg_Px->clear();
        Neg_Py->clear();
        Neg_Pz->clear();
        Neg_PdgCode->clear();
        Neg_IsTrue->clear();
        Neg_IsSignal->clear();
        Neg_IsSecondary->clear();
        Neg_ReactionID->clear();

        // pos //
        Pos_Entry->clear();
        Pos_Px->clear();
        Pos_Py->clear();
        Pos_Pz->clear();
        Pos_PdgCode->clear();
        Pos_IsTrue->clear();
        Pos_IsSignal->clear();
        Pos_IsSecondary->clear();
        Pos_ReactionID->clear();
    }
};

}  // namespace Tree2Secondaries::PackedEvents

#endif  // T2S_STRUCTS_PACKED_HXX
