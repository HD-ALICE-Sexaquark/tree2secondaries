#ifndef T2S_STRUCTS_PACKED_HXX
#define T2S_STRUCTS_PACKED_HXX

#include <cstddef>
#include <vector>

namespace Tree2Secondaries::PackedEvents {

struct Particles {
    std::vector<size_t>* Entry{nullptr};

    std::vector<float>* Xv{nullptr};
    std::vector<float>* Yv{nullptr};
    std::vector<float>* Zv{nullptr};
    std::vector<float>* Px{nullptr};
    std::vector<float>* Py{nullptr};
    std::vector<float>* Pz{nullptr};
    std::vector<float>* E{nullptr};

    std::vector<float>* SigmaX2{nullptr};
    std::vector<float>* SigmaXY{nullptr};
    std::vector<float>* SigmaY2{nullptr};
    std::vector<float>* SigmaXZ{nullptr};
    std::vector<float>* SigmaYZ{nullptr};
    std::vector<float>* SigmaZ2{nullptr};
    std::vector<float>* SigmaXPx{nullptr};
    std::vector<float>* SigmaYPx{nullptr};
    std::vector<float>* SigmaZPx{nullptr};
    std::vector<float>* SigmaPx2{nullptr};
    std::vector<float>* SigmaXPy{nullptr};
    std::vector<float>* SigmaYPy{nullptr};
    std::vector<float>* SigmaZPy{nullptr};
    std::vector<float>* SigmaPxPy{nullptr};
    std::vector<float>* SigmaPy2{nullptr};
    std::vector<float>* SigmaXPz{nullptr};
    std::vector<float>* SigmaYPz{nullptr};
    std::vector<float>* SigmaZPz{nullptr};
    std::vector<float>* SigmaPxPz{nullptr};
    std::vector<float>* SigmaPyPz{nullptr};
    std::vector<float>* SigmaPz2{nullptr};
    std::vector<float>* SigmaXE{nullptr};
    std::vector<float>* SigmaYE{nullptr};
    std::vector<float>* SigmaZE{nullptr};
    std::vector<float>* SigmaPxE{nullptr};
    std::vector<float>* SigmaPyE{nullptr};
    std::vector<float>* SigmaPzE{nullptr};
    std::vector<float>* SigmaE2{nullptr};

    void ClearParticles() {
        Entry->clear();

        Xv->clear();
        Yv->clear();
        Zv->clear();
        Px->clear();
        Py->clear();
        Pz->clear();
        E->clear();

        SigmaX2->clear();
        SigmaXY->clear();
        SigmaY2->clear();
        SigmaXZ->clear();
        SigmaYZ->clear();
        SigmaZ2->clear();
        SigmaXPx->clear();
        SigmaYPx->clear();
        SigmaZPx->clear();
        SigmaPx2->clear();
        SigmaXPy->clear();
        SigmaYPy->clear();
        SigmaZPy->clear();
        SigmaPxPy->clear();
        SigmaPy2->clear();
        SigmaXPz->clear();
        SigmaYPz->clear();
        SigmaZPz->clear();
        SigmaPxPz->clear();
        SigmaPyPz->clear();
        SigmaPz2->clear();
        SigmaXE->clear();
        SigmaYE->clear();
        SigmaZE->clear();
        SigmaPxE->clear();
        SigmaPyE->clear();
        SigmaPzE->clear();
        SigmaE2->clear();
    }
};

struct MC : public Particles {
    std::vector<int>* PdgCode{nullptr};
    std::vector<int>* MotherEntry{nullptr};
    std::vector<bool>* IsSignal{nullptr};
    std::vector<bool>* IsPrimary{nullptr};
    std::vector<bool>* IsSecFromMat{nullptr};
    std::vector<bool>* IsSecFromWeak{nullptr};
    std::vector<int>* ReactionID{nullptr};

    void ClearMC() {
        ClearParticles();
        PdgCode->clear();
        MotherEntry->clear();
        IsSignal->clear();
        IsPrimary->clear();
        IsSecFromMat->clear();
        IsSecFromWeak->clear();
        ReactionID->clear();
    }
};

struct Tracks : public Particles {
    MC True;
    void Clear(bool is_mc = false) {
        ClearParticles();
        if (is_mc) True.ClearMC();
    }
};

struct V0s : public Particles {
    std::vector<size_t>* Neg_Entry{nullptr};
    std::vector<float>* Neg_Xv{nullptr};
    std::vector<float>* Neg_Yv{nullptr};
    std::vector<float>* Neg_Zv{nullptr};
    std::vector<float>* Neg_Px{nullptr};
    std::vector<float>* Neg_Py{nullptr};
    std::vector<float>* Neg_Pz{nullptr};
    std::vector<size_t>* Pos_Entry{nullptr};
    std::vector<float>* Pos_Xv{nullptr};
    std::vector<float>* Pos_Yv{nullptr};
    std::vector<float>* Pos_Zv{nullptr};
    std::vector<float>* Pos_Px{nullptr};
    std::vector<float>* Pos_Py{nullptr};
    std::vector<float>* Pos_Pz{nullptr};
    MC True;
    void Clear(bool is_mc = false) {
        ClearParticles();
        Neg_Entry->clear();
        Neg_Xv->clear();
        Neg_Yv->clear();
        Neg_Zv->clear();
        Neg_Px->clear();
        Neg_Py->clear();
        Neg_Pz->clear();
        Pos_Entry->clear();
        Pos_Xv->clear();
        Pos_Yv->clear();
        Pos_Zv->clear();
        Pos_Px->clear();
        Pos_Py->clear();
        Pos_Pz->clear();
        if (is_mc) True.ClearMC();
    }
};

}  // namespace Tree2Secondaries::PackedEvents

#endif  // T2S_STRUCTS_PACKED_HXX
