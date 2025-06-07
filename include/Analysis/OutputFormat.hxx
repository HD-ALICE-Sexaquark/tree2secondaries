#ifndef T2S_OUTPUT_FORMAT_HXX
#define T2S_OUTPUT_FORMAT_HXX

#include <cstddef>
#include <vector>

namespace Tree2Secondaries::OutputSOA {

struct Particle {
    std::vector<size_t>* Entry{nullptr};

    std::vector<float>* Xv{nullptr};
    std::vector<float>* Yv{nullptr};
    std::vector<float>* Zv{nullptr};
    std::vector<float>* Px{nullptr};
    std::vector<float>* Py{nullptr};
    std::vector<float>* Pz{nullptr};
    std::vector<float>* E{nullptr};

    std::vector<float>* SigmaX2{nullptr};
    std::vector<float>* SigmaYX{nullptr};
    std::vector<float>* SigmaY2{nullptr};
    std::vector<float>* SigmaZX{nullptr};
    std::vector<float>* SigmaZY{nullptr};
    std::vector<float>* SigmaZ2{nullptr};
    std::vector<float>* SigmaPxX{nullptr};
    std::vector<float>* SigmaPxY{nullptr};
    std::vector<float>* SigmaPxZ{nullptr};
    std::vector<float>* SigmaPx2{nullptr};
    std::vector<float>* SigmaPyX{nullptr};
    std::vector<float>* SigmaPyY{nullptr};
    std::vector<float>* SigmaPyZ{nullptr};
    std::vector<float>* SigmaPyPx{nullptr};
    std::vector<float>* SigmaPy2{nullptr};
    std::vector<float>* SigmaPzX{nullptr};
    std::vector<float>* SigmaPzY{nullptr};
    std::vector<float>* SigmaPzZ{nullptr};
    std::vector<float>* SigmaPzPx{nullptr};
    std::vector<float>* SigmaPzPy{nullptr};
    std::vector<float>* SigmaPz2{nullptr};
    std::vector<float>* SigmaEX{nullptr};
    std::vector<float>* SigmaEY{nullptr};
    std::vector<float>* SigmaEZ{nullptr};
    std::vector<float>* SigmaEPx{nullptr};
    std::vector<float>* SigmaEPy{nullptr};
    std::vector<float>* SigmaEPz{nullptr};
    std::vector<float>* SigmaE2{nullptr};

    void ClearParticle() {
        Entry->clear();

        Xv->clear();
        Yv->clear();
        Zv->clear();
        Px->clear();
        Py->clear();
        Pz->clear();
        E->clear();

        SigmaX2->clear();
        SigmaYX->clear();
        SigmaY2->clear();
        SigmaZX->clear();
        SigmaZY->clear();
        SigmaZ2->clear();
        SigmaPxX->clear();
        SigmaPxY->clear();
        SigmaPxZ->clear();
        SigmaPx2->clear();
        SigmaPyX->clear();
        SigmaPyY->clear();
        SigmaPyZ->clear();
        SigmaPyPx->clear();
        SigmaPy2->clear();
        SigmaPzX->clear();
        SigmaPzY->clear();
        SigmaPzZ->clear();
        SigmaPzPx->clear();
        SigmaPzPy->clear();
        SigmaPz2->clear();
        SigmaEX->clear();
        SigmaEY->clear();
        SigmaEZ->clear();
        SigmaEPx->clear();
        SigmaEPy->clear();
        SigmaEPz->clear();
        SigmaE2->clear();
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
        ClearParticle();
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
        ClearParticle();
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

}  // namespace Tree2Secondaries::OutputSOA

#endif  // T2S_OUTPUT_FORMAT_HXX
