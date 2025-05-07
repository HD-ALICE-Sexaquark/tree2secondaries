#ifndef T2S_OUTPUT_FORMAT_HXX
#define T2S_OUTPUT_FORMAT_HXX

#include <vector>

namespace Tree2Secondaries::Output {

struct V0s {
    std::vector<int>* Index{nullptr};
    std::vector<int>* PID{nullptr};  // PID Hypothesis
    std::vector<int>* Neg_Entry{nullptr};
    std::vector<int>* Pos_Entry{nullptr};

    std::vector<float>* Xv{nullptr};
    std::vector<float>* Yv{nullptr};
    std::vector<float>* Zv{nullptr};
    std::vector<float>* Px{nullptr};
    std::vector<float>* Py{nullptr};
    std::vector<float>* Pz{nullptr};
    std::vector<float>* E{nullptr};

    std::vector<float>* Neg_Xv{nullptr};
    std::vector<float>* Neg_Yv{nullptr};
    std::vector<float>* Neg_Zv{nullptr};
    std::vector<float>* Neg_Px{nullptr};
    std::vector<float>* Neg_Py{nullptr};
    std::vector<float>* Neg_Pz{nullptr};

    std::vector<float>* Pos_Xv{nullptr};
    std::vector<float>* Pos_Yv{nullptr};
    std::vector<float>* Pos_Zv{nullptr};
    std::vector<float>* Pos_Px{nullptr};
    std::vector<float>* Pos_Py{nullptr};
    std::vector<float>* Pos_Pz{nullptr};

    /*
    // PENDING
    unsigned int reaction_id{};
    unsigned int index{};
    int pdg_code{};
    bool is_signal{};
    bool is_hybrid{};
    */

    void Clear() {
        Index->clear();
        PID->clear();  // PID Hypothesis
        Neg_Entry->clear();
        Pos_Entry->clear();

        Xv->clear();
        Yv->clear();
        Zv->clear();
        Px->clear();
        Py->clear();
        Pz->clear();
        E->clear();

        Neg_Xv->clear();
        Neg_Yv->clear();
        Neg_Zv->clear();
        Neg_Px->clear();
        Neg_Py->clear();
        Neg_Pz->clear();

        Pos_Xv->clear();
        Pos_Yv->clear();
        Pos_Zv->clear();
        Pos_Px->clear();
        Pos_Py->clear();
        Pos_Pz->clear();
    }
};

struct TypeA {
    std::vector<float>* Xv{nullptr};
    std::vector<float>* Yv{nullptr};
    std::vector<float>* Zv{nullptr};
    std::vector<float>* Px{nullptr};
    std::vector<float>* Py{nullptr};
    std::vector<float>* Pz{nullptr};
    std::vector<float>* Energy{nullptr};
    std::vector<float>* EnergyAsDecay{nullptr};

    std::vector<int>* Lambda_Index{nullptr};
    std::vector<float>* Lambda_Xv{nullptr};
    std::vector<float>* Lambda_Yv{nullptr};
    std::vector<float>* Lambda_Zv{nullptr};
    std::vector<float>* Lambda_Px{nullptr};
    std::vector<float>* Lambda_Py{nullptr};
    std::vector<float>* Lambda_Pz{nullptr};
    std::vector<float>* Lambda_E{nullptr};

    std::vector<int>* KaonZero_Index{nullptr};
    std::vector<float>* KaonZero_Xv{nullptr};
    std::vector<float>* KaonZero_Yv{nullptr};
    std::vector<float>* KaonZero_Zv{nullptr};
    std::vector<float>* KaonZero_Px{nullptr};
    std::vector<float>* KaonZero_Py{nullptr};
    std::vector<float>* KaonZero_Pz{nullptr};
    std::vector<float>* KaonZero_E{nullptr};

    void Clear() {
        Xv->clear();
        Yv->clear();
        Zv->clear();
        Px->clear();
        Py->clear();
        Pz->clear();
        Energy->clear();
        EnergyAsDecay->clear();

        Lambda_Index->clear();
        Lambda_Xv->clear();
        Lambda_Yv->clear();
        Lambda_Zv->clear();
        Lambda_Px->clear();
        Lambda_Py->clear();
        Lambda_Pz->clear();
        Lambda_E->clear();

        KaonZero_Index->clear();
        KaonZero_Xv->clear();
        KaonZero_Yv->clear();
        KaonZero_Zv->clear();
        KaonZero_Px->clear();
        KaonZero_Py->clear();
        KaonZero_Pz->clear();
        KaonZero_E->clear();
    }
};

struct TypeD {
    std::vector<float>* Xv{nullptr};
    std::vector<float>* Yv{nullptr};
    std::vector<float>* Zv{nullptr};
    std::vector<float>* Px{nullptr};
    std::vector<float>* Py{nullptr};
    std::vector<float>* Pz{nullptr};
    std::vector<float>* Energy{nullptr};
    std::vector<float>* EnergyAsDecay{nullptr};

    std::vector<int>* Lambda_Index{nullptr};
    std::vector<float>* Lambda_Xv{nullptr};
    std::vector<float>* Lambda_Yv{nullptr};
    std::vector<float>* Lambda_Zv{nullptr};
    std::vector<float>* Lambda_Px{nullptr};
    std::vector<float>* Lambda_Py{nullptr};
    std::vector<float>* Lambda_Pz{nullptr};
    std::vector<float>* Lambda_E{nullptr};

    std::vector<int>* Kaon_Entry{nullptr};
    std::vector<float>* Kaon_Xv{nullptr};
    std::vector<float>* Kaon_Yv{nullptr};
    std::vector<float>* Kaon_Zv{nullptr};
    std::vector<float>* Kaon_Px{nullptr};
    std::vector<float>* Kaon_Py{nullptr};
    std::vector<float>* Kaon_Pz{nullptr};
    std::vector<float>* Kaon_E{nullptr};

    void Clear() {
        Xv->clear();
        Yv->clear();
        Zv->clear();
        Px->clear();
        Py->clear();
        Pz->clear();
        Energy->clear();
        EnergyAsDecay->clear();

        Lambda_Index->clear();
        Lambda_Xv->clear();
        Lambda_Yv->clear();
        Lambda_Zv->clear();
        Lambda_Px->clear();
        Lambda_Py->clear();
        Lambda_Pz->clear();
        Lambda_E->clear();

        Kaon_Entry->clear();
        Kaon_Xv->clear();
        Kaon_Yv->clear();
        Kaon_Zv->clear();
        Kaon_Px->clear();
        Kaon_Py->clear();
        Kaon_Pz->clear();
        Kaon_E->clear();
    }
};

}  // namespace Tree2Secondaries::Output

#endif  // T2S_OUTPUT_FORMAT_HXX
