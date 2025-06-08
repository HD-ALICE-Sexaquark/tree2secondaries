#ifndef T2S_STRUCTS_CHANNEL_D_HXX
#define T2S_STRUCTS_CHANNEL_D_HXX

#include <vector>

namespace Tree2Secondaries::Found {

struct ChannelD {
    std::vector<float>* X{nullptr};
    std::vector<float>* Y{nullptr};
    std::vector<float>* Z{nullptr};
    std::vector<float>* Px{nullptr};
    std::vector<float>* Py{nullptr};
    std::vector<float>* Pz{nullptr};
    std::vector<float>* E{nullptr};
    std::vector<float>* E_NoNucleon{nullptr};

    std::vector<int>* Lambda_Index{nullptr};
    std::vector<float>* Lambda_X{nullptr};
    std::vector<float>* Lambda_Y{nullptr};
    std::vector<float>* Lambda_Z{nullptr};
    std::vector<float>* Lambda_Px{nullptr};
    std::vector<float>* Lambda_Py{nullptr};
    std::vector<float>* Lambda_Pz{nullptr};
    std::vector<float>* Lambda_E{nullptr};

    std::vector<int>* Kaon_Entry{nullptr};
    std::vector<float>* Kaon_X{nullptr};
    std::vector<float>* Kaon_Y{nullptr};
    std::vector<float>* Kaon_Z{nullptr};
    std::vector<float>* Kaon_Px{nullptr};
    std::vector<float>* Kaon_Py{nullptr};
    std::vector<float>* Kaon_Pz{nullptr};
    std::vector<float>* Kaon_E{nullptr};

    void Clear(bool is_mc = false) {
        X->clear();
        Y->clear();
        Z->clear();
        Px->clear();
        Py->clear();
        Pz->clear();
        E->clear();
        E_NoNucleon->clear();

        Lambda_Index->clear();
        Lambda_X->clear();
        Lambda_Y->clear();
        Lambda_Z->clear();
        Lambda_Px->clear();
        Lambda_Py->clear();
        Lambda_Pz->clear();
        Lambda_E->clear();

        Kaon_Entry->clear();
        Kaon_X->clear();
        Kaon_Y->clear();
        Kaon_Z->clear();
        Kaon_Px->clear();
        Kaon_Py->clear();
        Kaon_Pz->clear();
        Kaon_E->clear();
    }
};

}  // namespace Tree2Secondaries::Found

#endif  // T2S_STRUCTS_CHANNEL_D_HXX
