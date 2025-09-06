#ifndef T2S_STRUCTS_EVENTS_HXX
#define T2S_STRUCTS_EVENTS_HXX

#include <vector>

namespace Tree2Secondaries::Events {

struct alignas(32) Event {
    unsigned int RunNumber{0};
    unsigned int DirNumber{0};
    unsigned int DirNumberB{0};
    unsigned int EventNumber{0};
    float Centrality{0.};
    float MagneticField{0.};
    float PV_Xv{0.};
    float PV_Yv{0.};
    float PV_Zv{0.};

    float MC_PV_Xv{0.};
    float MC_PV_Yv{0.};
    float MC_PV_Zv{0.};
};

struct alignas(32) Injected {
    std::vector<int> *ReactionID{nullptr};
    std::vector<float> *X{nullptr};
    std::vector<float> *Y{nullptr};
    std::vector<float> *Z{nullptr};
    std::vector<float> *Px{nullptr};
    std::vector<float> *Py{nullptr};
    std::vector<float> *Pz{nullptr};
    std::vector<float> *Nucleon_Px{nullptr};
    std::vector<float> *Nucleon_Py{nullptr};
    std::vector<float> *Nucleon_Pz{nullptr};

    void Clear() {
        ReactionID->clear();
        Px->clear();
        Py->clear();
        Pz->clear();
        Nucleon_Px->clear();
        Nucleon_Py->clear();
        Nucleon_Pz->clear();
    }
};

struct alignas(32) MC {
    std::vector<float> *X{nullptr};
    std::vector<float> *Y{nullptr};
    std::vector<float> *Z{nullptr};
    std::vector<float> *Px{nullptr};
    std::vector<float> *Py{nullptr};
    std::vector<float> *Pz{nullptr};
    std::vector<float> *E{nullptr};

    std::vector<int> *PdgCode{nullptr};
    std::vector<int> *MotherEntry{nullptr};
    std::vector<int> *Status{nullptr};
    std::vector<int> *Generator{nullptr};
    std::vector<char> *IsPrimary{nullptr};
    std::vector<char> *IsSecFromMat{nullptr};
    std::vector<char> *IsSecFromWeak{nullptr};
};

struct alignas(32) Tracks {
    std::vector<float> *X{nullptr};
    std::vector<float> *Y{nullptr};
    std::vector<float> *Z{nullptr};
    std::vector<float> *Px{nullptr};
    std::vector<float> *Py{nullptr};
    std::vector<float> *Pz{nullptr};
    std::vector<int> *Charge{nullptr};
    std::vector<float> *NSigmaPion{nullptr};
    std::vector<float> *NSigmaKaon{nullptr};
    std::vector<float> *NSigmaProton{nullptr};

    std::vector<float> *SigmaX2{nullptr};
    std::vector<float> *SigmaXY{nullptr};
    std::vector<float> *SigmaY2{nullptr};
    std::vector<float> *SigmaXZ{nullptr};
    std::vector<float> *SigmaYZ{nullptr};
    std::vector<float> *SigmaZ2{nullptr};
    std::vector<float> *SigmaXPx{nullptr};
    std::vector<float> *SigmaYPx{nullptr};
    std::vector<float> *SigmaZPx{nullptr};
    std::vector<float> *SigmaPx2{nullptr};
    std::vector<float> *SigmaXPy{nullptr};
    std::vector<float> *SigmaYPy{nullptr};
    std::vector<float> *SigmaZPy{nullptr};
    std::vector<float> *SigmaPxPy{nullptr};
    std::vector<float> *SigmaPy2{nullptr};
    std::vector<float> *SigmaXPz{nullptr};
    std::vector<float> *SigmaYPz{nullptr};
    std::vector<float> *SigmaZPz{nullptr};
    std::vector<float> *SigmaPxPz{nullptr};
    std::vector<float> *SigmaPyPz{nullptr};
    std::vector<float> *SigmaPz2{nullptr};

    std::vector<int> *McEntry{nullptr};
};

}  // namespace Tree2Secondaries::Events

#endif  // T2S_STRUCTS_EVENTS_HXX
