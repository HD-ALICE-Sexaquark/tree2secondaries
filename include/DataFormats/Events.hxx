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
    std::vector<long> *MotherEntry{nullptr};
    std::vector<unsigned int> *Status{nullptr};
    std::vector<int> *Generator{nullptr};
    std::vector<bool> *IsPrimary{nullptr};
    std::vector<bool> *IsSecFromMat{nullptr};
    std::vector<bool> *IsSecFromWeak{nullptr};
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

    std::vector<float> *Alpha{nullptr};
    std::vector<float> *Snp{nullptr};
    std::vector<float> *Tgl{nullptr};
    std::vector<float> *Signed1Pt{nullptr};
    std::vector<float> *SigmaY2{nullptr};
    std::vector<float> *SigmaZY{nullptr};
    std::vector<float> *SigmaZ2{nullptr};
    std::vector<float> *SigmaSnpY{nullptr};
    std::vector<float> *SigmaSnpZ{nullptr};
    std::vector<float> *SigmaSnp2{nullptr};
    std::vector<float> *SigmaTglY{nullptr};
    std::vector<float> *SigmaTglZ{nullptr};
    std::vector<float> *SigmaTglSnp{nullptr};
    std::vector<float> *SigmaTgl2{nullptr};
    std::vector<float> *Sigma1PtY{nullptr};
    std::vector<float> *Sigma1PtZ{nullptr};
    std::vector<float> *Sigma1PtSnp{nullptr};
    std::vector<float> *Sigma1PtTgl{nullptr};
    std::vector<float> *Sigma1Pt2{nullptr};

    std::vector<long> *McEntry{nullptr};
};

}  // namespace Tree2Secondaries::Events

#endif  // T2S_STRUCTS_EVENTS_HXX
