#ifndef T2S_INPUT_FORMAT_HXX
#define T2S_INPUT_FORMAT_HXX

#include <vector>

namespace Tree2Secondaries::Struct {

struct Event {
    unsigned int RunNumber{0};
    unsigned int DirNumber{0};
    unsigned int EventNumber{0};
    float Centrality{0.};
    float MagneticField{0.};
    float PV_Xv{0.};
    float PV_Yv{0.};
    float PV_Zv{0.};

    float PV_TrueXv{0.};
    float PV_TrueYv{0.};
    float PV_TrueZv{0.};
};

}  // namespace Tree2Secondaries::Struct

namespace Tree2Secondaries::Input {

struct Injected {
    std::vector<int> *ReactionID{nullptr};
    std::vector<float> *Px{nullptr};
    std::vector<float> *Py{nullptr};
    std::vector<float> *Pz{nullptr};
    std::vector<float> *Nucleon_Px{nullptr};
    std::vector<float> *Nucleon_Py{nullptr};
    std::vector<float> *Nucleon_Pz{nullptr};
};

struct MC {
    std::vector<int> *PdgCode{nullptr};
    std::vector<int> *Mother_McEntry{nullptr};
    std::vector<int> *Status{nullptr};
    std::vector<int> *Generator{nullptr};
};

struct Tracks {
    std::vector<float> *Px{nullptr};
    std::vector<float> *Py{nullptr};
    std::vector<float> *Pz{nullptr};
    std::vector<float> *X{nullptr};
    std::vector<float> *Y{nullptr};
    std::vector<float> *Z{nullptr};
    std::vector<int> *Charge{nullptr};
    std::vector<float> *NSigmaPion{nullptr};
    std::vector<float> *NSigmaKaon{nullptr};
    std::vector<float> *NSigmaProton{nullptr};

    std::vector<int> *McEntry{nullptr};
};

}  // namespace Tree2Secondaries::Input

#endif  // T2S_INPUT_FORMAT_HXX
