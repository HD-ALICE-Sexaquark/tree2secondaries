#ifndef T2S_DF_EVENTS_HXX
#define T2S_DF_EVENTS_HXX

#include <vector>

#include <TTree.h>

#include "App/Utilities.hxx"
#include "DataFormats/DataFormats.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::DF {

namespace Flat {
struct alignas(T2S_SIMD_ALIGN) Event {
    Flat::Coordinates PV{};
    Flat::Coordinates MC_PV{};  // NOTE: not used when analyzing RD
    unsigned int RunNumber{0};
    unsigned int DirNumber{0};
    unsigned int DirNumberB{0};  // NOTE: not used when analyzing MC
    unsigned int EventNumber{0};
    float Centrality{0.};
    float MagneticField{0.};

    void CreateBranches_Event(TTree *tree, bool is_mc) {
        Utils::CreateBranch(tree, "RunNumber", &RunNumber);
        Utils::CreateBranch(tree, "DirNumber", &DirNumber);
        if (!is_mc) Utils::CreateBranch(tree, "DirNumberB", &DirNumberB);
        Utils::CreateBranch(tree, "EventNumber", &EventNumber);
        Utils::CreateBranch(tree, "Centrality", &Centrality);
        Utils::CreateBranch(tree, "MagneticField", &MagneticField);
        // primary vertex
        Utils::CreateBranch(tree, "PV_Xv", &PV.X);
        Utils::CreateBranch(tree, "PV_Yv", &PV.Y);
        Utils::CreateBranch(tree, "PV_Zv", &PV.Z);
        if (is_mc) {
            Utils::CreateBranch(tree, "MC_PV_Xv", &MC_PV.X);
            Utils::CreateBranch(tree, "MC_PV_Yv", &MC_PV.Y);
            Utils::CreateBranch(tree, "MC_PV_Zv", &MC_PV.Z);
        }
    }
    void ReadBranches_Event(TTree *tree, bool is_mc) {
        Utils::ReadBranch(tree, "RunNumber", &RunNumber);
        Utils::ReadBranch(tree, "DirNumber", &DirNumber);
        if (!is_mc) Utils::ReadBranch(tree, "DirNumberB", &DirNumberB);
        Utils::ReadBranch(tree, "EventNumber", &EventNumber);
        Utils::ReadBranch(tree, "Centrality", &Centrality);
        Utils::ReadBranch(tree, "MagneticField", &MagneticField);
        // primary vertex
        Utils::ReadBranch(tree, "PV_Xv", &PV.X);
        Utils::ReadBranch(tree, "PV_Yv", &PV.Y);
        Utils::ReadBranch(tree, "PV_Zv", &PV.Z);
        if (is_mc) {
            Utils::ReadBranch(tree, "MC_PV_Xv", &MC_PV.X);
            Utils::ReadBranch(tree, "MC_PV_Yv", &MC_PV.Y);
            Utils::ReadBranch(tree, "MC_PV_Zv", &MC_PV.Z);
        }
    }
};
}  // namespace Flat

namespace SOV {
struct alignas(T2S_SIMD_ALIGN) MC_Particles : SOV::States {
    std::vector<int> *PdgCode{nullptr};
    std::vector<int> *MotherEntry{nullptr};
    std::vector<int> *Status{nullptr};
    std::vector<int> *Generator{nullptr};
    std::vector<char> *IsPrimary{nullptr};
    std::vector<char> *IsSecFromMat{nullptr};
    std::vector<char> *IsSecFromWeak{nullptr};

    void ReadBranches_MCParticle(TTree *tree) {
        ReadBranches_States(tree, "MC");
        Utils::ReadBranch(tree, "MC_PdgCode", &PdgCode);
        Utils::ReadBranch(tree, "MC_Mother_McEntry", &MotherEntry);  // MAYBE: search a different name
        Utils::ReadBranch(tree, "MC_Status", &Status);
        Utils::ReadBranch(tree, "MC_Generator", &Generator);
        Utils::ReadBranch(tree, "MC_IsPrimary", &IsPrimary);
        Utils::ReadBranch(tree, "MC_IsSecFromMat", &IsSecFromMat);
        Utils::ReadBranch(tree, "MC_IsSecFromWeak", &IsSecFromWeak);
    }
};

struct alignas(T2S_SIMD_ALIGN) Tracks : SOV::States_NoE, SOV::CovMatrices_NoE {
    std::vector<int> *Charge{nullptr};
    std::vector<float> *DCAxy{nullptr};
    std::vector<float> *DCAz{nullptr};
    // pid info
    std::vector<float> *TPCSignal{nullptr};
    std::vector<float> *NSigmaPion{nullptr};
    std::vector<float> *NSigmaKaon{nullptr};
    std::vector<float> *NSigmaProton{nullptr};
    // mc info
    std::vector<int> *McEntry{nullptr};  // NOTE: not used when analyzing RD

    void ReadBranches_MCParticle(TTree *tree, bool is_mc) {
        ReadBranches_States_NoE(tree, "Track");
        ReadBranches_CovMatrices_NoE(tree, "Track");
        Utils::ReadBranch(tree, "Track_Charge", &Charge);
        Utils::ReadBranch(tree, "Track_DCAxy", &DCAxy);
        Utils::ReadBranch(tree, "Track_DCAz", &DCAz);
        // -- pid info
        Utils::ReadBranch(tree, "Track_TPCSignal", &TPCSignal);
        Utils::ReadBranch(tree, "Track_NSigmaPion", &NSigmaPion);
        Utils::ReadBranch(tree, "Track_NSigmaKaon", &NSigmaKaon);
        Utils::ReadBranch(tree, "Track_NSigmaProton", &NSigmaProton);
        // -- mc info
        if (is_mc) Utils::ReadBranch(tree, "Track_McEntry", &McEntry);
    }
};
}  // namespace SOV

}  // namespace Tree2Secondaries::DF

#endif  // T2S_DF_EVENTS_HXX
