#pragma once

#include "App/Utilities.hxx"
#include "Math/Constants.hxx"
#include "Storage/Vector/BaseVector.hxx"

namespace Tree2Secondaries::Storage::Vector {

// `Vector::States` +
// `PdgCode` + `MotherEntry` + `Status` + `Generator` + `IsPrimary` + `IsSecFromMat` + `IsSecFromWeak`.
struct alignas(T2S_SIMD_ALIGN) MCParticles : Vector::States {
    std::vector<int> *PdgCode{nullptr};
    std::vector<int> *MotherEntry{nullptr};
    std::vector<int> *Status{nullptr};
    std::vector<int> *Generator{nullptr};
    std::vector<char> *IsPrimary{nullptr};
    std::vector<char> *IsSecFromMat{nullptr};
    std::vector<char> *IsSecFromWeak{nullptr};

    void ReadBranches_VectorMCParticles(TTree *tree) {
        ReadBranches_VectorStates(tree, "MC");
        Utils::ReadBranch(tree, "MC_PdgCode", &PdgCode);
        Utils::ReadBranch(tree, "MC_Mother_McEntry", &MotherEntry);
        Utils::ReadBranch(tree, "MC_Status", &Status);
        Utils::ReadBranch(tree, "MC_Generator", &Generator);
        Utils::ReadBranch(tree, "MC_IsPrimary", &IsPrimary);
        Utils::ReadBranch(tree, "MC_IsSecFromMat", &IsSecFromMat);
        Utils::ReadBranch(tree, "MC_IsSecFromWeak", &IsSecFromWeak);
    }
};

}  // namespace Tree2Secondaries::Storage::Vector
