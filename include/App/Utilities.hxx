#ifndef T2S_UTILS_HXX
#define T2S_UTILS_HXX

#include <string>

#include <TTree.h>

namespace Tree2Secondaries::Utils {

template <typename T>
inline void ConnectBranch(TTree* tree, const std::string& branch_name, T address) {
    tree->SetBranchStatus(branch_name.c_str(), true);
    tree->SetBranchAddress(branch_name.c_str(), address);
}

template <typename T>
inline void CreateBranch(TTree* tree, const std::string& branch_name, T address) {
    tree->Branch(branch_name.c_str(), address);
}

}  // namespace Tree2Secondaries::Utils

#endif  // T2S_UTILS_HXX
