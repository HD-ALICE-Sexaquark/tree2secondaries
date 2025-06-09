#ifndef T2S_UTILS_HXX
#define T2S_UTILS_HXX

#include <string>

#include "TTree.h"

namespace Tree2Secondaries::Utils {

template <typename T>
inline void ConnectBranch(TTree* tree, const std::string& branch_name, T address, const std::string& prefix = "") {
    tree->SetBranchStatus((prefix + branch_name).c_str(), true);
    tree->SetBranchAddress((prefix + branch_name).c_str(), address);
}

}  // namespace Tree2Secondaries::Utils

#endif  // T2S_UTILS_HXX
