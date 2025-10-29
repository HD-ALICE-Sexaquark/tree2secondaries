#pragma once

#include <string>

#include <TTree.h>

namespace Tree2Secondaries::Utils {

template <typename T>
inline void ReadBranch(TTree* tree, const std::string& branch_name, T address) {
    tree->SetBranchStatus(branch_name.c_str(), true);
    tree->SetBranchAddress(branch_name.c_str(), address);
}

}  // namespace Tree2Secondaries::Utils
