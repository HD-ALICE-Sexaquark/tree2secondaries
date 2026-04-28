#pragma once

#include <array>
#include <cstddef>
#include <format>
#include <string_view>

#include <TTree.h>

#include <Eigen/Core>

// # Print Utilities # //

// For std::arrays //

template <typename T, std::size_t N>
struct std::formatter<std::array<T, N>> {
    constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }
    auto format(const std::array<T, N>& arr, std::format_context& ctx) const {
        auto out = ctx.out();
        out = std::format_to(out, "(");
        for (std::size_t i = 0; i < N; ++i) {
            out = std::format_to(out, "{:13.6e}", arr[i]);
            if (i < N - 1) out = std::format_to(out, ", ");
        }
        out = std::format_to(out, ")");
        return out;
    }
};

// For Eigen Matrices //

template <typename T, int NRows, int NCols, int Options, int MaxRows, int MaxCols>
struct std::formatter<Eigen::Matrix<T, NRows, NCols, Options, MaxRows, MaxCols>> {

    constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

    auto format(const Eigen::Matrix<T, NRows, NCols, Options, MaxRows, MaxCols>& mat, std::format_context& ctx) const {
        auto out = ctx.out();
        if constexpr (NCols == 1) {
            // Vector: print as a row "(x, y, z)"
            out = std::format_to(out, "(");
            for (int i = 0; i < NRows; ++i) {
                out = std::format_to(out, "{:13.6e}", mat(i));
                if (i < NRows - 1) out = std::format_to(out, ", ");
            }
            out = std::format_to(out, ")");
        } else {
            // Matrix: print row-by-row
            out = std::format_to(out, "\n");
            for (int i = 0; i < NRows; ++i) {
                for (int j = 0; j < NCols; ++j) {
                    out = std::format_to(out, "{:13.6e}", mat(i, j));
                    if (j < NCols - 1)
                        out = std::format_to(out, "   ");
                    else if (i < NRows - 1)
                        out = std::format_to(out, "\n");
                }
            }
        }
        return out;
    }
};

// Utilities //

namespace Tree2Secondaries::Utils {

template <typename T>
inline void ReadBranch(TTree* tree, std::string_view branch_name, T* address) {
    tree->SetBranchStatus(std::string{branch_name}.c_str(), true);
    tree->SetBranchAddress(std::string{branch_name}.c_str(), address);
}

template <typename T>
inline void CreateBranch(TTree* tree, std::string_view branch_name, T* address) {
    tree->Branch(std::string{branch_name}.c_str(), address);
}

}  // namespace Tree2Secondaries::Utils
