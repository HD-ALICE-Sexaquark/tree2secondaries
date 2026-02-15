#pragma once

#include <cstdlib>

#include "Math/Constants.hxx"

namespace Tree2Secondaries::View {

template <typename SourceType, typename IndexType>
struct Base {
    [[nodiscard]] size_t EntryAsSize() const { return static_cast<size_t>(Entry); }

    const SourceType* Source{};
    IndexType Entry{};
};

template <typename SourceType>
bool IsValid(const View::Base<SourceType, int>& view) {
    return view.Entry > Const::DummyInt;
}

}  // namespace Tree2Secondaries::View
