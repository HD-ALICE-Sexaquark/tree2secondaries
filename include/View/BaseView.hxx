#pragma once

#include "Math/Constants.hxx"

namespace Tree2Secondaries::View {

template <typename SourceType, typename IndexType>
struct Base {
    [[nodiscard]] unsigned int EntryAsSize() const { return static_cast<unsigned int>(Entry); }

    const SourceType* Source{};
    IndexType Entry{};
};

template <typename SourceType>
bool IsValid(const View::Base<SourceType, int>& view) {
    return view.Entry > Const::DummyInt;
}

}  // namespace Tree2Secondaries::View
