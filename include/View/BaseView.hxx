#pragma once

#include "Math/Constants.hxx"

namespace Tree2Secondaries::View {

template <typename SourceType>
struct Base {
    const SourceType* Source{};
    int Entry{};

    [[nodiscard]] bool IsValid() const { return Entry > Const::DummyInt; }
};

}  // namespace Tree2Secondaries::View
