#ifndef ALICE_MATH_HXX
#define ALICE_MATH_HXX

#include <cmath>
#include <utility>

namespace ALICE::Math {

// Based on https://stackoverflow.com/a/64247207
template <class S>
inline std::pair<S, S> sincos(S arg) {
    return {std::sin(arg), std::cos(arg)};
}

}  // namespace ALICE::Math

#endif  // ALICE_MATH_HXX
