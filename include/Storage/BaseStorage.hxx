#pragma once

#include <cstddef>
#include <vector>

#include "Math/Constants.hxx"

namespace Tree2Secondaries::Storage {

// Flat //

struct Coordinates {
    float x{Const::DummyFloat};
    float y{Const::DummyFloat};
    float z{Const::DummyFloat};
};

struct Mom3 {
    float px{Const::DummyFloat};
    float py{Const::DummyFloat};
    float pz{Const::DummyFloat};
};

struct Mom4 {
    float px{Const::DummyFloat};
    float py{Const::DummyFloat};
    float pz{Const::DummyFloat};
    float energy{Const::DummyFloat};
};

struct State6 {
    float x{Const::DummyFloat};
    float y{Const::DummyFloat};
    float z{Const::DummyFloat};
    float px{Const::DummyFloat};
    float py{Const::DummyFloat};
    float pz{Const::DummyFloat};
};

// Vector //

struct VecCoordinates {
    std::vector<float>* x{nullptr};
    std::vector<float>* y{nullptr};
    std::vector<float>* z{nullptr};
};

struct VecMom3 {
    std::vector<float>* px{nullptr};
    std::vector<float>* py{nullptr};
    std::vector<float>* pz{nullptr};
};

struct VecMom4 {
    std::vector<float>* px{nullptr};
    std::vector<float>* py{nullptr};
    std::vector<float>* pz{nullptr};
    std::vector<float>* energy{nullptr};
};

struct VecState6 {
    std::vector<float>* x{nullptr};
    std::vector<float>* y{nullptr};
    std::vector<float>* z{nullptr};
    std::vector<float>* px{nullptr};
    std::vector<float>* py{nullptr};
    std::vector<float>* pz{nullptr};
};

template <std::size_t N>
struct VecCovMatrix {
    static constexpr std::size_t n_elements = N * (N + 1) / 2;
    std::vector<float>* mat{nullptr};
};

}  // namespace Tree2Secondaries::Storage
