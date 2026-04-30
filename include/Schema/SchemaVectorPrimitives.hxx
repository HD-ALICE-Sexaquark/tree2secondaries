#pragma once

#include <cstddef>
#include <vector>

#include "App/Utilities.hxx"
#include "Math/Constants.hxx"

class TTree;

namespace Tree2Secondaries::Schema {

struct VecCoordinates {

    // ## Storage ## //

    void CreateBranches(TTree* t, std::string_view prefix) {
        Utils::CreateBranch(t, std::format("{}_X", prefix), &x);
        Utils::CreateBranch(t, std::format("{}_Y", prefix), &y);
        Utils::CreateBranch(t, std::format("{}_Z", prefix), &z);
    }

    void ReadBranches(TTree* t, std::string_view prefix) {
        Utils::ReadBranch(t, std::format("{}_X", prefix), &x);
        Utils::ReadBranch(t, std::format("{}_Y", prefix), &y);
        Utils::ReadBranch(t, std::format("{}_Z", prefix), &z);
    }

    void PushDummy() {
        x->push_back(Const::DummyFloat);
        y->push_back(Const::DummyFloat);
        z->push_back(Const::DummyFloat);
    }

    void ClearBranches() {
        x->clear();
        y->clear();
        z->clear();
    }

    // ## Member Variables ## //

    std::vector<float>* x{nullptr};
    std::vector<float>* y{nullptr};
    std::vector<float>* z{nullptr};
};

struct VecMom3 {

    // ## Storage ## //

    void CreateBranches(TTree* t, std::string_view prefix) {
        Utils::CreateBranch(t, std::format("{}_Px", prefix), &px);
        Utils::CreateBranch(t, std::format("{}_Py", prefix), &py);
        Utils::CreateBranch(t, std::format("{}_Pz", prefix), &pz);
    }

    void ReadBranches(TTree* t, std::string_view prefix) {
        Utils::ReadBranch(t, std::format("{}_Px", prefix), &px);
        Utils::ReadBranch(t, std::format("{}_Py", prefix), &py);
        Utils::ReadBranch(t, std::format("{}_Pz", prefix), &pz);
    }

    void PushDummy() {
        px->push_back(Const::DummyFloat);
        py->push_back(Const::DummyFloat);
        pz->push_back(Const::DummyFloat);
    }

    void ClearBranches() {
        px->clear();
        py->clear();
        pz->clear();
    }

    // ## Member Variables ## //

    std::vector<float>* px{nullptr};
    std::vector<float>* py{nullptr};
    std::vector<float>* pz{nullptr};
};

struct VecMom4 {

    // ## Storage ## //

    void CreateBranches(TTree* t, std::string_view prefix) {
        Utils::CreateBranch(t, std::format("{}_Px", prefix), &px);
        Utils::CreateBranch(t, std::format("{}_Py", prefix), &py);
        Utils::CreateBranch(t, std::format("{}_Pz", prefix), &pz);
        Utils::CreateBranch(t, std::format("{}_E", prefix), &energy);
    }

    void ReadBranches(TTree* t, std::string_view prefix) {
        Utils::ReadBranch(t, std::format("{}_Px", prefix), &px);
        Utils::ReadBranch(t, std::format("{}_Py", prefix), &py);
        Utils::ReadBranch(t, std::format("{}_Pz", prefix), &pz);
        Utils::ReadBranch(t, std::format("{}_E", prefix), &energy);
    }

    void PushDummy() {
        px->push_back(Const::DummyFloat);
        py->push_back(Const::DummyFloat);
        pz->push_back(Const::DummyFloat);
        energy->push_back(Const::DummyFloat);
    }

    void ClearBranches() {
        px->clear();
        py->clear();
        pz->clear();
        energy->clear();
    }

    // ## Member Variables ## //

    std::vector<float>* px{nullptr};
    std::vector<float>* py{nullptr};
    std::vector<float>* pz{nullptr};
    std::vector<float>* energy{nullptr};
};

struct VecState6 {

    // ## Storage ## //

    void CreateBranches(TTree* t, std::string_view prefix) {
        Utils::CreateBranch(t, std::format("{}_X", prefix), &x);
        Utils::CreateBranch(t, std::format("{}_Y", prefix), &y);
        Utils::CreateBranch(t, std::format("{}_Z", prefix), &z);
        Utils::CreateBranch(t, std::format("{}_Px", prefix), &px);
        Utils::CreateBranch(t, std::format("{}_Py", prefix), &py);
        Utils::CreateBranch(t, std::format("{}_Pz", prefix), &pz);
    }

    void ReadBranches(TTree* t, std::string_view prefix) {
        Utils::ReadBranch(t, std::format("{}_X", prefix), &x);
        Utils::ReadBranch(t, std::format("{}_Y", prefix), &y);
        Utils::ReadBranch(t, std::format("{}_Z", prefix), &z);
        Utils::ReadBranch(t, std::format("{}_Px", prefix), &px);
        Utils::ReadBranch(t, std::format("{}_Py", prefix), &py);
        Utils::ReadBranch(t, std::format("{}_Pz", prefix), &pz);
    }

    void PushDummy() {
        x->push_back(Const::DummyFloat);
        y->push_back(Const::DummyFloat);
        z->push_back(Const::DummyFloat);
        px->push_back(Const::DummyFloat);
        py->push_back(Const::DummyFloat);
        pz->push_back(Const::DummyFloat);
    }

    void ClearBranches() {
        x->clear();
        y->clear();
        z->clear();
        px->clear();
        py->clear();
        pz->clear();
    }

    // ## Member Variables ## //

    std::vector<float>* x{nullptr};
    std::vector<float>* y{nullptr};
    std::vector<float>* z{nullptr};
    std::vector<float>* px{nullptr};
    std::vector<float>* py{nullptr};
    std::vector<float>* pz{nullptr};
};

template <std::size_t N>
struct VecCovMatrix {

    // ## Storage ## //

    void CreateBranches(TTree* t, std::string_view prefix) { Utils::CreateBranch(t, std::format("{}_CovMatrix", prefix), &mat); }

    void ReadBranches(TTree* t, std::string_view prefix) { Utils::ReadBranch(t, std::format("{}_CovMatrix", prefix), &mat); }

    void Push(const std::array<float, N*(N + 1) / 2>& cov_matrix) {
        for (std::size_t i = 0; i < n_elements; ++i) {
            mat->push_back(cov_matrix[i]);
        }
    }

    void ClearBranches() { mat->clear(); }

    // ## Member Variables ## //

    static constexpr std::size_t n_elements = N * (N + 1) / 2;
    std::vector<float>* mat{nullptr};
};

}  // namespace Tree2Secondaries::Schema
