#pragma once

#include "App/Utilities.hxx"
#include "Math/Constants.hxx"

class TTree;

namespace Tree2Secondaries::Schema {

struct Coordinates {
    // Storage //
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

    // Member Variables //
    float x{Const::DummyFloat};
    float y{Const::DummyFloat};
    float z{Const::DummyFloat};
};

struct Mom3 {
    // Storage //
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

    // Member Variables //
    float px{Const::DummyFloat};
    float py{Const::DummyFloat};
    float pz{Const::DummyFloat};
};

struct Mom4 {
    // Storage //
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

    // Member Variables //
    float px{Const::DummyFloat};
    float py{Const::DummyFloat};
    float pz{Const::DummyFloat};
    float energy{Const::DummyFloat};
};

struct State6 {
    // Storage //
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

    // Member Variables //
    float x{Const::DummyFloat};
    float y{Const::DummyFloat};
    float z{Const::DummyFloat};
    float px{Const::DummyFloat};
    float py{Const::DummyFloat};
    float pz{Const::DummyFloat};
};

}  // namespace Tree2Secondaries::Schema
