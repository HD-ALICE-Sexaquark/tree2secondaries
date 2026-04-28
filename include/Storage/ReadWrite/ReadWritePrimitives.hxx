#pragma once

#include <cstddef>
#include <format>
#include <string_view>

#include <TTree.h>

#include "App/Utilities.hxx"
#include "Storage/BaseStorage.hxx"

namespace Tree2Secondaries::IO {

// ## Flat ## //

// Coordinates //

inline void CreateBranches(TTree* t, Storage::Coordinates& c, std::string_view prefix) {
    Utils::CreateBranch(t, std::format("{}_X", prefix), &c.x);
    Utils::CreateBranch(t, std::format("{}_Y", prefix), &c.y);
    Utils::CreateBranch(t, std::format("{}_Z", prefix), &c.z);
}

inline void ReadBranches(TTree* t, Storage::Coordinates& c, std::string_view prefix) {
    Utils::ReadBranch(t, std::format("{}_X", prefix), &c.x);
    Utils::ReadBranch(t, std::format("{}_Y", prefix), &c.y);
    Utils::ReadBranch(t, std::format("{}_Z", prefix), &c.z);
}

// Mom3 //

inline void CreateBranches(TTree* t, Storage::Mom3& m, std::string_view prefix) {
    Utils::CreateBranch(t, std::format("{}_Px", prefix), &m.px);
    Utils::CreateBranch(t, std::format("{}_Py", prefix), &m.py);
    Utils::CreateBranch(t, std::format("{}_Pz", prefix), &m.pz);
}

inline void ReadBranches(TTree* t, Storage::Mom3& m, std::string_view prefix) {
    Utils::ReadBranch(t, std::format("{}_Px", prefix), &m.px);
    Utils::ReadBranch(t, std::format("{}_Py", prefix), &m.py);
    Utils::ReadBranch(t, std::format("{}_Pz", prefix), &m.pz);
}

// Mom4 //

inline void CreateBranches(TTree* t, Storage::Mom4& m, std::string_view prefix) {
    Utils::CreateBranch(t, std::format("{}_Px", prefix), &m.px);
    Utils::CreateBranch(t, std::format("{}_Py", prefix), &m.py);
    Utils::CreateBranch(t, std::format("{}_Pz", prefix), &m.pz);
    Utils::CreateBranch(t, std::format("{}_E", prefix), &m.energy);
}

inline void ReadBranches(TTree* t, Storage::Mom4& m, std::string_view prefix) {
    Utils::ReadBranch(t, std::format("{}_Px", prefix), &m.px);
    Utils::ReadBranch(t, std::format("{}_Py", prefix), &m.py);
    Utils::ReadBranch(t, std::format("{}_Pz", prefix), &m.pz);
    Utils::ReadBranch(t, std::format("{}_E", prefix), &m.energy);
}

// State6 //

inline void CreateBranches(TTree* t, Storage::State6& s, std::string_view prefix) {
    Utils::CreateBranch(t, std::format("{}_X", prefix), &s.x);
    Utils::CreateBranch(t, std::format("{}_Y", prefix), &s.y);
    Utils::CreateBranch(t, std::format("{}_Z", prefix), &s.z);
    Utils::CreateBranch(t, std::format("{}_Px", prefix), &s.px);
    Utils::CreateBranch(t, std::format("{}_Py", prefix), &s.py);
    Utils::CreateBranch(t, std::format("{}_Pz", prefix), &s.pz);
}

inline void ReadBranches(TTree* t, Storage::State6& s, std::string_view prefix) {
    Utils::ReadBranch(t, std::format("{}_X", prefix), &s.x);
    Utils::ReadBranch(t, std::format("{}_Y", prefix), &s.y);
    Utils::ReadBranch(t, std::format("{}_Z", prefix), &s.z);
    Utils::ReadBranch(t, std::format("{}_Px", prefix), &s.px);
    Utils::ReadBranch(t, std::format("{}_Py", prefix), &s.py);
    Utils::ReadBranch(t, std::format("{}_Pz", prefix), &s.pz);
}

// ## Vectors ## //

// VecCoordinates //

inline void CreateBranches(TTree* t, Storage::VecCoordinates& c, std::string_view prefix) {
    Utils::CreateBranch(t, std::format("{}_X", prefix), &c.x);
    Utils::CreateBranch(t, std::format("{}_Y", prefix), &c.y);
    Utils::CreateBranch(t, std::format("{}_Z", prefix), &c.z);
}

inline void ReadBranches(TTree* t, Storage::VecCoordinates& c, std::string_view prefix) {
    Utils::ReadBranch(t, std::format("{}_X", prefix), &c.x);
    Utils::ReadBranch(t, std::format("{}_Y", prefix), &c.y);
    Utils::ReadBranch(t, std::format("{}_Z", prefix), &c.z);
}

inline void ClearBranches(Storage::VecCoordinates& c) {
    c.x->clear();
    c.y->clear();
    c.z->clear();
}

// VecMom3 //

inline void CreateBranches(TTree* t, Storage::VecMom3& m, std::string_view prefix) {
    Utils::CreateBranch(t, std::format("{}_Px", prefix), &m.px);
    Utils::CreateBranch(t, std::format("{}_Py", prefix), &m.py);
    Utils::CreateBranch(t, std::format("{}_Pz", prefix), &m.pz);
}

inline void ReadBranches(TTree* t, Storage::VecMom3& m, std::string_view prefix) {
    Utils::ReadBranch(t, std::format("{}_Px", prefix), &m.px);
    Utils::ReadBranch(t, std::format("{}_Py", prefix), &m.py);
    Utils::ReadBranch(t, std::format("{}_Pz", prefix), &m.pz);
}

inline void ClearBranches(Storage::VecMom3& m) {
    m.px->clear();
    m.py->clear();
    m.pz->clear();
}

// VecMom4 //

inline void CreateBranches(TTree* t, Storage::VecMom4& m, std::string_view prefix) {
    Utils::CreateBranch(t, std::format("{}_Px", prefix), &m.px);
    Utils::CreateBranch(t, std::format("{}_Py", prefix), &m.py);
    Utils::CreateBranch(t, std::format("{}_Pz", prefix), &m.pz);
    Utils::CreateBranch(t, std::format("{}_E", prefix), &m.energy);
}

inline void ReadBranches(TTree* t, Storage::VecMom4& m, std::string_view prefix) {
    Utils::ReadBranch(t, std::format("{}_Px", prefix), &m.px);
    Utils::ReadBranch(t, std::format("{}_Py", prefix), &m.py);
    Utils::ReadBranch(t, std::format("{}_Pz", prefix), &m.pz);
    Utils::ReadBranch(t, std::format("{}_E", prefix), &m.energy);
}

inline void ClearBranches(Storage::VecMom4& m) {
    m.px->clear();
    m.py->clear();
    m.pz->clear();
    m.energy->clear();
}

// VecState6 //

inline void CreateBranches(TTree* t, Storage::VecState6& s, std::string_view prefix) {
    Utils::CreateBranch(t, std::format("{}_X", prefix), &s.x);
    Utils::CreateBranch(t, std::format("{}_Y", prefix), &s.y);
    Utils::CreateBranch(t, std::format("{}_Z", prefix), &s.z);
    Utils::CreateBranch(t, std::format("{}_Px", prefix), &s.px);
    Utils::CreateBranch(t, std::format("{}_Py", prefix), &s.py);
    Utils::CreateBranch(t, std::format("{}_Pz", prefix), &s.pz);
}

inline void ReadBranches(TTree* t, Storage::VecState6& s, std::string_view prefix) {
    Utils::ReadBranch(t, std::format("{}_X", prefix), &s.x);
    Utils::ReadBranch(t, std::format("{}_Y", prefix), &s.y);
    Utils::ReadBranch(t, std::format("{}_Z", prefix), &s.z);
    Utils::ReadBranch(t, std::format("{}_Px", prefix), &s.px);
    Utils::ReadBranch(t, std::format("{}_Py", prefix), &s.py);
    Utils::ReadBranch(t, std::format("{}_Pz", prefix), &s.pz);
}

inline void ClearBranches(Storage::VecState6& s) {
    s.x->clear();
    s.y->clear();
    s.z->clear();
    s.px->clear();
    s.py->clear();
    s.pz->clear();
}

// Cov Matrix //

template <std::size_t N>
inline void CreateBranches(TTree* t, Storage::VecCovMatrix<N>& cov, std::string_view prefix) {
    Utils::CreateBranch(t, std::format("{}_CovMatrix", prefix), &cov.mat);
}

template <std::size_t N>
inline void ReadBranches(TTree* t, Storage::VecCovMatrix<N>& cov, std::string_view prefix) {
    Utils::ReadBranch(t, std::format("{}_CovMatrix", prefix), &cov.mat);
}

template <std::size_t N>
inline void ClearBranches(Storage::VecCovMatrix<N>& cov) {
    cov.mat->clear();
}

}  // namespace Tree2Secondaries::IO
