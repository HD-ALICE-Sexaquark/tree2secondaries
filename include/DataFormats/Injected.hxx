#ifndef T2S_DF_INJECTED_HXX
#define T2S_DF_INJECTED_HXX

#include <vector>

#include <TTree.h>

#include "App/Utilities.hxx"
#include "DataFormats/DataFormats.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::DF {

namespace SOV {
struct alignas(T2S_SIMD_ALIGN) Injected {
    std::vector<int> *ReactionID{nullptr};
    std::vector<float> *X{nullptr};  // NOTE: not read but written by packager
    std::vector<float> *Y{nullptr};  // NOTE: not read but written by packager
    std::vector<float> *Z{nullptr};  // NOTE: not read but written by packager
    std::vector<float> *Px{nullptr};
    std::vector<float> *Py{nullptr};
    std::vector<float> *Pz{nullptr};
    std::vector<float> *Nucleon_Px{nullptr};
    std::vector<float> *Nucleon_Py{nullptr};
    std::vector<float> *Nucleon_Pz{nullptr};

    void Clear_SOV_Injected(bool include_coord) {
        ReactionID->clear();
        if (include_coord) {
            X->clear();
            Y->clear();
            Z->clear();
        }
        Px->clear();
        Py->clear();
        Pz->clear();
        Nucleon_Px->clear();
        Nucleon_Py->clear();
        Nucleon_Pz->clear();
    }
    void ReadBranches_SOV_Injected(TTree *tree, bool include_coord) {
        Utils::ReadBranch(tree, "ReactionID", &ReactionID);
        if (include_coord) {
            Utils::ReadBranch(tree, "SV_X", &X);
            Utils::ReadBranch(tree, "SV_Y", &Y);
            Utils::ReadBranch(tree, "SV_Z", &Z);
        }
        Utils::ReadBranch(tree, "Sexaquark_Px", &Px);
        Utils::ReadBranch(tree, "Sexaquark_Py", &Py);
        Utils::ReadBranch(tree, "Sexaquark_Pz", &Pz);
        Utils::ReadBranch(tree, "Nucleon_Px", &Nucleon_Px);
        Utils::ReadBranch(tree, "Nucleon_Py", &Nucleon_Py);
        Utils::ReadBranch(tree, "Nucleon_Pz", &Nucleon_Pz);
    }
    void CreateBranches_SOV_Injected(TTree *tree, bool include_coord) {
        tree->Branch("ReactionID", &ReactionID);
        if (include_coord) {
            tree->Branch("SV_X", &X);
            tree->Branch("SV_Y", &Y);
            tree->Branch("SV_Z", &Z);
        }
        tree->Branch("Sexaquark_Px", &Px);
        tree->Branch("Sexaquark_Py", &Py);
        tree->Branch("Sexaquark_Pz", &Pz);
        tree->Branch("Nucleon_Px", &Nucleon_Px);
        tree->Branch("Nucleon_Py", &Nucleon_Py);
        tree->Branch("Nucleon_Pz", &Nucleon_Pz);
    }
};
}  // namespace SOV

namespace Flat {
struct alignas(T2S_SIMD_ALIGN) Injected : Flat::Coordinates, Flat::LorentzVector {
    Flat::LorentzVector Nucleon;
    unsigned int RunNumber{};
    unsigned int DirNumber{};
    unsigned int EventNumber{};
    int ReactionID{};

    void CreateBranches_Flat_Injected(TTree *tree) {
        CreateBranches_Coordinates(tree, "SV");
        CreateBranches_LorentzVector(tree, "Sexaquark");
        Nucleon.CreateBranches_LorentzVector(tree, "Nucleon");
        tree->Branch("RunNumber", &RunNumber);
        tree->Branch("DirNumber", &DirNumber);
        tree->Branch("EventNumber", &EventNumber);
        tree->Branch("ReactionID", &ReactionID);
    }
};
}  // namespace Flat

}  // namespace Tree2Secondaries::DF

#endif  // T2S_DF_INJECTED_HXX
