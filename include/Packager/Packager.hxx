#ifndef T2S_PACKAGER_HXX
#define T2S_PACKAGER_HXX

#include <memory>
#include <utility>

#include "TChain.h"
#include "TFile.h"

#include "App/Settings.hxx"
#include "Math/Constants.hxx"
#include "Math/KFWrapper.hxx"
#include "Packager/TruthHandler.hxx"
#include "Structures/Events.hxx"
#include "Structures/PackedEvents.hxx"

namespace Tree2Secondaries {

// Pack secondary V0s and tracks.
class Packager {
   public:
    Packager(const Packager &) = delete;
    Packager(Packager &&) = delete;
    Packager &operator=(const Packager &) = delete;
    Packager &operator=(Packager &&) = delete;
    ~Packager() = default;

    explicit Packager(Settings settings) : fSettings{std::move(settings)} {}

    ReactionChannel GetReactionChannel() const { return fSettings.Channel; }

    bool Initialize();
    void ConnectInputBranches();

    void ConnectBranchesEvents();
    void ConnectBranchesInjected();
    void ConnectBranchesMC();
    void ConnectBranchesTracks();

    bool PrepareOutputFile();
    bool PrepareOutputTree();
    void CreateOutputBranches();

    void CreateOutputBranchesEvents();
    void CreateOutputBranchesInjected();
    void CreateOutputBranchesV0s(const std::string &name_v0, PackedEvents::V0s &out_branches);
    void CreateOutputBranchesTracks(const std::string &name_part, PackedEvents::Tracks &out_branches);

    long long NumberEventsToRead() const { return fSettings.LimitToNEvents ? fSettings.LimitToNEvents : fEventsTree->GetEntries(); }
    bool IsMC() const { return fSettings.IsMC; }
    bool IsSignalMC() const { return fSettings.IsSignalMC; }
    void GetEvent(long long i_event) { fEventsTree->GetEntry(i_event); }

    size_t NumberInjected() const { return fInput_Injected.ReactionID->size(); }
    size_t NumberTracks() const { return fInput_Tracks.Px->size(); }

    void ProcessEvent();
    void ProcessInjected();
    void ProcessMC();
    void ProcessTracks();

    void PackTracks(PdgCode pdg_code);

    void FindV0s(PdgCode pdg_code);
    bool PassesCuts(const KF::V0 &v0, PdgCode pdg_code) const {
        switch (pdg_code) {
            case PdgCode::AntiLambda:
            case PdgCode::Lambda:
                return PassesCuts_Lambda(v0);
            case PdgCode::KaonZeroShort:
                return PassesCuts_KaonZeroShort(v0);
            default:
                return false;
        }
    }

    void EndOfEvent();
    void EndOfAnalysis();

   private:
    std::array<double, 6> PackParams_KF(size_t esd_idx) { return KF::PackParams(fInput_Tracks, esd_idx); }
    std::array<float, 5> PackParams_ALICE(size_t esd_idx) { return ALICE::PackParams(fInput_Tracks, esd_idx); }
    std::array<float, 15> PackCovMatrix_ALICE(size_t esd_idx) { return ALICE::PackCovMatrix(fInput_Tracks, esd_idx); }

    void Store(const KF::Track &track, PackedEvents::Tracks &out_branches);

    bool PassesCuts_Lambda(const KF::V0 &v0) const;
    bool PassesCuts_KaonZeroShort(const KF::V0 &v0) const;
    void Store(const KF::V0 &v0, PackedEvents::V0s &out_branches);

    Settings fSettings;
    std::unique_ptr<TChain> fEventsTree;

    std::unique_ptr<TFile> fOutputFile;
    std::unique_ptr<TTree> fOutputTree;

    // input structs //

    Events::Event fInput_Event;
    Events::Injected fInput_Injected;
    Events::Tracks fInput_Tracks;

    // helpers //

    TruthHandler fTruthHandler;

    // indices //

    std::vector<size_t> fVec_AntiProtons;
    std::vector<size_t> fVec_Protons;
    std::vector<size_t> fVec_NegKaons;
    std::vector<size_t> fVec_PosKaons;
    std::vector<size_t> fVec_PiMinus;
    std::vector<size_t> fVec_PiPlus;

    // output structs //

    Events::Event fOutput_Event;
    Events::Injected fOutput_Injected;

    PackedEvents::V0s fOutput_AntiLambdas;
    PackedEvents::V0s fOutput_Lambdas;
    PackedEvents::V0s fOutput_KaonsZeroShort;

    PackedEvents::Tracks fOutput_NegKaons;
    PackedEvents::Tracks fOutput_PosKaons;
    PackedEvents::Tracks fOutput_PiMinus;
    PackedEvents::Tracks fOutput_PiPlus;
};

}  // namespace Tree2Secondaries

#endif  // T2S_PACKAGER_HXX
