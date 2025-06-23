#ifndef T2S_PACKAGER_HXX
#define T2S_PACKAGER_HXX

#include <memory>
#include <utility>

#include "TChain.h"
#include "TFile.h"

#include "App/Settings.hxx"
#include "DataFormats/Events.hxx"
#include "DataFormats/PackedEvents.hxx"
#include "KF/Extensions.hxx"
#include "MC/Particles.hxx"
#include "Math/Constants.hxx"

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

    void ConnectBranches_Events();
    void ConnectBranches_Injected();
    void ConnectBranches_MC();
    void ConnectBranches_Tracks();

    bool PrepareOutputFile();
    bool PrepareOutputTree();
    void CreateOutputBranches();

    void CreateOutputBranches_Events();
    void CreateOutputBranches_Injected();
    void CreateOutputBranches_V0s(const std::string &name_v0, PackedEvents::V0s &sov);
    void CreateOutputBranches_Tracks(const std::string &name_part, PackedEvents::Tracks &sov);
    void CreateOutputBranches_MC_V0s(const std::string &name_v0, PackedEvents::MC_V0s &sov);
    void CreateOutputBranches_MC_Tracks(const std::string &name_part, PackedEvents::MC_Tracks &sov);

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
    void Store(const KF::Track &track, PackedEvents::Tracks &sov);
    void StoreMC(const MC::Track &track, PackedEvents::MC_Tracks &sov);

    bool PassesCuts_Lambda(const KF::V0 &v0) const;
    bool PassesCuts_KaonZeroShort(const KF::V0 &v0) const;
    void Store(const KF::V0 &v0, PackedEvents::V0s &sov);
    void StoreMC(const MC::V0 &v0, PackedEvents::MC_V0s &sov);

    Settings fSettings;
    std::unique_ptr<TChain> fEventsTree;

    std::unique_ptr<TFile> fOutputFile;
    std::unique_ptr<TTree> fOutputTree;

    // input branches //

    Events::Event fInput_Event;
    Events::Injected fInput_Injected;

    Events::MC fInput_MC;
    Events::Tracks fInput_Tracks;

    // indices -- temporary containers //

    std::vector<size_t> fVec_AntiProtons;
    std::vector<size_t> fVec_Protons;
    std::vector<size_t> fVec_NegKaons;
    std::vector<size_t> fVec_PosKaons;
    std::vector<size_t> fVec_PiMinus;
    std::vector<size_t> fVec_PiPlus;

    // output branches //

    Events::Event fOutput_Event;
    Events::Injected fOutput_Injected;

    PackedEvents::MC_V0s fOutput_MC_AntiLambdas;
    PackedEvents::MC_V0s fOutput_MC_Lambdas;
    PackedEvents::MC_V0s fOutput_MC_KaonsZeroShort;

    PackedEvents::MC_Tracks fOutput_MC_NegKaons;
    PackedEvents::MC_Tracks fOutput_MC_PosKaons;
    PackedEvents::MC_Tracks fOutput_MC_PiMinus;
    PackedEvents::MC_Tracks fOutput_MC_PiPlus;

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
