#ifndef T2S_PACKAGER_HXX
#define T2S_PACKAGER_HXX

#include <memory>
#include <utility>

#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

#include "ALICE/ESD.hxx"
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

    EReactionChannel GetReactionChannel() const { return fSettings.Channel; }

    bool Initialize();
    void ConnectInputBranches();

    void ConnectBranches_Events();
    void ConnectBranches_Injected();
    void ConnectBranches_MC();
    void ConnectBranches_Tracks();

    bool PrepareOutputFile();
    bool PrepareOutputTree();
    void CreateOutputBranches();
    void CreateCutFlowHistograms();

    void CreateOutputBranches_Events();
    void CreateOutputBranches_Injected();
    void CreateOutputBranches_V0s(EParticle pid, PackedEvents::V0s &sov);
    void CreateOutputBranches_Tracks(EParticle pid, PackedEvents::Tracks &sov);
    void CreateOutputBranches_MC_V0s(EParticle pid, PackedEvents::MC_V0s &sov);
    void CreateOutputBranches_MC_Tracks(EParticle pid, PackedEvents::MC_Tracks &sov);

    int NumberEventsToRead() const { return fSettings.LimitToNEvents ? fSettings.LimitToNEvents : static_cast<int>(fEventsTree->GetEntries()); }
    bool IsMC() const { return fSettings.IsMC; }
    void GetEvent(int i_event) { fEventsTree->GetEntry(i_event); }

    int NumberMC() const { return static_cast<int>(fInput_MC.Px->size()); }
    int NumberInjected() const { return static_cast<int>(fInput_Injected.ReactionID->size()); }
    int NumberTracks() const { return static_cast<int>(fInput_Tracks.Px->size()); }

    void ProcessEvent();

    void Injected_GetSecondaryVertex();
    void Injected_Store();
    void ProcessInjected() {
        Injected_GetSecondaryVertex();
        Injected_Store();
    }

    void ProcessTracks();

    void PackTracks(EParticle pid);

    void FindV0s(EParticle pid);
    bool PassesCuts(const KF::V0 &v0, EParticle pid) const {
        switch (pid) {
            case EParticle::AntiLambda:
            case EParticle::Lambda:
                return PassesCuts_Lambda(v0);
            case EParticle::KaonZeroShort:
                return PassesCuts_KaonZeroShort(v0);
            default:
                return false;
        }
    }

    void EndOfEvent();
    void EndOfAnalysis();

   private:
    void Store(const KF::Track &track, PackedEvents::Tracks &sov);
    void StoreMC(const MC::Track &mc_track, PackedEvents::MC_Tracks &sov);

    bool PassesCuts_Lambda(const KF::V0 &v0) const;
    bool PassesCuts_KaonZeroShort(const KF::V0 &v0) const;
    void Store(const KF::V0 &v0, PackedEvents::V0s &sov);
    void StoreMC(const MC::V0 &mc_v0, PackedEvents::MC_V0s &sov);

    void Store(const ALICE::V0 &v0, PackedEvents::V0s &sov);

    Settings fSettings;
    std::unique_ptr<TChain> fEventsTree;

    std::unique_ptr<TFile> fOutputFile;
    std::unique_ptr<TTree> fOutputTree;

    std::unique_ptr<TH1D> fCutFlowHist_AntiLambdas;
    std::unique_ptr<TH1D> fCutFlowHist_Lambdas;
    std::unique_ptr<TH1D> fCutFlowHist_KaonsZeroShort;

    // input branches //

    Events::Event fInput_Event;
    Events::Injected fInput_Injected;

    Events::MC fInput_MC;
    Events::Tracks fInput_Tracks;

    // temporary containers, cleaned after event loop //

    std::vector<float> fVec_SV_X;
    std::vector<float> fVec_SV_Y;
    std::vector<float> fVec_SV_Z;

    std::vector<int> fVec_AntiProtons;
    std::vector<int> fVec_Protons;
    std::vector<int> fVec_NegKaons;
    std::vector<int> fVec_PosKaons;
    std::vector<int> fVec_PiMinus;
    std::vector<int> fVec_PiPlus;

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
