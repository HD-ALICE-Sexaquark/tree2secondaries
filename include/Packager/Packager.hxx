#ifndef T2S_PACKAGER_HXX
#define T2S_PACKAGER_HXX

#include <memory>

#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

#include "App/Settings.hxx"
#include "DataFormats/Events.hxx"
#include "DataFormats/Injected.hxx"
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

    [[nodiscard]] EReactionChannel GetReactionChannel() const { return fSettings.Channel; }

    bool Initialize();
    void ReadInputBranches();

    void ReadBranches_Events();
    void ReadBranches_Injected();
    void ReadBranches_MC();
    void ReadBranches_Tracks();

    bool PrepareOutputFile();
    bool PrepareOutputTree();
    void CreateOutputBranches();
    void CreateCutFlowHistograms();

    void CreateOutputBranches_Events();
    void CreateOutputBranches_Injected();
    void CreateOutputBranches_V0s(EParticle pid, DF::Packed::V0s &df);
    void CreateOutputBranches_Tracks(EParticle pid, DF::Packed::Tracks &df);
    void CreateOutputBranches_LinkedV0s(EParticle pid, DF::Packed::LinkedV0s &df);
    void CreateOutputBranches_LinkedTracks(EParticle pid, DF::Packed::LinkedTracks &df);

    [[nodiscard]] int NumberEventsToRead() const {
        return fSettings.LimitToNEvents ? fSettings.LimitToNEvents : static_cast<int>(fInputChain_Events->GetEntries());
    }
    [[nodiscard]] bool IsMC() const { return fSettings.IsMC; }
    void GetEvent(int i_event) { fInputChain_Events->GetEntry(i_event); }

    [[nodiscard]] int NumberMC() const { return static_cast<int>(fInput_MC.Px->size()); }
    [[nodiscard]] int NumberInjected() const { return static_cast<int>(fInput_Injected.ReactionID->size()); }
    [[nodiscard]] int NumberTracks() const { return static_cast<int>(fInput_Tracks.Px->size()); }

    void ProcessEvent();
    void ProcessMCEvent();

    void Injected_GetSecondaryVertex();
    void Injected_Store();
    void ProcessInjected() {
        Injected_GetSecondaryVertex();
        Injected_Store();
    }

    void ProcessTracks();

    void PackTracks(EParticle pid);

    void FindV0s(EParticle pid);
    [[nodiscard]] bool PassesCuts(const KF::V0 &v0, EParticle pid) const {
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
    void Store(const KF::Track &track, DF::Packed::Tracks &df);
    void StoreMC(const MC::Track &mc_track, DF::Packed::LinkedTracks &df);

    [[nodiscard]] bool PassesCuts_Lambda(const KF::V0 &v0) const;
    [[nodiscard]] bool PassesCuts_KaonZeroShort(const KF::V0 &v0) const;
    void Store(const KF::V0 &v0, DF::Packed::V0s &df);
    void StoreMC(const MC::V0 &mc_v0, DF::Packed::LinkedV0s &df);

    Settings fSettings;
    std::unique_ptr<TChain> fInputChain_Events;

    std::unique_ptr<TFile> fOutputFile;
    std::unique_ptr<TTree> fOutputTree;

    std::unique_ptr<TH1D> fCutFlowHist_AntiLambdas;
    std::unique_ptr<TH1D> fCutFlowHist_Lambdas;
    std::unique_ptr<TH1D> fCutFlowHist_KaonsZeroShort;

    // input branches //

    DF::Flat::Event fInput_Event;
    DF::SOV::Injected fInput_Injected;

    DF::SOV::MC_Particles fInput_MC;
    DF::SOV::Tracks fInput_Tracks;

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

    DF::Flat::Event fOutput_Event;
    DF::SOV::Injected fOutput_Injected;

    DF::Packed::LinkedV0s fOutput_Linked_AntiLambdas;
    DF::Packed::LinkedV0s fOutput_Linked_Lambdas;
    DF::Packed::LinkedV0s fOutput_Linked_KaonsZeroShort;

    DF::Packed::LinkedTracks fOutput_Linked_NegKaons;
    DF::Packed::LinkedTracks fOutput_Linked_PosKaons;
    DF::Packed::LinkedTracks fOutput_Linked_PiMinus;
    DF::Packed::LinkedTracks fOutput_Linked_PiPlus;

    DF::Packed::V0s fOutput_AntiLambdas;
    DF::Packed::V0s fOutput_Lambdas;
    DF::Packed::V0s fOutput_KaonsZeroShort;

    DF::Packed::Tracks fOutput_NegKaons;
    DF::Packed::Tracks fOutput_PosKaons;
    DF::Packed::Tracks fOutput_PiMinus;
    DF::Packed::Tracks fOutput_PiPlus;
};

}  // namespace Tree2Secondaries

#endif  // T2S_PACKAGER_HXX
