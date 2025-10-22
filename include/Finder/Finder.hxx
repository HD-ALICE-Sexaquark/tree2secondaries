#ifndef T2S_FINDER_HXX
#define T2S_FINDER_HXX

#include <memory>
#include <utility>

#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

#include "App/Settings.hxx"
#include "DataFormats/ChannelA.hxx"
#include "DataFormats/ChannelD.hxx"
#include "DataFormats/Events.hxx"
#include "DataFormats/PackedEvents.hxx"
#include "KF/Extensions.hxx"
#include "MC/Particles.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries {

// Read packed events and find anti-sexaquarks reactions.
class Finder {
   public:
    Finder(const Finder &) = delete;
    Finder(Finder &&) = delete;
    Finder &operator=(const Finder &) = delete;
    Finder &operator=(Finder &&) = delete;
    ~Finder() = default;

    explicit Finder(Settings settings) : fSettings{std::move(settings)} {}

    [[nodiscard]] EReactionChannel GetReactionChannel() const { return fSettings.Channel; }

    bool Initialize();
    void ReadInputBranches();

    void ReadBranches_Events();
    void ReadBranches_Injected();
    void ReadBranches_V0s(EParticle pid, DF::Packed::V0s &df);
    void ReadBranches_Tracks(EParticle pid, DF::Packed::Tracks &df);
    void ReadBranches_LinkedV0s(EParticle pid, DF::Packed::LinkedV0s &df);
    void ReadBranches_LinkedTracks(EParticle pid, DF::Packed::LinkedTracks &df);

    bool PrepareOutputFile();
    bool PrepareOutputTree();
    void CreateCutFlowHistogram();

    void CreateOutputBranches(DF::ChannelA &df);
    void CreateOutputBranches(DF::ChannelD &df);
    void CreateOutputBranches(DF::MC_ChannelA &df);
    void CreateOutputBranches(DF::MC_ChannelD &df);
    void CreateOutputBranches() {
        switch (GetReactionChannel()) {
            // standard channels //
            case EReactionChannel::A:
            case EReactionChannel::AntiA:
                CreateOutputBranches(fOutput_ChannelA);
                if (IsMC()) CreateOutputBranches(fOutput_MC_ChannelA);
                break;
            case EReactionChannel::D:
            case EReactionChannel::AntiD:
                CreateOutputBranches(fOutput_ChannelD);
                if (IsMC()) CreateOutputBranches(fOutput_MC_ChannelD);
                break;
            default:
                break;
        }  // end of switch statement
    }

    bool Injected_PrepareOutputTree();
    void Injected_CreateOutputBranches();
    void Injected_FlattenAndStore();

    [[nodiscard]] int NumberEventsToRead() const {
        return fSettings.LimitToNEvents ? fSettings.LimitToNEvents : static_cast<int>(fInputChain_PackedEvents->GetEntries());
    }
    [[nodiscard]] bool IsMC() const { return fSettings.IsMC; }
    void GetEvent(int i_event) { fInputChain_PackedEvents->GetEntry(i_event); }

    void FindSexaquarks_ChannelA(bool anti_channel);
    void FindSexaquarks_ChannelD(bool anti_channel);
    void Find(EReactionChannel reaction_channel) {
        switch (reaction_channel) {
            case EReactionChannel::A:
            case EReactionChannel::AntiA:
                FindSexaquarks_ChannelA(reaction_channel == EReactionChannel::AntiA);
                break;
            case EReactionChannel::D:
            case EReactionChannel::AntiD:
                FindSexaquarks_ChannelD(reaction_channel == EReactionChannel::AntiD);
                break;
            default:
                return;
        }
    }

    [[nodiscard]] bool PassesCuts(const KF::ChannelA &sexa) const;
    [[nodiscard]] bool PassesCuts(const KF::ChannelD &sexa) const;

    void Store(const KF::ChannelA &sexa);
    void Store(const KF::ChannelD &sexa);
    void StoreMC(const MC::ChannelA &sexa);
    void StoreMC(const MC::ChannelD &sexa);

    void EndOfAnalysis();

   private:
    Settings fSettings;
    std::unique_ptr<TChain> fInputChain_PackedEvents;

    std::unique_ptr<TFile> fOutputFile;
    std::unique_ptr<TTree> fOutputTree;
    std::unique_ptr<TTree> fOutputTree_Injected;

    std::unique_ptr<TH1D> fCutFlowHist;

    // input //

    DF::Flat::Event fInput_Event;
    DF::SOV::Injected fInput_Injected;

    DF::Packed::V0s fInput_AntiLambdas;
    DF::Packed::V0s fInput_Lambdas;
    DF::Packed::V0s fInput_KaonsZeroShort;

    DF::Packed::Tracks fInput_NegKaons;
    DF::Packed::Tracks fInput_PosKaons;
    DF::Packed::Tracks fInput_PiMinus;
    DF::Packed::Tracks fInput_PiPlus;

    DF::Packed::LinkedV0s fInput_Linked_AntiLambdas;
    DF::Packed::LinkedV0s fInput_Linked_Lambdas;
    DF::Packed::LinkedV0s fInput_Linked_KaonsZeroShort;

    DF::Packed::LinkedTracks fInput_Linked_NegKaons;
    DF::Packed::LinkedTracks fInput_Linked_PosKaons;
    DF::Packed::LinkedTracks fInput_Linked_PiMinus;
    DF::Packed::LinkedTracks fInput_Linked_PiPlus;

    // output //

    DF::ChannelA fOutput_ChannelA;
    DF::ChannelD fOutput_ChannelD;

    DF::Flat::Injected fOutput_Injected;

    DF::MC_ChannelA fOutput_MC_ChannelA;
    DF::MC_ChannelD fOutput_MC_ChannelD;
};

}  // namespace Tree2Secondaries

#endif  // T2S_FINDER_HXX
