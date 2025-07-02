#ifndef T2S_FINDER_HXX
#define T2S_FINDER_HXX

#include <memory>
#include <utility>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include "App/Settings.hxx"
#include "DataFormats/Events.hxx"
#include "DataFormats/Found.hxx"
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

    ReactionChannel GetReactionChannel() const { return fSettings.Channel; }

    bool Initialize();
    void ConnectInputBranches();

    void ConnectBranches_Events();
    void ConnectBranches_Injected();
    void ConnectBranches_V0s(const std::string &name_v0, PackedEvents::V0s &vec_v0s);
    void ConnectBranches_Tracks(const std::string &name_part, PackedEvents::Tracks &vec_tracks);
    void ConnectBranches_MC_V0s(const std::string &name_v0, PackedEvents::MC_V0s &sov);
    void ConnectBranches_MC_Tracks(const std::string &name_part, PackedEvents::MC_Tracks &sov);

    bool PrepareOutputFile();
    bool PrepareOutputTree();

    void CreateOutputBranches(Found::ChannelA &out_branches);
    void CreateOutputBranches(Found::ChannelD &out_branches);
    void CreateOutputBranches(Found::ChannelE &out_branches);
    void CreateOutputBranches(Found::ChannelH &out_branches);
    void CreateOutputBranches(Found::MC_ChannelA &sov);
    void CreateOutputBranches(Found::MC_ChannelD &sov);
    void CreateOutputBranches() {
        switch (GetReactionChannel()) {
            // standard channels //
            case ReactionChannel::A:
            case ReactionChannel::AntiA:
                CreateOutputBranches(fOutput_ChannelA);
                if (IsMC()) CreateOutputBranches(fOutput_MC_ChannelA);
                break;
            case ReactionChannel::D:
            case ReactionChannel::AntiD:
                CreateOutputBranches(fOutput_ChannelD);
                if (IsMC()) CreateOutputBranches(fOutput_MC_ChannelD);
                break;
            case ReactionChannel::E:
            case ReactionChannel::AntiE:
                CreateOutputBranches(fOutput_ChannelE);
                break;
            case ReactionChannel::H:
            case ReactionChannel::AntiH:
                CreateOutputBranches(fOutput_ChannelH);
                break;
            default:
                break;
        }  // end of switch statement
    }

    bool Injected_PrepareOutputTree();
    void Injected_CreateOutputBranches();
    void Injected_FlattenAndStore();

    int NumberEventsToRead() const {
        return fSettings.LimitToNEvents ? fSettings.LimitToNEvents : static_cast<int>(fTree_PackedEvents->GetEntries());
    }
    bool IsMC() const { return fSettings.IsMC; }
    bool IsSignalMC() const { return fSettings.IsSignalMC; }
    void GetEvent(int i_event) { fTree_PackedEvents->GetEntry(i_event); }

    void FindSexaquarks_ChannelA(bool anti_channel);
    void FindSexaquarks_ChannelD(bool anti_channel);
    void FindSexaquarks_ChannelE(bool anti_channel);
    void FindSexaquarks_ChannelH(bool anti_channel);
    void Find(ReactionChannel reaction_channel) {
        switch (reaction_channel) {
            case ReactionChannel::A:
            case ReactionChannel::AntiA:
                FindSexaquarks_ChannelA(reaction_channel == ReactionChannel::AntiA);
                break;
            case ReactionChannel::D:
            case ReactionChannel::AntiD:
                FindSexaquarks_ChannelD(reaction_channel == ReactionChannel::AntiD);
                break;
            case ReactionChannel::E:
            case ReactionChannel::AntiE:
                FindSexaquarks_ChannelE(reaction_channel == ReactionChannel::AntiE);
                break;
            case ReactionChannel::H:
            case ReactionChannel::AntiH:
                FindSexaquarks_ChannelH(reaction_channel == ReactionChannel::AntiH);
                break;
            default:
                return;
        }
    }

    bool PassesCuts(const KF::ChannelA &sexa) const;
    bool PassesCuts(const KF::ChannelD &sexa) const;
    bool PassesCuts(const KF::ChannelE &sexa) const;
    bool PassesCuts(const KF::ChannelH &sexa) const;

    void Store(const KF::ChannelA &sexa);
    void Store(const KF::ChannelD &sexa);
    void Store(const KF::ChannelE &sexa);
    void Store(const KF::ChannelH &sexa);
    void StoreMC(const MC::ChannelA &sexa);
    void StoreMC(const MC::ChannelD &sexa);

    void EndOfAnalysis();

   private:
    Settings fSettings;
    std::unique_ptr<TChain> fTree_PackedEvents;

    std::unique_ptr<TFile> fOutputFile;
    std::unique_ptr<TTree> fOutputTree;
    std::unique_ptr<TTree> fOutputTree_Injected;

    // input structs //

    Events::Event fInput_Event;
    Events::Injected fInput_Injected;

    PackedEvents::V0s fPacked_AntiLambdas;
    PackedEvents::V0s fPacked_Lambdas;
    PackedEvents::V0s fPacked_KaonsZeroShort;

    PackedEvents::Tracks fPacked_NegKaons;
    PackedEvents::Tracks fPacked_PosKaons;
    PackedEvents::Tracks fPacked_PiMinus;
    PackedEvents::Tracks fPacked_PiPlus;

    PackedEvents::MC_V0s fPacked_MC_AntiLambdas;
    PackedEvents::MC_V0s fPacked_MC_Lambdas;
    PackedEvents::MC_V0s fPacked_MC_KaonsZeroShort;

    PackedEvents::MC_Tracks fPacked_MC_NegKaons;
    PackedEvents::MC_Tracks fPacked_MC_PosKaons;
    PackedEvents::MC_Tracks fPacked_MC_PiMinus;
    PackedEvents::MC_Tracks fPacked_MC_PiPlus;

    // output structs //

    Found::ChannelA fOutput_ChannelA;
    Found::ChannelD fOutput_ChannelD;
    Found::ChannelE fOutput_ChannelE;
    Found::ChannelH fOutput_ChannelH;

    Found::MC_ChannelA fOutput_MC_ChannelA;
    Found::MC_ChannelD fOutput_MC_ChannelD;

    Found::Injected fOutput_Injected;
};

}  // namespace Tree2Secondaries

#endif  // T2S_FINDER_HXX
