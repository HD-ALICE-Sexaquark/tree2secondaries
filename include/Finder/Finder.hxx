#pragma once

#include <memory>
#include <utility>

#include <TFile.h>
#include <TH1.h>

#include <Math/Point3D.h>
#include <ROOT/RNTupleReader.hxx>

#include "common/Constants.hpp"
#include "common/DB_ReactionChannels.hpp"
#include "common/Framework.hpp"
#include "common/POD_InjectedSexa.hpp"
#include "common/Schema_FoundChannelA.hpp"
#include "common/Schema_FoundChannelD.hpp"
#include "common/Schema_FoundChannelH.hpp"
#include "common/Schema_FoundSexaquark.hpp"
#include "common/Schema_PackedEvents.hpp"

#include "App/Settings.hxx"

// forward declarations //
// clang-format off
namespace POD {
    struct InjectedSexa;
    struct Track;
    struct V0;
}
namespace KF {
    struct ChannelA;
    struct ChannelD;
    struct ChannelH;
}

namespace R2DS {

namespace Seeder { struct Seed; }
// clang-format on

class Finder {
   public:
    Finder(const Finder &) = delete;
    Finder(Finder &&) = delete;
    Finder &operator=(const Finder &) = delete;
    Finder &operator=(Finder &&) = delete;
    ~Finder() = default;

    explicit Finder(const Settings &settings)
        : fSettings{settings},
          fReactionChannel{DB::ReactionChannels::FindReactionChannel(fSettings.ReactionChannel)},
          fName_FoundRNT{std::format("FoundChannel{}", fReactionChannel.name)},
          // input
          fInput_File{std::make_unique<TFile>(fSettings.PathInputFile.c_str(), "READ")},
          fInput{},
          fReader{nullptr},
          // output
          fOutput_File{std::make_unique<TFile>(fSettings.PathOutputFile.c_str(), "RECREATE")},
          fOutput_ChannelA{},
          fOutput_ChannelD{},
          fOutput_ChannelH{},
          fOutput_Base{},
          fWriter{nullptr} {

        // further reaction channel rules //
        Framework::Model model;
        if (fReactionChannel.name == 'A') {
            model = fOutput_ChannelA.CreateModel(fSettings.IsMC);
            fOutput_Base = &fOutput_ChannelA;
        } else if (fReactionChannel.name == 'D') {
            model = fOutput_ChannelD.CreateModel(fSettings.IsMC);
            fOutput_Base = &fOutput_ChannelD;
        } else if (fReactionChannel.name == 'H') {
            model = fOutput_ChannelH.CreateModel(fSettings.IsMC);
            fOutput_Base = &fOutput_ChannelH;
        } else {
            Logger::Error(__FUNCTION__, "Invalid reaction channel.");  // NOTE: `Parser` protects against this
            return;
        }

        fReader = std::make_unique<Framework::Reader>(fInput.CreateModel(fSettings.IsMC), R2DS::Name_PackedRNT, *fInput_File);
        fWriter = std::make_unique<Framework::Writer>(std::move(model), fName_FoundRNT, *fOutput_File);

        PrepareOutputHistograms();

        Logger::Info(__FUNCTION__, "Finder initialized successfully.");
    }

    void PrepareOutputHistograms();

    void Load(ROOT::NTupleSize_t entry_id) { fReader->Load(entry_id); }
    [[nodiscard]] ROOT::NTupleSize_t NumberEventsToRead() {
        auto total = fReader->Iter()->GetNEntries();
        return fSettings.LimitToNEvents.has_value() ? std::min(fSettings.LimitToNEvents.value(), total) : total;
    }

    void ProcessEvent();
    void ProcessInjected();

    void Find() {
        if (fReactionChannel.name == 'A') {
            FindSexaquarks_ChannelA(false);
            FindSexaquarks_ChannelA(true);
            return;
        } else if (fReactionChannel.name == 'D') {
            FindSexaquarks_ChannelD(false);
            FindSexaquarks_ChannelD(true);
            return;
        } else if (fReactionChannel.name == 'H') {
            FindSexaquarks_ChannelH(false);
            FindSexaquarks_ChannelH(true);
            return;
        } else {
            return;
        }
    }

    void EndOfEvent();
    void EndOfAnalysis();

   private:
    // mc sexaquark //
    POD::Linked::InjectedSexa BuildMcSexaquark(const POD::Extended::McParticle &mc_dau1, const POD::Extended::McParticle &mc_dau2);

    // channel A //
    void FindSexaquarks_ChannelA(bool anti_channel);
    [[nodiscard]] bool FastCuts_ChannelA(const Seeder::Seed &seed_v0a, const Seeder::Seed &seed_v0b, TH1D *cut_flow_hist) const;
    [[nodiscard]] bool SlowCuts_ChannelA(const KF::ChannelA &kf_sexa, TH1D *cut_flow_hist) const;
    POD::Sexaquark Create_ChannelA(const KF::ChannelA &kf_sexa, bool anti_channel);

    // channel D //
    void FindSexaquarks_ChannelD(bool anti_channel);
    [[nodiscard]] bool FastCuts_ChannelD(const Seeder::Seed &seed_ka, const Seeder::Seed &seed_v0, TH1D *cut_flow_hist) const;
    [[nodiscard]] bool SlowCuts_ChannelD(const KF::ChannelD &kf_sexa, TH1D *cut_flow_hist) const;
    POD::Sexaquark Create_ChannelD(const KF::ChannelD &kf_sexa, bool anti_channel);

    // channel H //
    void FindSexaquarks_ChannelH(bool anti_channel);
    [[nodiscard]] bool FastCuts_ChannelH(const Seeder::Seed &seed_kaon1, const Seeder::Seed &seed_kaon2, TH1D *cut_flow_hist) const;
    [[nodiscard]] bool SlowCuts_ChannelH(const KF::ChannelH &kf_sexa, TH1D *cut_flow_hist) const;
    POD::Sexaquark Create_ChannelH(const KF::ChannelH &kf_sexa, bool anti_channel);

    // member variables //

    const Settings &fSettings;
    DB::ReactionChannels::Definition fReactionChannel;
    std::string fName_FoundRNT;

    // input //

    std::unique_ptr<TFile> fInput_File;
    Schema::PackedEvents fInput;
    std::unique_ptr<Framework::Reader> fReader;
    // -- cached
    ROOT::Math::XYZPoint fPrimaryVertex;

    // output //

    std::unique_ptr<TFile> fOutput_File;
    Schema::FoundChannelA fOutput_ChannelA;
    Schema::FoundChannelD fOutput_ChannelD;
    Schema::FoundChannelH fOutput_ChannelH;
    Schema::FoundSexaquark *fOutput_Base;  // NOTE: already initialized in constructor

    std::unique_ptr<Framework::Writer> fWriter;
    std::unique_ptr<TH1D> fHist_EventCounter;

    std::unique_ptr<TH1D> fHist_CutFlow;
    std::unique_ptr<TH1D> fHist_CutFlow_AntiChannel;
};

}  // namespace R2DS
