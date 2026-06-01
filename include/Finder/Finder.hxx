#pragma once

#include <memory>

#include <TFile.h>
#include <TH1.h>

#include <Math/Point3D.h>
#include <ROOT/RNTupleReader.hxx>

#include "common/Framework.hpp"
#include "common/Schema_FoundSexa.hpp"
#include "common/Schema_PackedEvents.hpp"

#include "App/Settings.hxx"

// forward declarations //
// clang-format off
namespace POD {
    struct ChannelA;
    struct ChannelD;
    struct ChannelH;
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
          // input
          fInput_File{std::make_unique<TFile>(fSettings.PathInputFile.c_str(), "READ")},
          fReader{E2R::Name_OutputRNT, *fInput_File},
          fInput{fReader.Data()},
          // output
          fOutput_File{std::make_unique<TFile>(fSettings.PathOutputFile.c_str(), "RECREATE")},
          fWriter{R2DS::Name_FoundSexaRNT, *fOutput_File},
          fOutput{fWriter.Data()} {

        PrepareOutputHistograms();

        Logger::Info(__FUNCTION__, "Finder initialized successfully.");
    }

    void PrepareOutputHistograms();

    void Load(ROOT::NTupleSize_t entry_id) { fReader.Load(entry_id); }
    [[nodiscard]] ROOT::NTupleSize_t NumberEventsToRead() {
        auto total = fReader.Iter()->GetNEntries();
        return fSettings.LimitToNEvents.has_value() ? std::min(fSettings.LimitToNEvents.value(), total) : total;
    }

    void ProcessEvent();
    void ProcessInjected();

    void Find() {
        FindSexaquarks_ChannelA(false);
        FindSexaquarks_ChannelA(true);
        FindSexaquarks_ChannelD(false);
        FindSexaquarks_ChannelD(true);
        FindSexaquarks_ChannelH(false);
        FindSexaquarks_ChannelH(true);
    }

    void EndOfEvent();
    void EndOfAnalysis();

   private:
    // channel A //

    void FindSexaquarks_ChannelA(bool control_channel);

    [[nodiscard]] bool FastCuts_ChannelA(const Seeder::Seed &seed_v0a, const Seeder::Seed &seed_v0b, TH1D *cut_flow_hist) const;
    [[nodiscard]] bool SlowCuts(const KF::ChannelA &kf_sexa, TH1D *cut_flow_hist) const;

    void BuildMcInfo(POD::ChannelA &new_sexa, bool control_channel);
    void BuildRecInfo(POD::ChannelA &new_sexa, const KF::ChannelA &kf_sexa, bool control_channel);

    // channel D //

    void FindSexaquarks_ChannelD(bool control_channel);

    [[nodiscard]] bool FastCuts_ChannelD(const Seeder::Seed &seed_ka, const Seeder::Seed &seed_v0, TH1D *cut_flow_hist) const;
    [[nodiscard]] bool SlowCuts(const KF::ChannelD &kf_sexa, TH1D *cut_flow_hist) const;

    void BuildMcInfo(POD::ChannelD &new_sexa, bool control_channel);
    void BuildRecInfo(POD::ChannelD &new_sexa, const KF::ChannelD &kf_sexa, bool control_channel);

    // channel H //

    void FindSexaquarks_ChannelH(bool control_channel);

    [[nodiscard]] bool FastCuts_ChannelH(const Seeder::Seed &seed_kaon1, const Seeder::Seed &seed_kaon2, TH1D *cut_flow_hist) const;
    [[nodiscard]] bool SlowCuts(const KF::ChannelH &kf_sexa, TH1D *cut_flow_hist) const;

    void BuildMcInfo(POD::ChannelH &new_sexa, bool control_channel);
    void BuildRecInfo(POD::ChannelH &new_sexa, const KF::ChannelH &kf_sexa, bool control_channel);

    // member variables //

    const Settings &fSettings;

    std::unique_ptr<TH1D> fHist_EventCounter;

    std::unique_ptr<TH1D> fHist_CutFlow;
    std::unique_ptr<TH1D> fHist_CutFlow_AntiChannel;

    // input //

    std::unique_ptr<TFile> fInput_File;
    Framework::Reader<Schema::PackedEvents> fReader;
    Schema::PackedEvents &fInput;
    // -- cached
    ROOT::Math::XYZPoint fPrimaryVertex;

    // output //

    std::unique_ptr<TFile> fOutput_File;
    Framework::Writer<Schema::FoundSexa> fWriter;
    Schema::FoundSexa &fOutput;
};

}  // namespace R2DS
