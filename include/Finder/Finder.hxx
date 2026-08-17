#pragma once

#include <exception>
#include <memory>
#include <string>
#include <string_view>
#include <utility>

#include <TFile.h>
#include <TH1.h>

#include <Math/Point3D.h>

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
namespace POD { struct InjectedSexa; struct Track; struct V0; }
namespace Cached { struct Sexaquark; }

namespace T2DS {

namespace Seeder { struct PCA; }
namespace KF { struct Particle; struct FitResult; }
// clang-format on

class Finder {
    // PENDING for author: missing fit constraints + cuts enums+structs

   public:
    Finder() = delete;
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

        fWriter = std::make_unique<Framework::Writer>(std::move(model), fName_FoundRNT, *fOutput_File);

        PrepareOutputHistograms();

        Logger::Info(__FUNCTION__, "Finder initialized successfully.");
    }

    [[nodiscard]] bool OpenInput(std::string_view path) {
        fReader.reset();
        try {
            fReader = std::make_unique<Framework::Reader>(fInput.CreateModel(fSettings.IsMC), T2DS::Name_PackedRNT, path);
        } catch (const std::exception &exc) {
            Logger::Error(__FUNCTION__, "Couldn't read {} ({}) -- skipping it.", path, exc.what());
            return false;
        }
        return true;
    }

    void PrepareOutputHistograms();

    void Load(ROOT::NTupleSize_t entry_id) { fReader->Load(entry_id); }
    [[nodiscard]] unsigned long NumberEventsToRead() { return fReader->Iter()->GetNEntries(); }

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
    bool EndOfAnalysis();

   private:
    // mc sexaquark //
    POD::Linked::InjectedSexa BuildMcSexaquark(const POD::Extended::McParticle &mc_dau1, const POD::Extended::McParticle &mc_dau2);

    // channel A //
    void FindSexaquarks_ChannelA(bool is_bkg_channel);
    [[nodiscard]] bool FastCuts_ChannelA(const Seeder::PCA &pca_v0a, const Seeder::PCA &pca_v0b, TH1D *hist_cut_flow) const;
    [[nodiscard]] bool SlowCuts_ChannelA(const Cached::Sexaquark &c_sexa, TH1D *hist_cut_flow) const;
    POD::Sexaquark Create_ChannelA(const KF::FitResult &fit, const Seeder::PCA &pca_v0a, const Seeder::PCA &pca_v0b, bool is_bkg_channel);

    // channel D //
    void FindSexaquarks_ChannelD(bool is_bkg_channel);
    [[nodiscard]] bool FastCuts_ChannelD(const Seeder::PCA &pca_v0, const Seeder::PCA &pca_ka, TH1D *hist_cut_flow) const;
    [[nodiscard]] bool SlowCuts_ChannelD(const Cached::Sexaquark &c_sexa, TH1D *hist_cut_flow) const;
    POD::Sexaquark Create_ChannelD(const KF::FitResult &fit, const Seeder::PCA &pca_v0, const Seeder::PCA &pca_ka, bool is_bkg_channel);

    // channel H //
    void FindSexaquarks_ChannelH(bool is_bkg_channel);
    [[nodiscard]] bool FastCuts_ChannelH(const Seeder::PCA &pca_kaon1, const Seeder::PCA &pca_kaon2, TH1D *hist_cut_flow) const;
    [[nodiscard]] bool SlowCuts_ChannelH(const Cached::Sexaquark &c_sexa, TH1D *hist_cut_flow) const;
    POD::Sexaquark Create_ChannelH(const KF::FitResult &fit, const Seeder::PCA &pca_kaon1, const Seeder::PCA &pca_kaon2, bool is_bkg_channel);

    // member variables //

    const Settings &fSettings;
    DB::ReactionChannels::Definition fReactionChannel;
    std::string fName_FoundRNT;

    // input //

    Schema::PackedEvents fInput;
    std::unique_ptr<Framework::Reader> fReader;
    // -- cached
    ROOT::Math::XYZPoint fPrimaryVertex;

    // output //

    std::unique_ptr<TFile> fOutput_File;  // single file, kept alive across every input file, if multiple
    Schema::FoundChannelA fOutput_ChannelA;
    Schema::FoundChannelD fOutput_ChannelD;
    Schema::FoundChannelH fOutput_ChannelH;
    Schema::FoundSexaquark *fOutput_Base;
    std::unique_ptr<Framework::Writer> fWriter;

    // histograms
    // -- event counter
    std::unique_ptr<TH1D> fHist_EventCounter;
    // -- cut flow for anti-sexaquarks + bkg. sexaquarks
    std::unique_ptr<TH1D> fHist_CutFlow_AntiSexaquark;
    std::unique_ptr<TH1D> fHist_CutFlow_BkgCandidates;
};

}  // namespace T2DS
