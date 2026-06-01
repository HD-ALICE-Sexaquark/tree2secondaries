#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include <TFile.h>
#include <TH1.h>

#include <Math/Point3D.h>

#include <ROOT/RNTupleReader.hxx>

#include "App/Logger.hxx"
#include "common/DB_Particles.hpp"
#include "common/Framework.hpp"
#include "common/Schema_Events.hpp"
#include "common/Schema_PackedEvents.hpp"

#include "App/Settings.hxx"

// forward declarations //
// clang-format off
namespace POD {
    struct Track;
    struct V0;
}
namespace KF { struct V0; }

namespace R2DS {

namespace Seeder{ struct Seed; }
// clang-format on

// Pack secondary V0s and tracks.
class Packager {
   public:
    Packager(const Packager &) = delete;
    Packager(Packager &&) = delete;
    Packager &operator=(const Packager &) = delete;
    Packager &operator=(Packager &&) = delete;
    ~Packager() = default;

    explicit Packager(const Settings &settings)
        : fSettings{settings},
          // input
          fInput_File{std::make_unique<TFile>(fSettings.PathInputFile.c_str(), "READ")},
          fReader{E2R::Name_OutputRNT, *fInput_File},
          fInput{fReader.Data()},
          // output
          fOutput_File{std::make_unique<TFile>(fSettings.PathOutputFile.c_str(), "RECREATE")},
          fWriter{R2DS::Name_PackedRNT, *fOutput_File},
          fOutput{fWriter.Data()} {

        PrepareOutputHistograms();

        Logger::Info(__FUNCTION__, "Packager initialized successfully.");
    }

    [[nodiscard]] bool CheckArguments() const {
        // -- if data kind is real data, all good
        if (!fInput.McParticle.has_value()) return true;
        // -- if data kind is MC, it requires
        bool state = fSettings.ReactionChannel.has_value() && fSettings.SexaquarkMass.has_value();
        if (!state) Logger::Error(__FUNCTION__, "\"Packager\" is reading MC, but \"--channel <channel> --mass <mass>\" are missing.");
        return state;
    }

    void PrepareOutputHistograms();

    void Load(ROOT::NTupleSize_t entry_id) { fReader.Load(entry_id); }
    [[nodiscard]] ROOT::NTupleSize_t NumberEventsToRead() {
        auto total = fReader.Iter()->GetNEntries();
        return fSettings.LimitToNEvents.has_value() ? std::min(fSettings.LimitToNEvents.value(), total) : total;
    }

    void ProcessEvent();
    void ProcessInjected();
    void ProcessTracks();

    void Pack() {
        FindV0s(DB::Particles::Particle("AntiLambda"));
        FindV0s(DB::Particles::Particle("Lambda"));
        FindV0s(DB::Particles::Particle("KaonZeroShort"));
        PackTracks(DB::Particles::Particle("NegKaon"));
        PackTracks(DB::Particles::Particle("PosKaon"));
    }

    [[nodiscard]] bool FastCuts(const Seeder::Seed &pca_neg, const Seeder::Seed &pca_pos, const DB::Particles::Definition &pid) const {
        switch (pid.pdg_code) {
            case DB::Particles::Particle("AntiLambda").pdg_code: {
                return FastCuts_Lambda(pca_neg, pca_pos, fHist_CutFlow_AntiLambda.get());
            }
            case DB::Particles::Particle("Lambda").pdg_code: {
                return FastCuts_Lambda(pca_neg, pca_pos, fHist_CutFlow_Lambda.get());
            }
            case DB::Particles::Particle("KaonZeroShort").pdg_code: {
                return FastCuts_KaonZeroShort(pca_neg, pca_pos, fHist_CutFlow_KaonZeroShort.get());
            }
            default:
                return false;
        }
    }

    [[nodiscard]] bool SlowCuts(const KF::V0 &v0, const DB::Particles::Definition &pid) const {
        switch (pid.pdg_code) {
            case DB::Particles::Particle("AntiLambda").pdg_code: {
                return SlowCuts_Lambda(v0, fHist_CutFlow_AntiLambda.get());
            }
            case DB::Particles::Particle("Lambda").pdg_code: {
                return SlowCuts_Lambda(v0, fHist_CutFlow_Lambda.get());
            }
            case DB::Particles::Particle("KaonZeroShort").pdg_code: {
                return SlowCuts_KaonZeroShort(v0, fHist_CutFlow_KaonZeroShort.get());
            }
            default:
                return false;
        }
    }

    void EndOfEvent();
    void EndOfAnalysis();

   private:
    // tracks //

    void PackTracks(const DB::Particles::Definition &pid);

    bool PassesProtonCuts(const POD::Track &track, TH1D *cut_flow_hist) const;
    bool PassesKaonCuts(const POD::Track &track, TH1D *cut_flow_hist) const;
    bool PassesPionCuts(const POD::Track &track, TH1D *cut_flow_hist) const;

    void BuildMcInfo(POD::Track &new_track, int pdg_code_hypothesis, bool include_gm);

    // V0s //

    void FindV0s(const DB::Particles::Definition &pid);

    bool FastCuts_Lambda(const Seeder::Seed &seed_neg, const Seeder::Seed &seed_pos, TH1D *cut_flow_hist) const;
    bool FastCuts_KaonZeroShort(const Seeder::Seed &seed_neg, const Seeder::Seed &seed_pos, TH1D *cut_flow_hist) const;

    bool SlowCuts_Lambda(const KF::V0 &kf_v0, TH1D *cut_flow_hist) const;
    bool SlowCuts_KaonZeroShort(const KF::V0 &kf_v0, TH1D *cut_flow_hist) const;

    void BuildMcInfo(POD::V0 &new_v0, const DB::Particles::Definition &pid_hypothesis);
    void BuildRecInfo(POD::V0 &new_v0, const KF::V0 &kf_v0);

    // member variables //

    const Settings &fSettings;

    std::unique_ptr<TH1D> fHist_EventCounter;  // event counter

    std::unique_ptr<TH1D> fHist_CutFlow_AntiProton;
    std::unique_ptr<TH1D> fHist_CutFlow_Proton;
    std::unique_ptr<TH1D> fHist_CutFlow_NegKaon;
    std::unique_ptr<TH1D> fHist_CutFlow_PosKaon;
    std::unique_ptr<TH1D> fHist_CutFlow_PiMinus;
    std::unique_ptr<TH1D> fHist_CutFlow_PiPlus;

    std::unique_ptr<TH1D> fHist_CutFlow_AntiLambda;
    std::unique_ptr<TH1D> fHist_CutFlow_Lambda;
    std::unique_ptr<TH1D> fHist_CutFlow_KaonZeroShort;

    // temporary entries containers, cleaned after event loop //

    std::vector<std::size_t> fEntries_AntiProton;
    std::vector<std::size_t> fEntries_Proton;
    std::vector<std::size_t> fEntries_NegKaon;
    std::vector<std::size_t> fEntries_PosKaon;
    std::vector<std::size_t> fEntries_PiMinus;
    std::vector<std::size_t> fEntries_PiPlus;

    // input //

    std::unique_ptr<TFile> fInput_File;
    Framework::Reader<Schema::Events> fReader;
    Schema::Events &fInput;
    // -- cached
    ROOT::Math::XYZPoint fPrimaryVertex;

    // output //

    std::unique_ptr<TFile> fOutput_File;
    Framework::Writer<Schema::PackedEvents> fWriter;
    Schema::PackedEvents &fOutput;
};

}  // namespace R2DS
