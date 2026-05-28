#pragma once

#include <memory>

#include <TFile.h>
#include <TH1.h>

#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleWriter.hxx>

#include "common/FL_ChannelA.hpp"
#include "common/FL_ChannelD.hpp"
#include "common/FL_ChannelH.hpp"
#include "common/FL_Event.hpp"
#include "common/FL_InjectedSexa.hpp"
#include "common/FL_Track.hpp"
#include "common/FL_V0.hpp"
#include "common/VC_InjectedSexa.hpp"
#include "common/VC_InjectedSexaView.hpp"
#include "common/VC_Track.hpp"
#include "common/VC_TrackView.hpp"
#include "common/VC_V0.hpp"
#include "common/VC_V0View.hpp"

#include "App/Settings.hxx"
#include "KalmanFitter/KalmanFitterChannelA.hxx"
#include "KalmanFitter/KalmanFitterChannelD.hxx"
#include "KalmanFitter/KalmanFitterChannelH.hxx"

class TTree;

namespace R2DS {

// Read Vector events and find anti-sexaquarks reactions.
class Finder {
   public:
    Finder(const Finder &) = delete;
    Finder(Finder &&) = delete;
    Finder &operator=(const Finder &) = delete;
    Finder &operator=(Finder &&) = delete;
    ~Finder() = default;

    explicit Finder(const Settings &settings) : fSettings{settings} {}

    bool Initialize();
    void ReadInputBranches();

    bool PrepareOutputFile();
    void PrepareOutputHistograms();

    void ProcessEvent();

    bool Injected_PrepareOutputTree();
    void ProcessInjected();

    [[nodiscard]] unsigned long NumberEventsToRead() const {
        unsigned long total_entries = fReader->GetNEntries();
        return 0 < fSettings.LimitToNEvents && fSettings.LimitToNEvents < total_entries ? fSettings.LimitToNEvents : total_entries;
    }

    void Find() {
        switch (fSettings.ReactionChannel.name) {
            case 'A': {
                FindSexaquarks_ChannelA(false);
                FindSexaquarks_ChannelA(true);
                break;
            }
            case 'D': {
                FindSexaquarks_ChannelD(false);
                FindSexaquarks_ChannelD(true);
                break;
            }
            case 'H': {
                FindSexaquarks_ChannelH(false);
                FindSexaquarks_ChannelH(true);
                break;
            }
            default:
                return;
        }
    }

    [[nodiscard]] bool FastCuts_ChannelA(const Seeder::Seed &seed_v0a, const Seeder::Seed &seed_v0b, TH1D *cut_flow_hist) const;
    [[nodiscard]] bool FastCuts_ChannelD(const Seeder::Seed &seed_ka, const Seeder::Seed &seed_v0, TH1D *cut_flow_hist) const;
    [[nodiscard]] bool FastCuts_ChannelH(const Seeder::Seed &seed_kaon1, const Seeder::Seed &seed_kaon2, TH1D *cut_flow_hist) const;

    [[nodiscard]] bool SlowCuts(const KalmanFitter::ChannelA &sexa, TH1D *cut_flow_hist) const;
    [[nodiscard]] bool SlowCuts(const KalmanFitter::ChannelD &sexa, TH1D *cut_flow_hist) const;
    [[nodiscard]] bool SlowCuts(const KalmanFitter::ChannelH &sexa, TH1D *cut_flow_hist) const;

    void EndOfAnalysis();

    std::unique_ptr<ROOT::RNTupleModel> fInput_Model;
    std::unique_ptr<ROOT::RNTupleReader> fReader;

    std::unique_ptr<ROOT::RNTupleModel> fOutput_Model;
    std::unique_ptr<ROOT::RNTupleWriter> fWriter;

    std::unique_ptr<ROOT::RNTupleModel> fOutput_Model_InjectedSexa;
    std::unique_ptr<ROOT::RNTupleWriter> fWriter_InjectedSexa;

   private:
    void Assign_Event(Flat::Sexaquark &out);
    void Assign_Candidate(const KalmanFitter::Particle &sexa, Flat::Sexaquark &out, bool anti_channel);
    void Assign(const Vector::V0View &v0, Flat::V0 &out);
    void Assign(const Vector::TrackView &track, Flat::Track &out, bool include_gm);
    void Assign_InjectedSexa(const std::optional<Vector::InjectedSexaView> &injected, Flat::InjectedSexa &out, bool embedded_to_rec);

    void FindSexaquarks_ChannelA(bool anti_channel);
    std::optional<Vector::InjectedSexaView> LinkInjectedSignal(const KalmanFitter::ChannelA &sexa);
    void Assign(const KalmanFitter::ChannelA &sexa, bool anti_channel);

    void FindSexaquarks_ChannelD(bool anti_channel);
    std::optional<Vector::InjectedSexaView> LinkInjectedSignal(const KalmanFitter::ChannelD &sexa);
    void Assign(const KalmanFitter::ChannelD &sexa, bool anti_channel);

    void FindSexaquarks_ChannelH(bool anti_channel);
    std::optional<Vector::InjectedSexaView> LinkInjectedSignal(const KalmanFitter::ChannelH &sexa);
    void Assign(const KalmanFitter::ChannelH &sexa, bool anti_channel);

    const Settings &fSettings;
    std::string fName_FoundRNT;

    std::unique_ptr<TFile> fOutput_File;

    std::unique_ptr<TH1D> fHist_EventCounter;

    std::unique_ptr<TH1D> fHist_CutFlow;
    std::unique_ptr<TH1D> fHist_CutFlow_AntiChannel;

    // input //

    Flat::Event fInput_Event;
    ROOT::Math::XYZPoint fPrimaryVertex;

    Vector::InjectedSexa fInput_InjectedSexa;

    Vector::V0 fInput_AntiLambda;
    Vector::V0 fInput_Lambda;
    Vector::V0 fInput_KaonZeroShort;

    Vector::Track fInput_NegKaon;
    Vector::Track fInput_PosKaon;

    // output //

    Flat::InjectedSexa fOutput_InjectedSexa;

    Flat::ChannelA fOutput_ChannelA;
    Flat::ChannelD fOutput_ChannelD;
    Flat::ChannelH fOutput_ChannelH;
};

}  // namespace R2DS
