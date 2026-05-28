#include "Verifier/Verifier.hxx"

#include <filesystem>
#include <memory>

#include "common/Constants.hpp"
#include "common/DB_Particles.hpp"
#include "common/R2DS_Cuts.hpp"
#include "common/VC_McParticleView.hpp"
#include "common/VC_OnTheFlyLambdaView.hpp"

#include "KalmanFitter/KF_LambdaPair.hxx"
#include "Seeder/SeederLineLine.hxx"

namespace R2DS {

namespace KF = KalmanFitter;

bool Verifier::Initialize() {

    // Prepare Input Reader //
    // PENDING: could refactor into a function

    fInput_Model = ROOT::RNTupleModel::Create();

    fInput_Event.AddFieldsTo(fInput_Model.get(), fSettings.IsMC);

    if (fSettings.IsMC) fInput_McParticle.AddFieldsTo(fInput_Model.get());

    fInput_OnTheFlyLambda.AddFieldsTo(fInput_Model.get(), fSettings.IsMC);

    fReader =
        ROOT::RNTupleReader::Open(std::move(fInput_Model), E2R::Name_OutputRNT, fSettings.PathInputFiles.front());  // PENDING: handle multiple files?

    // Prepare Output //

    if (!PrepareOutputFile()) return false;

    PrepareOutputHistograms();

    // Prepare Output Writers //
    // PENDING: could refactor into a function

    // -- reconstructed
    fOutput_Model = ROOT::RNTupleModel::Create();
    fOutput_LambdaPair.AddFieldsTo(fOutput_Model.get(), fSettings.IsMC);
    fWriter = ROOT::RNTupleWriter::Append(std::move(fOutput_Model), R2DS::Name_LambdaPairRNT, *fOutput_File);
    // -- mc
    if (fSettings.IsMC) {
        fOutput_Model_InjectedHdib = ROOT::RNTupleModel::Create();
        fOutput_InjectedHdib.AddFieldsTo(fOutput_Model_InjectedHdib.get(), false);
        fWriter_InjectedHdib = ROOT::RNTupleWriter::Append(std::move(fOutput_Model_InjectedHdib), R2DS::Name_InjectedHdibRNT, *fOutput_File);
    }

    Logger::Info(__FUNCTION__, "Verifier initialized successfully.");

    return true;
}

// ## OUTPUT ZONE ## //

bool Verifier::PrepareOutputFile() {

    const std::filesystem::path output_path(fSettings.PathOutputFile);
    if (output_path.has_parent_path()) std::filesystem::create_directories(output_path.parent_path());

    fOutput_File = std::unique_ptr<TFile>(TFile::Open(fSettings.PathOutputFile.c_str(), "RECREATE"));
    if (!fOutput_File || fOutput_File->IsZombie()) {
        Logger::Error(__FUNCTION__, "Couldn't create TFile {}", fSettings.PathOutputFile);
        return false;
    }

    return true;
}

void Verifier::PrepareOutputHistograms() {

    // event counter //
    fHist_EventCounter = std::make_unique<TH1D>("N_Events", ";;N_Events", 1, 0, 1);

    constexpr int x_nbins = 20;
    constexpr float x_min = 0.;
    constexpr float x_max = 20.;
    constexpr const char* hist_title = ";Cut N;N Passed Cut";

    fHist_CutFlow = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", DB::Particles::Particle("AntiProton").acronym).c_str(), hist_title, x_nbins, x_min, x_max);
    fHist_CutFlow_AntiChannel = std::make_unique<TH1D>(  //
        std::format("CutFlow_{}", DB::Particles::Particle("AntiLambda").acronym).c_str(), hist_title, x_nbins, x_min, x_max);
}

// ## Helpers ## //

void Verifier::Assign_Event() {
    *fOutput_LambdaPair.RunNumber = *fInput_Event.RunNumber;
    *fOutput_LambdaPair.DirNumber = *fInput_Event.DirNumber;
    if (!fSettings.IsMC) *fOutput_LambdaPair.DirNumberB = *fInput_Event.DirNumberB;
    *fOutput_LambdaPair.EventNumber = *fInput_Event.EventNumber;
    *fOutput_LambdaPair.Centrality = *fInput_Event.Centrality;
    *fOutput_LambdaPair.MagneticField = *fInput_Event.MagneticField;
    *fOutput_LambdaPair.PV_X = *fInput_Event.PV_X;
    *fOutput_LambdaPair.PV_Y = *fInput_Event.PV_Y;
    *fOutput_LambdaPair.PV_Z = *fInput_Event.PV_Z;
    if (fSettings.IsMC) {
        *fOutput_LambdaPair.MC_PV_X = *fInput_Event.MC_PV_X;
        *fOutput_LambdaPair.MC_PV_Y = *fInput_Event.MC_PV_Y;
        *fOutput_LambdaPair.MC_PV_Z = *fInput_Event.MC_PV_Z;
    }
}

void Verifier::Assign_Candidate(const KF::LambdaPair& hdibaryon, bool anti_channel) {
    *fOutput_LambdaPair.Decay_X = static_cast<float>(hdibaryon.X());
    *fOutput_LambdaPair.Decay_Y = static_cast<float>(hdibaryon.Y());
    *fOutput_LambdaPair.Decay_Z = static_cast<float>(hdibaryon.Z());
    *fOutput_LambdaPair.Px = static_cast<float>(hdibaryon.Px());
    *fOutput_LambdaPair.Py = static_cast<float>(hdibaryon.Py());
    *fOutput_LambdaPair.Pz = static_cast<float>(hdibaryon.Pz());
    *fOutput_LambdaPair.Energy = static_cast<float>(hdibaryon.E());
    *fOutput_LambdaPair.AntiChannel = anti_channel;
}

void Verifier::Assign(const Vector::OnTheFlyLambdaView& lambda, const Vector::McParticleView* mc_lambda, const Vector::McParticleView* mc_neg,
                      const Vector::McParticleView* mc_pos, Flat::OnTheFlyLambda& out) {
    // -- (anti)lambda info
    *out.OnTheFlyEntry = lambda.OnTheFlyEntry();
    *out.Decay_X = lambda.Decay_X();
    *out.Decay_Y = lambda.Decay_Y();
    *out.Decay_Z = lambda.Decay_Z();
    *out.Px = lambda.Px();
    *out.Py = lambda.Py();
    *out.Pz = lambda.Pz();
    *out.Energy = static_cast<float>(lambda.Energy());
    *out.DcaV0Daughters = lambda.DcaV0Daughters();
    if (fSettings.IsMC) {
        if (mc_lambda != nullptr) {
            *out.McEntry = static_cast<int>(mc_lambda->Entry);
            *out.PdgCode = mc_lambda->PdgCode();
            *out.MC_Px = mc_lambda->Px();
            *out.MC_Py = mc_lambda->Py();
            *out.MC_Pz = mc_lambda->Pz();
            *out.MC_Energy = mc_lambda->Energy();
            *out.Origin_X = mc_lambda->Origin_X();
            *out.Origin_Y = mc_lambda->Origin_Y();
            *out.Origin_Z = mc_lambda->Origin_Z();
            *out.MC_Decay_X = static_cast<float>(mc_lambda->Decay_X());
            *out.MC_Decay_Y = static_cast<float>(mc_lambda->Decay_Y());
            *out.MC_Decay_Z = static_cast<float>(mc_lambda->Decay_Z());
            *out.IsTrue = mc_lambda->IsTrue();
            *out.IsSignal = mc_lambda->IsSignal();
            *out.IsHybrid = false;  // PENDING!
            *out.InjectionID = mc_lambda->InjectionID();
        } else {
            *out.McEntry = Common::DummyInt;
            *out.PdgCode = Common::DummyInt;
            *out.MC_Px = Common::DummyFloat;
            *out.MC_Py = Common::DummyFloat;
            *out.MC_Pz = Common::DummyFloat;
            *out.MC_Energy = Common::DummyFloat;
            *out.Origin_X = Common::DummyFloat;
            *out.Origin_Y = Common::DummyFloat;
            *out.Origin_Z = Common::DummyFloat;
            *out.MC_Decay_X = Common::DummyFloat;
            *out.MC_Decay_Y = Common::DummyFloat;
            *out.MC_Decay_Z = Common::DummyFloat;
            *out.IsTrue = false;
            *out.IsSignal = false;
            *out.IsHybrid = false;
            *out.InjectionID = false;
        }
    }
    // -- negative daughter
    *out.Neg_EsdEntry = lambda.Neg_EsdEntry();
    *out.Neg_PCAwrtV0_Px = lambda.Neg_PCAwrtV0_Px();
    *out.Neg_PCAwrtV0_Py = lambda.Neg_PCAwrtV0_Py();
    *out.Neg_PCAwrtV0_Pz = lambda.Neg_PCAwrtV0_Pz();
    *out.Neg_PreDCAxy = lambda.Neg_PreDCAxy();
    *out.Neg_PreDCAz = lambda.Neg_PreDCAz();
    *out.Neg_NSigmaProton = lambda.Neg_NSigmaProton();
    *out.Neg_NSigmaKaon = lambda.Neg_NSigmaKaon();
    *out.Neg_NSigmaPion = lambda.Neg_NSigmaPion();
    if (fSettings.IsMC) {
        // NOTE: `mc_neg` cannot be `nullptr`, by construction
        *out.Neg_McEntry = lambda.Neg_McEntry();
        *out.Neg_PdgCode = mc_neg->PdgCode();
        *out.Neg_Origin_X = mc_neg->Origin_X();
        *out.Neg_Origin_Y = mc_neg->Origin_Y();
        *out.Neg_Origin_Z = mc_neg->Origin_Z();
        *out.Neg_MC_Px = mc_neg->Px();
        *out.Neg_MC_Py = mc_neg->Py();
        *out.Neg_MC_Pz = mc_neg->Pz();
        *out.Neg_MC_Energy = mc_neg->Energy();
        *out.Neg_IsTrue = mc_neg->IsTrue();
        *out.Neg_IsSignal = mc_neg->IsSignal();
        *out.Neg_InjectionID = mc_neg->InjectionID();
        *out.Neg_Mother_McEntry = mc_neg->Mother_McEntry();
        *out.Neg_Mother_PdgCode = mc_neg->Mother_PdgCode();
    }
    // -- positive daughter
    *out.Pos_EsdEntry = lambda.Pos_EsdEntry();
    *out.Pos_PCAwrtV0_Px = lambda.Pos_PCAwrtV0_Px();
    *out.Pos_PCAwrtV0_Py = lambda.Pos_PCAwrtV0_Py();
    *out.Pos_PCAwrtV0_Pz = lambda.Pos_PCAwrtV0_Pz();
    *out.Pos_PreDCAxy = lambda.Pos_PreDCAxy();
    *out.Pos_PreDCAz = lambda.Pos_PreDCAz();
    *out.Pos_NSigmaProton = lambda.Pos_NSigmaProton();
    *out.Pos_NSigmaKaon = lambda.Pos_NSigmaKaon();
    *out.Pos_NSigmaPion = lambda.Pos_NSigmaPion();
    if (fSettings.IsMC) {
        // NOTE: `mc_pos` cannot be `nullptr`, by construction
        *out.Pos_McEntry = lambda.Pos_McEntry();
        *out.Pos_PdgCode = mc_pos->PdgCode();
        *out.Pos_Origin_X = mc_pos->Origin_X();
        *out.Pos_Origin_Y = mc_pos->Origin_Y();
        *out.Pos_Origin_Z = mc_pos->Origin_Z();
        *out.Pos_MC_Px = mc_pos->Px();
        *out.Pos_MC_Py = mc_pos->Py();
        *out.Pos_MC_Pz = mc_pos->Pz();
        *out.Pos_MC_Energy = mc_pos->Energy();
        *out.Pos_IsTrue = mc_pos->IsTrue();
        *out.Pos_IsSignal = mc_pos->IsSignal();
        *out.Pos_InjectionID = mc_pos->InjectionID();
        *out.Pos_Mother_McEntry = mc_pos->Mother_McEntry();
        *out.Pos_Mother_PdgCode = mc_pos->Mother_PdgCode();
    }
}

void Verifier::Assign_InjectedHdib(const Vector::McParticleView* mc, Flat::InjectedHdib& out, bool embedded_to_rec) {
    if (!embedded_to_rec) {
        *out.RunNumber = *fInput_Event.RunNumber;
        *out.DirNumber = *fInput_Event.DirNumber;
        *out.EventNumber = *fInput_Event.EventNumber;
    }
    if (mc != nullptr) {
        *out.InjectionID = static_cast<int>(mc->InjectionID());
        *out.PdgCode = static_cast<int>(mc->PdgCode());
        *out.Decay_X = static_cast<float>(mc->Decay_X());
        *out.Decay_Y = static_cast<float>(mc->Decay_Y());
        *out.Decay_Z = static_cast<float>(mc->Decay_Z());
        *out.Px = mc->Px();
        *out.Py = mc->Py();
        *out.Pz = mc->Pz();
        *out.Energy = static_cast<float>(mc->Energy());
    } else {
        *out.InjectionID = Common::DummyInt;
        *out.PdgCode = Common::DummyInt;
        *out.Decay_X = Common::DummyFloat;
        *out.Decay_Y = Common::DummyFloat;
        *out.Decay_Z = Common::DummyFloat;
        *out.Px = Common::DummyFloat;
        *out.Py = Common::DummyFloat;
        *out.Pz = Common::DummyFloat;
        *out.Energy = Common::DummyFloat;
    }
}

// ## Event ZONE ## //

void Verifier::ProcessEvent() {  //
    fHist_EventCounter->Fill(0.);
    fPrimaryVertex.SetCoordinates(*fInput_Event.PV_X, *fInput_Event.PV_Y, *fInput_Event.PV_Z);
}

// ## Injected ZONE ## //

void Verifier::ProcessInjected() {

    // prepare view //
    Vector::McParticleView mc_view(&fInput_McParticle);

    // loop over mc particles //
    for (std::size_t entry_mc = 0; entry_mc < mc_view.Size(); ++entry_mc) {
        mc_view.Entry = entry_mc;  // NOTE: prevent `CacheCalculations(...)` machinery
        // select only injected h-dibaryons //
        if (std::abs(mc_view.PdgCode()) != DB::Particles::Particle("Hdibaryon").pdg_code) continue;
        // cache `decay_{x/y/z}` //
        mc_view.CacheDescendantsInfo();
        // assign values //
        Assign_InjectedHdib(&mc_view, fOutput_InjectedHdib, false);
        // fill rntuple //
        fWriter_InjectedHdib->Fill();
    }
}

// ## H-dibaryon ZONE ## //

void Verifier::VerifyLambdaPair(bool anti_channel) {

    // determine rules based on reconstruction of lambdas or anti-lambdas //
    DB::Particles::Definition pid_lambda = anti_channel ? DB::Particles::Particle("AntiLambda") : DB::Particles::Particle("Lambda");
    DB::Particles::Definition pid_neg = anti_channel ? DB::Particles::Particle("AntiProton") : DB::Particles::Particle("PiMinus");
    DB::Particles::Definition pid_pos = anti_channel ? DB::Particles::Particle("PiPlus") : DB::Particles::Particle("Proton");
    TH1D* hist_cut_flow = anti_channel ? fHist_CutFlow_AntiChannel.get() : fHist_CutFlow.get();

    // prepare views //
    Vector::OnTheFlyLambdaView lambda1(&fInput_OnTheFlyLambda);
    Vector::OnTheFlyLambdaView lambda2(&fInput_OnTheFlyLambda);
    const std::size_t n_lambdas = lambda1.Size();
    // -- mc views
    std::unique_ptr<Vector::McParticleView> mc_lambda1_neg = nullptr;
    std::unique_ptr<Vector::McParticleView> mc_lambda1_pos = nullptr;
    std::unique_ptr<Vector::McParticleView> mc_lambda1 = nullptr;
    std::unique_ptr<Vector::McParticleView> mc_lambda2_neg = nullptr;
    std::unique_ptr<Vector::McParticleView> mc_lambda2_pos = nullptr;
    std::unique_ptr<Vector::McParticleView> mc_lambda2 = nullptr;
    if (fSettings.IsMC) {
        mc_lambda1_neg = std::make_unique<Vector::McParticleView>(&fInput_McParticle);
        mc_lambda1_pos = std::make_unique<Vector::McParticleView>(&fInput_McParticle);
        mc_lambda2_neg = std::make_unique<Vector::McParticleView>(&fInput_McParticle);
        mc_lambda2_pos = std::make_unique<Vector::McParticleView>(&fInput_McParticle);
    }

    // loop over all possible pairs of (anti)lambdas //
    for (std::size_t entry_lambda1 = 0; entry_lambda1 + 1 < n_lambdas; ++entry_lambda1) {
        // cache //
        lambda1.CacheCalculations(entry_lambda1, fPrimaryVertex, pid_neg, pid_pos);

        if (fSettings.IsMC) {
            // -- lambda1's negative daughter
            auto mc_entry_lambda1_neg = lambda1.Neg_McEntry();  // NOTE: cannot be invalid, by construction
            mc_lambda1_neg->Entry = mc_entry_lambda1_neg;       // NOTE: prevent full `CacheCalculations()` machinery
            mc_lambda1_neg->CacheAscendantsInfo();
            mc_lambda1_neg->CacheDescendantsInfo();
            mc_lambda1_neg->Classify_AsInHdibaryonSimulations(anti_channel, pid_neg.pdg_code);
            // -- lambda1's positive daughter
            auto mc_entry_lambda1_pos = lambda1.Pos_McEntry();  // NOTE: cannot be invalid, by construction
            mc_lambda1_pos->Entry = mc_entry_lambda1_pos;       // NOTE: prevent full `CacheCalculations()` machinery
            mc_lambda1_pos->CacheAscendantsInfo();
            mc_lambda1_pos->CacheDescendantsInfo();
            mc_lambda1_pos->Classify_AsInHdibaryonSimulations(anti_channel, pid_pos.pdg_code);
            // -- lambda1
            mc_lambda1 =
                mc_lambda1_neg->Mother_McEntry() == mc_lambda1_pos->Mother_McEntry()
                    ? std::make_unique<Vector::McParticleView>(&fInput_McParticle, static_cast<std::size_t>(mc_lambda1_neg->Mother_McEntry()))
                    : nullptr;
            if (mc_lambda1 != nullptr) {
                mc_lambda1->CacheAscendantsInfo();
                mc_lambda1->CacheDescendantsInfo();
                mc_lambda1->Classify_AsInHdibaryonSimulations(anti_channel, pid_lambda.pdg_code);
            }
        }

        for (std::size_t entry_lambda2 = entry_lambda1; entry_lambda2 < n_lambdas; ++entry_lambda2) {
            // NOTE: sanity check is not needed, because loops don't intersect
            // cache //
            lambda2.CacheCalculations(entry_lambda2, fPrimaryVertex, pid_neg, pid_pos);

            if (fSettings.IsMC) {
                // -- lambda2's negative daughter
                auto mc_entry_lambda2_neg = lambda2.Neg_McEntry();  // NOTE: cannot be invalid, by construction
                mc_lambda2_neg->Entry = mc_entry_lambda2_neg;       // NOTE: prevent full `CacheCalculations()` machinery
                mc_lambda2_neg->CacheAscendantsInfo();
                mc_lambda2_neg->CacheDescendantsInfo();
                mc_lambda2_neg->Classify_AsInHdibaryonSimulations(anti_channel, pid_neg.pdg_code);
                // -- lambda2's positive daughter
                auto mc_entry_lambda2_pos = lambda2.Pos_McEntry();  // NOTE: cannot be invalid, by construction
                mc_lambda2_pos->Entry = mc_entry_lambda2_pos;       // NOTE: prevent full `CacheCalculations()` machinery
                mc_lambda2_pos->CacheAscendantsInfo();
                mc_lambda2_pos->CacheDescendantsInfo();
                mc_lambda2_pos->Classify_AsInHdibaryonSimulations(anti_channel, pid_pos.pdg_code);
                // -- lambda2
                mc_lambda2 =
                    mc_lambda2_neg->Mother_McEntry() == mc_lambda2_pos->Mother_McEntry()
                        ? std::make_unique<Vector::McParticleView>(&fInput_McParticle, static_cast<std::size_t>(mc_lambda2_neg->Mother_McEntry()))
                        : nullptr;
                if (mc_lambda2 != nullptr) {
                    mc_lambda2->CacheAscendantsInfo();
                    mc_lambda2->CacheDescendantsInfo();
                    mc_lambda2->Classify_AsInHdibaryonSimulations(anti_channel, pid_lambda.pdg_code);
                }
            }

            // PCAs //
            auto [seed_lambda1, seed_lambda2] = Seeder::LineLine::FullPCAs(lambda1, lambda2);

            // create composite particle //
            KF::LambdaPair hdibaryon(lambda1, lambda2, seed_lambda1.pca, seed_lambda2.pca);

            // apply cuts //
            if (!Cuts(hdibaryon, hist_cut_flow)) continue;

            std::unique_ptr<MC::LambdaPair> mc_hdibaryon = nullptr;
            if (fSettings.IsMC) {
                mc_hdibaryon = std::make_unique<MC::LambdaPair>(mc_lambda1.get(), mc_lambda1_neg.get(), mc_lambda1_pos.get(), mc_lambda2.get(),
                                                                mc_lambda2_neg.get(), mc_lambda2_pos.get());
            }

            // store //
            Assign(hdibaryon, mc_hdibaryon.get(), anti_channel);

            // fill //
            fWriter->Fill();
        }  // end of loop over pos
    }  // end of loop over neg
}

bool Verifier::Cuts(const KF::LambdaPair& hdibaryon, TH1D* cut_flow_hist) const {
    cut_flow_hist->Fill(0.);

    if (Common::Math::SquaredDistance(hdibaryon.Lambda1_PCAwrtDV.GetXYZ_AsROOT(), hdibaryon.Lambda2_PCAwrtDV.GetXYZ_AsROOT()) >
        Cuts::LambdaPair::Max_DCAbtwDau * Cuts::LambdaPair::Max_DCAbtwDau)
        return false;
    cut_flow_hist->Fill(1.);

    auto mass = hdibaryon.M();  // cached
    if (mass < Cuts::LambdaPair::Min_Mass || mass > Cuts::LambdaPair::Max_Mass) return false;
    cut_flow_hist->Fill(2.);

    return true;
}

std::unique_ptr<Vector::McParticleView> Verifier::LinkInjectedSignal(const MC::LambdaPair& mc_hdibaryon) {
    if (!mc_hdibaryon.lambda1->IsSignal()) return nullptr;
    if (!mc_hdibaryon.lambda2->IsSignal()) return nullptr;
    if (mc_hdibaryon.lambda1->Mother_McEntry() != mc_hdibaryon.lambda2->Mother_McEntry()) return nullptr;
    if (mc_hdibaryon.lambda1->Mother_McEntry() > Common::DummyInt) return nullptr;
    return std::make_unique<Vector::McParticleView>(&fInput_McParticle, static_cast<std::size_t>(mc_hdibaryon.lambda1->Mother_McEntry()));
}

void Verifier::Assign(const KF::LambdaPair& hdibaryon, const MC::LambdaPair* mc_hdibaryon, bool anti_channel) {
    // -- event info
    Assign_Event();
    // -- candidate info
    Assign_Candidate(hdibaryon, anti_channel);
    // -- lambda1
    Assign(hdibaryon.Lambda1, mc_hdibaryon->lambda1, mc_hdibaryon->lambda1_neg, mc_hdibaryon->lambda1_pos, fOutput_LambdaPair.Lambda1);
    *fOutput_LambdaPair.Lambda1_PCAwrtDV_X = static_cast<float>(hdibaryon.Lambda1_PCAwrtDV.X());
    *fOutput_LambdaPair.Lambda1_PCAwrtDV_Y = static_cast<float>(hdibaryon.Lambda1_PCAwrtDV.Y());
    *fOutput_LambdaPair.Lambda1_PCAwrtDV_Z = static_cast<float>(hdibaryon.Lambda1_PCAwrtDV.Z());
    *fOutput_LambdaPair.Lambda1_PCAwrtDV_Px = static_cast<float>(hdibaryon.Lambda1_PCAwrtDV.Px());
    *fOutput_LambdaPair.Lambda1_PCAwrtDV_Py = static_cast<float>(hdibaryon.Lambda1_PCAwrtDV.Py());
    *fOutput_LambdaPair.Lambda1_PCAwrtDV_Pz = static_cast<float>(hdibaryon.Lambda1_PCAwrtDV.Pz());
    // -- lambda2
    Assign(hdibaryon.Lambda2, mc_hdibaryon->lambda2, mc_hdibaryon->lambda2_neg, mc_hdibaryon->lambda2_pos, fOutput_LambdaPair.Lambda2);
    *fOutput_LambdaPair.Lambda2_PCAwrtDV_X = static_cast<float>(hdibaryon.Lambda2_PCAwrtDV.X());
    *fOutput_LambdaPair.Lambda2_PCAwrtDV_Y = static_cast<float>(hdibaryon.Lambda2_PCAwrtDV.Y());
    *fOutput_LambdaPair.Lambda2_PCAwrtDV_Z = static_cast<float>(hdibaryon.Lambda2_PCAwrtDV.Z());
    *fOutput_LambdaPair.Lambda2_PCAwrtDV_Px = static_cast<float>(hdibaryon.Lambda2_PCAwrtDV.Px());
    *fOutput_LambdaPair.Lambda2_PCAwrtDV_Py = static_cast<float>(hdibaryon.Lambda2_PCAwrtDV.Py());
    *fOutput_LambdaPair.Lambda2_PCAwrtDV_Pz = static_cast<float>(hdibaryon.Lambda2_PCAwrtDV.Pz());
    if (!fSettings.IsMC) return;
    Assign_InjectedHdib(LinkInjectedSignal(*mc_hdibaryon).get(), fOutput_LambdaPair.MC, true);
}

// ## END OF CYCLES ## //

void Verifier::EndOfAnalysis() {

    Logger::Info(__FUNCTION__, "The following objects have been written into TFile {}:", fSettings.PathOutputFile);

    if (fSettings.IsMC) {
        Logger::Info(__FUNCTION__, "- RNTuple \"{}\"", R2DS::Name_InjectedHdibRNT);
    }
    Logger::Info(__FUNCTION__, "- RNTuple \"{}\"", R2DS::Name_LambdaPairRNT);

    // write histograms //

    // -- event counter
    fHist_EventCounter->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_EventCounter->GetName());

    // -- cut flows
    fHist_CutFlow->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow->GetName());
    fHist_CutFlow_AntiChannel->Write();
    Logger::Info(__FUNCTION__, "- TH1D  \"{}\"", fHist_CutFlow_AntiChannel->GetName());

    Logger::Info(__FUNCTION__, "All done.");
}

}  // namespace R2DS
