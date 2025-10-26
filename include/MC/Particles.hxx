#ifndef T2S_MC_PARTICLES_HXX
#define T2S_MC_PARTICLES_HXX

#include <cmath>

#include <KFParticle.hxx>
#include <KFParticle_Math.hxx>

#include "DataFormats/Events.hxx"
#include "DataFormats/Injected.hxx"
#include "DataFormats/PackedEvents.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::MC {

struct alignas(T2S_SIMD_ALIGN) Particle {
    // constructors //
    Particle() = default;
    Particle(int entry, int pdg_code, int mother_entry, int mother_pdg_code, float x, float y, float z, float px, float py, float pz, float energy)
        : Entry{entry},
          Mother_Entry{mother_entry},
          X{x},
          Y{y},
          Z{z},
          Px{px},
          Py{py},
          Pz{pz},
          Energy{energy},
          PdgCode{pdg_code},
          Mother_PdgCode{mother_pdg_code} {}

    // utilities //
    void FillBasic(const DF::SOV::MCParticles& df, int mc_idx) {
        if (mc_idx >= 0) {
            Entry = mc_idx;
            X = df.X->at(mc_idx);
            Y = df.Y->at(mc_idx);
            Z = df.Z->at(mc_idx);
            Px = df.Px->at(mc_idx);
            Py = df.Py->at(mc_idx);
            Pz = df.Pz->at(mc_idx);
            Energy = df.Energy->at(mc_idx);

            PdgCode = df.PdgCode->at(mc_idx);
            Mother_Entry = df.MotherEntry->at(mc_idx);
            if (Mother_Entry >= 0) Mother_PdgCode = df.PdgCode->at(Mother_Entry);
        }
    }

    // member vars //
    int Entry{Const::DummyInt};
    int Mother_Entry{Const::DummyInt};
    float X{Const::DummyFloat};  // origin
    float Y{Const::DummyFloat};  // origin
    float Z{Const::DummyFloat};  // origin
    float Px{Const::DummyFloat};
    float Py{Const::DummyFloat};
    float Pz{Const::DummyFloat};
    float Energy{Const::DummyFloat};
    int PdgCode{Const::DummyInt};
    int Mother_PdgCode{Const::DummyInt};
};

// `MC::Track`
struct alignas(T2S_SIMD_ALIGN) Track : MC::Particle {
    // constructors //
    Track() = default;
    Track(const DF::SOV::MCParticles& df, int mc_idx, EParticle pid_hypothesis) { Init(df, mc_idx, pid_hypothesis); }
    Track(const DF::Packed::LinkedTracks& df, int idx)
        : MC::Particle{df.Entry->at(idx),
                       df.PdgCode->at(idx),
                       df.Mother_Entry->at(idx),
                       df.Mother_PdgCode->at(idx),
                       df.X->at(idx),
                       df.Y->at(idx),
                       df.Z->at(idx),
                       df.Px->at(idx),
                       df.Py->at(idx),
                       df.Pz->at(idx),
                       df.Energy->at(idx)},
          GrandMother_Entry{df.GrandMother_Entry->at(idx)},
          GrandMother_PdgCode{df.GrandMother_PdgCode->at(idx)},
          ReactionID(df.ReactionID->at(idx)),
          IsTrue(df.IsTrue->at(idx)),
          IsSignal(df.IsSignal->at(idx)),
          IsSecondary(df.IsSecondary->at(idx)) {}
    Track(int entry, int pdg_code, float px, float py, float pz, bool is_true, bool is_signal, bool is_secondary, int reaction_id)
        : MC::Particle{entry, pdg_code, Const::DummyInt,  Const::DummyInt, Const::DummyFloat, Const::DummyFloat, Const::DummyFloat, px,
                       py,    pz,       Const::DummyFloat},
          ReactionID(reaction_id),
          IsTrue(is_true),
          IsSignal(is_signal),
          IsSecondary(is_secondary) {}

    // utilities //
    void Init(const DF::SOV::MCParticles& df, int mc_idx, EParticle pid_hypothesis) {
        FillBasic(df, mc_idx);
        FillDerived(df, mc_idx, pid_hypothesis);
    }
    void FillDerived(const DF::SOV::MCParticles& df, int mc_idx, EParticle pid_hypothesis) {
        if (mc_idx >= 0) {
            IsTrue = PdgCode == Tree2Secondaries::Particle::PdgCode[pid_hypothesis];
            IsSignal = IsTrue && df.Generator->at(mc_idx) == 2;
            IsSecondary = df.IsSecFromMat->at(mc_idx) || df.IsSecFromWeak->at(mc_idx) || IsSignal;
            if (Mother_Entry >= 0) {
                if (IsSignal) ReactionID = df.Status->at(Mother_Entry);
                GrandMother_Entry = df.MotherEntry->at(Mother_Entry);
                if (GrandMother_Entry >= 0) GrandMother_PdgCode = df.PdgCode->at(GrandMother_Entry);
            }
        }
    }

    // member vars //
    int GrandMother_Entry{Const::DummyInt};
    int GrandMother_PdgCode{Const::DummyInt};
    int ReactionID{Const::DummyInt};
    bool IsTrue{false};
    bool IsSignal{false};
    bool IsSecondary{false};
};

// `MC::V0`
struct alignas(T2S_SIMD_ALIGN) V0 : MC::Particle {
    // constructors //
    V0() = default;
    V0(const DF::SOV::MCParticles& df, int mc_neg, int mc_pos, EParticle v0_hypothesis, EParticle neg_hypothesis, EParticle pos_hypothesis) {
        if (mc_neg < 0 || mc_pos < 0) return;  // protection
        neg.Init(df, mc_neg, neg_hypothesis);
        pos.Init(df, mc_pos, pos_hypothesis);

        int mc_mother_neg{df.MotherEntry->at(mc_neg)};
        int mc_mother_pos{df.MotherEntry->at(mc_pos)};
        if (mc_mother_neg >= 0 && mc_mother_neg == mc_mother_pos) {
            int mc_v0{mc_mother_neg};
            FillBasic(df, mc_v0);

            // fill derived props //
            IsTrue = PdgCode == Tree2Secondaries::Particle::PdgCode[v0_hypothesis] &&
                     neg.PdgCode == Tree2Secondaries::Particle::PdgCode[neg_hypothesis] &&
                     pos.PdgCode == Tree2Secondaries::Particle::PdgCode[pos_hypothesis];
            IsSignal = IsTrue && df.Generator->at(mc_v0) == 2 && neg.ReactionID == pos.ReactionID;
            IsSecondary = df.IsSecFromMat->at(mc_v0) || df.IsSecFromWeak->at(mc_v0) || IsSignal;
            if (IsSignal)
                ReactionID = df.Status->at(mc_v0);
            else
                IsHybrid = (neg.IsSignal && !pos.IsSignal) || (!neg.IsSignal && pos.IsSignal) ||
                           (neg.IsSignal && pos.IsSignal && neg.Mother_Entry != pos.Mother_Entry);
        }
    }
    V0(const DF::Packed::LinkedV0s& df, int idx)
        : MC::Particle{df.Entry->at(idx),
                       df.PdgCode->at(idx),
                       df.Mother_Entry->at(idx),
                       df.Mother_PdgCode->at(idx),
                       df.X->at(idx),
                       df.Y->at(idx),
                       df.Z->at(idx),
                       df.Px->at(idx),
                       df.Py->at(idx),
                       df.Pz->at(idx),
                       df.Energy->at(idx)},
          neg{df.Neg.Entry->at(idx),
              df.Neg.PdgCode->at(idx),
              df.Neg.Px->at(idx),
              df.Neg.Py->at(idx),
              df.Neg.Pz->at(idx),
              static_cast<bool>(df.Neg.IsTrue->at(idx)),
              static_cast<bool>(df.Neg.IsSignal->at(idx)),
              static_cast<bool>(df.Neg.IsSecondary->at(idx)),
              df.Neg.ReactionID->at(idx)},
          pos{df.Pos.Entry->at(idx),
              df.Pos.PdgCode->at(idx),
              df.Pos.Px->at(idx),
              df.Pos.Py->at(idx),
              df.Pos.Pz->at(idx),
              static_cast<bool>(df.Pos.IsTrue->at(idx)),
              static_cast<bool>(df.Pos.IsSignal->at(idx)),
              static_cast<bool>(df.Pos.IsSecondary->at(idx)),
              df.Pos.ReactionID->at(idx)},
          ReactionID{df.ReactionID->at(idx)},
          IsTrue{static_cast<bool>(df.IsTrue->at(idx))},
          IsSignal{static_cast<bool>(df.IsSignal->at(idx))},
          IsSecondary{static_cast<bool>(df.IsSecondary->at(idx))},
          IsHybrid{static_cast<bool>(df.IsHybrid->at(idx))} {}

    // utilities //
    [[nodiscard]] float DecayX() const { return neg.Entry >= 0 ? neg.X : pos.X; };
    [[nodiscard]] float DecayY() const { return neg.Entry >= 0 ? neg.Y : pos.Y; };
    [[nodiscard]] float DecayZ() const { return neg.Entry >= 0 ? neg.Z : pos.Z; };

    // member vars //
    MC::Track neg;
    MC::Track pos;
    int ReactionID{Const::DummyInt};
    bool IsTrue{false};
    bool IsSignal{false};
    bool IsSecondary{false};
    bool IsHybrid{false};
};

// `MC::Sexaquark`
struct alignas(T2S_SIMD_ALIGN) Sexaquark {
    // constructors //
    Sexaquark() = default;
    Sexaquark(const DF::SOV::Injected& df, double mass_sexaquark, double mass_nucleon, int reaction_id) {
        FillInfo_BeforeReaction(df, mass_sexaquark, mass_nucleon, reaction_id);
    }

    // utilities //
    void FillInfo_BeforeReaction(const DF::SOV::Injected& df, double mass_sexaquark, double mass_nucleon, int reaction_id) {
        if (reaction_id >= 600 && reaction_id < 620) {
            int reaction_idx{reaction_id - 600};
            BeforePx = df.Px->at(reaction_idx);
            BeforePy = df.Py->at(reaction_idx);
            BeforePz = df.Pz->at(reaction_idx);
            BeforeE = std::sqrt(mass_sexaquark * mass_sexaquark + BeforePx * BeforePx + BeforePy * BeforePy + BeforePz * BeforePz);
            ReactionID = df.ReactionID->at(reaction_idx);
            NucleonPx = df.Nucleon_Px->at(reaction_idx);
            NucleonPy = df.Nucleon_Py->at(reaction_idx);
            NucleonPz = df.Nucleon_Pz->at(reaction_idx);
            NucleonE = std::sqrt(mass_nucleon * mass_nucleon + NucleonPx * NucleonPx + NucleonPy * NucleonPy + NucleonPz * NucleonPz);
        }
    }

    // member vars //
    double BeforePx{Const::DummyDouble};
    double BeforePy{Const::DummyDouble};
    double BeforePz{Const::DummyDouble};
    double BeforeE{Const::DummyDouble};
    double NucleonPx{Const::DummyDouble};
    double NucleonPy{Const::DummyDouble};
    double NucleonPz{Const::DummyDouble};
    double NucleonE{Const::DummyDouble};
    int ReactionID{Const::DummyInt};
};

// `MC::ChannelA`
// Reminder: AntiSexaquark + Neutron -> AntiLambda + KaonZeroShort
struct alignas(T2S_SIMD_ALIGN) ChannelA : MC::Sexaquark {
    // constructors //
    ChannelA() = default;
    ChannelA(const DF::SOV::Injected& df, double mass_sexaquark, const V0& v0a, const V0& v0b)
        : V0A{v0a},
          V0B{v0b},
          X{v0a.X},
          Y{v0a.Y},
          Z{v0a.Z},
          AfterPx{v0a.Px + v0b.Px},
          AfterPy{v0a.Py + v0b.Py},
          AfterPz{v0a.Pz + v0b.Pz},
          AfterE{v0a.Energy + v0b.Energy},
          IsSignal{v0a.IsSignal && v0b.IsSignal && v0a.ReactionID == v0b.ReactionID} {
        if (IsSignal) {
            FillInfo_BeforeReaction(df, mass_sexaquark, Tree2Secondaries::Particle::Mass[EParticle::Neutron], v0a.ReactionID);
        } else {
            IsHybrid = (v0a.IsSignal && !v0b.IsSignal) || (!v0a.IsSignal && v0b.IsSignal) ||
                       (v0a.IsSignal && v0b.IsSignal && v0a.ReactionID != v0b.ReactionID) || v0a.IsHybrid || v0b.IsHybrid;
        }
    }

    // member vars //
    MC::V0 V0A;
    MC::V0 V0B;
    double X{Const::DummyDouble};
    double Y{Const::DummyDouble};
    double Z{Const::DummyDouble};
    double AfterPx{Const::DummyDouble};
    double AfterPy{Const::DummyDouble};
    double AfterPz{Const::DummyDouble};
    double AfterE{Const::DummyDouble};
    bool IsSignal{false};
    bool IsHybrid{false};
};

// `MC::ChannelD`
// Reminder: AntiSexaquark + Proton -> AntiLambda + PosKaon
struct alignas(T2S_SIMD_ALIGN) ChannelD : MC::Sexaquark {
    // constructors //
    ChannelD() = default;
    ChannelD(const DF::SOV::Injected& df, double mass_sexaquark, const V0& v0, const Track& kaon)
        : V0{v0},
          Kaon{kaon},
          X{v0.X},
          Y{v0.Y},
          Z{v0.Z},
          AfterPx{v0.Px + kaon.Px},
          AfterPy{v0.Py + kaon.Py},
          AfterPz{v0.Pz + kaon.Pz},
          AfterE{v0.Energy + kaon.Energy},
          IsSignal{v0.IsSignal && kaon.IsSignal && v0.ReactionID == kaon.ReactionID} {
        if (IsSignal) {
            FillInfo_BeforeReaction(df, mass_sexaquark, Tree2Secondaries::Particle::Mass[EParticle::Proton], v0.ReactionID);
        } else {
            IsHybrid = (v0.IsSignal && !kaon.IsSignal) || (!v0.IsSignal && kaon.IsSignal) ||
                       (v0.IsSignal && kaon.IsSignal && v0.ReactionID != kaon.ReactionID) || v0.IsHybrid;
        }
    }

    // member vars //
    MC::V0 V0;
    MC::Track Kaon;
    double X{Const::DummyDouble};
    double Y{Const::DummyDouble};
    double Z{Const::DummyDouble};
    double AfterPx{Const::DummyDouble};
    double AfterPy{Const::DummyDouble};
    double AfterPz{Const::DummyDouble};
    double AfterE{Const::DummyDouble};
    bool IsSignal{false};
    bool IsHybrid{false};
};

}  // namespace Tree2Secondaries::MC

#endif  // T2S_MC_PARTICLES_HXX
