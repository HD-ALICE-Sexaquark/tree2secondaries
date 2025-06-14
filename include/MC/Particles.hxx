#ifndef T2S_MC_PARTICLES_HXX
#define T2S_MC_PARTICLES_HXX

#include <cmath>

#include <KFParticle.hxx>
#include <KFParticle_Math.hxx>

#include "DataFormats/Events.hxx"
#include "Math/Constants.hxx"

namespace Tree2Secondaries::MC {

struct alignas(32) Particle {
    // utilities //
    void FillBasic(const Events::MC& soa, long mc_idx) {
        if (mc_idx >= 0) {
            Entry = mc_idx;
            X = soa.X->at(mc_idx);
            Y = soa.Y->at(mc_idx);
            Z = soa.Z->at(mc_idx);
            Px = soa.Px->at(mc_idx);
            Py = soa.Py->at(mc_idx);
            Pz = soa.Pz->at(mc_idx);
            Energy = soa.E->at(mc_idx);

            PdgCode = soa.PdgCode->at(mc_idx);
            MotherEntry = soa.MotherEntry->at(mc_idx);
            if (MotherEntry >= 0) PdgCode_Mother = soa.PdgCode->at(MotherEntry);
        }
    }

    // member vars //
    long Entry{Const::DummyInt};
    long MotherEntry{Const::DummyInt};
    float X{Const::DummyFloat};
    float Y{Const::DummyFloat};
    float Z{Const::DummyFloat};
    float Px{Const::DummyFloat};
    float Py{Const::DummyFloat};
    float Pz{Const::DummyFloat};
    float Energy{Const::DummyFloat};
    int PdgCode{Const::DummyInt};
    int PdgCode_Mother{Const::DummyInt};
};

struct alignas(32) Track : Particle {
    // constructors //
    Track() = default;
    Track(const Events::MC& soa, long mc_idx, Tree2Secondaries::PdgCode hypothesis) { Init(soa, mc_idx, hypothesis); }

    // utilities //
    void Init(const Events::MC& soa, long mc_idx, Tree2Secondaries::PdgCode hypothesis) {
        FillBasic(soa, mc_idx);
        FillDerived(soa, mc_idx, hypothesis);
    }
    void FillDerived(const Events::MC& soa, long mc_idx, Tree2Secondaries::PdgCode hypothesis) {
        if (mc_idx >= 0) {
            IsTrue = PdgCode == static_cast<int>(hypothesis);
            IsSignal = IsTrue && soa.Generator->at(mc_idx) == 2;
            IsSecondary = soa.IsSecFromMat->at(mc_idx) || soa.IsSecFromWeak->at(mc_idx) || IsSignal;
            if (MotherEntry >= 0) {
                if (IsSignal) ReactionID = static_cast<int>(soa.Status->at(MotherEntry));
                GrandMotherEntry = soa.MotherEntry->at(MotherEntry);
                if (GrandMotherEntry >= 0) PdgCode_GrandMother = soa.PdgCode->at(GrandMotherEntry);
            }
        }
    }

    // member vars //
    long GrandMotherEntry{Const::DummyInt};
    int PdgCode_GrandMother{Const::DummyInt};
    int ReactionID{Const::DummyInt};
    bool IsTrue{false};
    bool IsSignal{false};
    bool IsSecondary{false};
};

struct alignas(32) V0 : Particle {
    // constructors //
    V0() = default;
    V0(const Events::MC& soa, long mc_neg, long mc_pos, Tree2Secondaries::PdgCode hyp_v0, Tree2Secondaries::PdgCode hyp_neg,
       Tree2Secondaries::PdgCode hyp_pos) {
        if (mc_neg < 0 || mc_pos < 0) return;  // protection
        neg.Init(soa, mc_neg, hyp_neg);
        pos.Init(soa, mc_pos, hyp_pos);

        long mc_mother_neg{soa.MotherEntry->at(mc_neg)};
        long mc_mother_pos{soa.MotherEntry->at(mc_pos)};
        if (mc_mother_neg >= 0 && mc_mother_neg == mc_mother_pos) {
            long mc_v0{mc_mother_neg};
            FillBasic(soa, mc_v0);

            // fill derived props //
            IsTrue = PdgCode == static_cast<int>(hyp_v0) && neg.PdgCode == static_cast<int>(hyp_neg) && pos.PdgCode == static_cast<int>(hyp_pos);
            IsSignal = IsTrue && soa.Generator->at(mc_v0) == 2 && neg.ReactionID == pos.ReactionID;
            IsSecondary = soa.IsSecFromMat->at(mc_v0) || soa.IsSecFromWeak->at(mc_v0) || IsSignal;
            if (IsSignal)
                ReactionID = static_cast<int>(soa.Status->at(mc_v0));
            else
                IsHybrid = (neg.IsSignal && !pos.IsSignal) || (!neg.IsSignal && pos.IsSignal);
        }
    }

    // utilities //
    float DecayX() const { return neg.Entry >= 0 ? neg.X : pos.X; };
    float DecayY() const { return neg.Entry >= 0 ? neg.Y : pos.Y; };
    float DecayZ() const { return neg.Entry >= 0 ? neg.Z : pos.Z; };

    // member vars //
    MC::Track neg;
    MC::Track pos;
    int ReactionID{Const::DummyInt};
    bool IsHybrid{false};
    bool IsTrue{false};
    bool IsSignal{false};
    bool IsSecondary{false};
};

}  // namespace Tree2Secondaries::MC

#endif  // T2S_MC_PARTICLES_HXX
