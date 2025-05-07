#ifndef T2S_SECONDARY_NEUTRAL_HXX
#define T2S_SECONDARY_NEUTRAL_HXX

#include <memory>
#include <utility>

#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

#include "Math/Common.hxx"
#include "Math/VtxrResults.hxx"
#include "Secondary/Charged.hxx"
#include "Secondary/Reconstructed.hxx"

namespace Tree2Secondaries {

//
class Neutral final : public Reconstructed {
   public:
    Neutral(const Neutral&) = delete;
    Neutral(Neutral&&) noexcept = default;
    Neutral& operator=(const Neutral&) = delete;
    Neutral& operator=(Neutral&&) = default;
    ~Neutral() final = default;

    Neutral(int charge, ROOT::Math::XYZPoint vtx, ROOT::Math::PxPyPzMVector momentum) : Reconstructed{charge, std::move(vtx), std::move(momentum)} {}

    Neutral(size_t id, const std::shared_ptr<Charged>& neg, const std::shared_ptr<Charged>& pos, const PartPartResults& res)
        : Reconstructed{neg->Charge() + pos->Charge(),
                        {0.5 * (res.q.pca.X() + res.t.pca.X()),  //
                         0.5 * (res.q.pca.Y() + res.t.pca.Y()),  //
                         0.5 * (res.q.pca.Z() + res.t.pca.Z())},
                        neg->PxPyPzM(res.q.ds) + pos->PxPyPzM(res.t.ds)},  //
          fNeg{neg},
          fNegRes{res.q},
          fPos{pos},
          fPosRes{res.t},
          fIndex{id} {}

    [[nodiscard]] size_t Index() const override { return fIndex; }
    [[nodiscard]] size_t NegIndex() const { return fNeg->Index(); }
    [[nodiscard]] size_t PosIndex() const { return fPos->Index(); }

    [[nodiscard]] std::shared_ptr<Charged> Neg() const { return fNeg; }
    [[nodiscard]] std::shared_ptr<Charged> Pos() const { return fPos; }
    [[nodiscard]] ROOT::Math::XYZPoint NegPCA() const { return fNegRes.pca; }
    [[nodiscard]] ROOT::Math::XYZPoint PosPCA() const { return fPosRes.pca; }
    [[nodiscard]] ROOT::Math::PxPyPzEVector NegMom() const { return fNeg->PxPyPzE(fNegRes.ds); }
    [[nodiscard]] ROOT::Math::PxPyPzEVector PosMom() const { return fPos->PxPyPzE(fPosRes.ds); }
    [[nodiscard]] ROOT::Math::XYZVector NegDir() const { return NegMom().Vect(); }
    [[nodiscard]] ROOT::Math::XYZVector PosDir() const { return PosMom().Vect(); }

    // [[nodiscard]] double ProductionX() const { return fProdVtx.X(); }
    // [[nodiscard]] double ProductionY() const { return fProdVtx.Y(); }
    // [[nodiscard]] double ProductionZ() const { return fProdVtx.Z(); }
    // [[nodiscard]] const ROOT::Math::XYZPoint ProductionVertex() const { return fProdVtx; }
    // [[nodiscard]] double ProductionRadius() const { return fProdVtx.Rho(); }

    [[nodiscard]] ROOT::Math::XYZPoint DecayVertex() const { return fMeasuredVtx; }
    [[nodiscard]] double DecayX() const { return fMeasuredVtx.X(); }
    [[nodiscard]] double DecayY() const { return fMeasuredVtx.Y(); }
    [[nodiscard]] double DecayZ() const { return fMeasuredVtx.Z(); }
    [[nodiscard]] double DecayRadius() const { return fMeasuredVtx.Rho(); }

    // Parametrization of neutral particles: (line) //
    //   X(s)  = X0 + Px s                          //
    //   Y(s)  = Y0 + Py s                          //
    //   Z(s)  = Z0 + Pz s                          //
    [[nodiscard]] double X(double s) const override { return DecayX() + Px(s) * s; }
    [[nodiscard]] double Y(double s) const override { return DecayY() + Py(s) * s; }
    [[nodiscard]] double Z(double s) const override { return DecayZ() + Pz(s) * s; }
    [[nodiscard]] double Px(double s [[maybe_unused]]) const override { return Px(); }  // PENDING: wtf
    [[nodiscard]] double Py(double s [[maybe_unused]]) const override { return Py(); }  // PENDING: wtf
    [[nodiscard]] double Pz(double s [[maybe_unused]]) const override { return Pz(); }  // PENDING: wtf
    [[nodiscard]] double Px() const { return fMeasuredMom.Px(); }
    [[nodiscard]] double Py() const { return fMeasuredMom.Py(); }
    [[nodiscard]] double Pz() const { return fMeasuredMom.Pz(); }

    [[nodiscard]] double CPAwrt(const ROOT::Math::XYZPoint& v) const { return Math::CosinePointingAngle(PxPyPz(), DecayVertex(), v); }
    // [[nodiscard]] double DecayLength() const { return (fDecayVtx - fProdVtx).R(); }
    [[nodiscard]] double ArmenterosAlpha() const { return Math::ArmenterosAlpha(PxPyPz(), NegDir(), PosDir()); }
    [[nodiscard]] double ArmenterosQt() const { return Math::ArmenterosQt(PxPyPz(), NegDir()); }
    [[nodiscard]] double DCANegWrtV0() const { return (NegPCA() - DecayVertex()).R(); }
    [[nodiscard]] double DCAPosWrtV0() const { return (PosPCA() - DecayVertex()).R(); }
    [[nodiscard]] double DCAbtwDaughters() const { return (NegPCA() - PosPCA()).R(); }

    // based on true info //
    [[nodiscard]] bool DaughtersHaveSameMcMother() const {
        return fNeg->LinkedMc()->MotherMcEntry() != -1 && fNeg->LinkedMc()->MotherMcEntry() == fPos->LinkedMc()->MotherMcEntry();
    }
    [[nodiscard]] bool IsTrue(int pdg_code_v0, int pdg_code_neg, int pdg_code_pos) const {
        return LinkedMc()->PdgCode() == pdg_code_v0 && fNeg->LinkedMc()->PdgCode() == pdg_code_neg && fPos->LinkedMc()->PdgCode() == pdg_code_pos;
    }
    [[nodiscard]] bool IsSignal(int pdg_code_v0, int pdg_code_neg, int pdg_code_pos) const {
        return IsTrue(pdg_code_v0, pdg_code_neg, pdg_code_pos) && LinkedMc()->IsSignal();
    }
    [[nodiscard]] bool IsHybrid(int pdg_code_v0, int pdg_code_neg, int pdg_code_pos) const {
        return !IsSignal(pdg_code_v0, pdg_code_neg, pdg_code_pos) &&  //
               ((fNeg->LinkedMc()->IsSignal() && !fPos->LinkedMc()->IsSignal()) || (!fNeg->LinkedMc()->IsSignal() && fPos->LinkedMc()->IsSignal()));
    }

    // find true info //
    // if (IsMC() && neg.LinkedMc()->MotherMcEntry() != -1 && neg.LinkedMc()->MotherMcEntry() == neg.LinkedMc()->MotherMcEntry()) {
    // V0.SetLinkedMc(fMCParticles[neg.LinkedMc()->MotherMcEntry()]);
    // }

   protected:
    std::shared_ptr<Charged> fNeg;
    VtxrResults fNegRes;
    std::shared_ptr<Charged> fPos;
    VtxrResults fPosRes;
    size_t fIndex{0};
};

}  // namespace Tree2Secondaries

#endif  // T2S_SECONDARY_NEUTRAL_HXX
