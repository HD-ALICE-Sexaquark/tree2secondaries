#ifndef T2S_SECONDARY_CHANNEL_D_HXX
#define T2S_SECONDARY_CHANNEL_D_HXX

#include <memory>

#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

#include "Math/Vertexer.hxx"
#include "Math/VtxrResults.hxx"
#include "Secondary/Charged.hxx"
#include "Secondary/Neutral.hxx"
#include "Secondary/Sexaquark.hxx"

namespace Tree2Secondaries {

class TypeD final : public Sexaquark {
   public:
    TypeD(const TypeD&) = delete;
    TypeD(TypeD&&) noexcept = default;
    TypeD& operator=(const TypeD&) = delete;
    TypeD& operator=(TypeD&&) = default;
    ~TypeD() final = default;

    TypeD(size_t idx, const std::shared_ptr<Neutral>& lambda, const std::shared_ptr<Charged>& kaon, const PartPartResults& res, double nucleon_mass)
        : Sexaquark{idx,                                      //
                    {0.5 * (res.q.pca.X() + res.t.pca.X()),   //
                     0.5 * (res.q.pca.Y() + res.t.pca.Y()),   //
                     0.5 * (res.q.pca.Z() + res.t.pca.Z())},  //
                    lambda->PxPyPzE(res.q.ds) + kaon->PxPyPzE(res.t.ds) - ROOT::Math::PxPyPzEVector{0., 0., 0., nucleon_mass},
                    nucleon_mass},  //
          fLambda{lambda},
          fLambdaRes{res.q},
          fKaon{kaon},
          fKaonRes{res.t} {}

    [[nodiscard]] std::shared_ptr<Neutral> Lambda() const { return fLambda; }
    [[nodiscard]] std::shared_ptr<Charged> Kaon() const { return fKaon; }
    [[nodiscard]] ROOT::Math::XYZPoint LambdaPCA() const { return fLambdaRes.pca; }
    [[nodiscard]] ROOT::Math::XYZPoint KaonPCA() const { return fKaonRes.pca; }
    [[nodiscard]] ROOT::Math::PxPyPzEVector LambdaMom() const { return fLambda->PxPyPzE(fLambdaRes.ds); }
    [[nodiscard]] ROOT::Math::PxPyPzEVector KaonMom() const { return fKaon->PxPyPzE(fKaonRes.ds); }
    [[nodiscard]] ROOT::Math::XYZVector LambdaDir() const { return LambdaMom().Vect(); }
    [[nodiscard]] ROOT::Math::XYZVector KaonDir() const { return KaonMom().Vect(); }

    [[nodiscard]] double DecayLengthLa() const { return (SecondaryVertex() - fLambda->DecayVertex()).R(); }
    [[nodiscard]] double DCALaSV() const { return (SecondaryVertex() - LambdaPCA()).R(); }
    [[nodiscard]] double DCAKaSV() const { return (SecondaryVertex() - KaonPCA()).R(); }
    [[nodiscard]] double DCAKaLa() const { return (LambdaPCA() - KaonPCA()).R(); }

    [[nodiscard]] double DCALaNegSV() const {
        VtxrResults res{Vertexer::MinimizeDistanceHelixVertex(*fLambda->Neg(), SecondaryVertex())};
        return (SecondaryVertex() - res.pca).R();
    }
    [[nodiscard]] double DCALaPosSV() const {
        VtxrResults res{Vertexer::MinimizeDistanceHelixVertex(*fLambda->Pos(), SecondaryVertex())};
        return (SecondaryVertex() - res.pca).R();
    }
    /*
    [[nodiscard]] double DCALaNegKa() const {
        PartPartResults res{Vertexer::MinimizeDistanceHelixHelix(*fLambda->Neg(), *fKaon)};
        return (res.q.pca - res.t.pca).R();
    }
    */

   protected:
    std::shared_ptr<Neutral> fLambda;
    VtxrResults fLambdaRes;
    std::shared_ptr<Charged> fKaon;
    VtxrResults fKaonRes;
};

}  // namespace Tree2Secondaries

#endif  // T2S_SECONDARY_CHANNEL_D_HXX
