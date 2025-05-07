#ifndef T2S_SECONDARY_CHANNEL_A_HXX
#define T2S_SECONDARY_CHANNEL_A_HXX

#include <memory>

#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

#include "Math/Vertexer.hxx"
#include "Math/VtxrResults.hxx"
#include "Secondary/Neutral.hxx"
#include "Secondary/Sexaquark.hxx"

namespace Tree2Secondaries {

class TypeA final : public Sexaquark {
   public:
    TypeA(const TypeA&) = delete;
    TypeA(TypeA&&) noexcept = default;
    TypeA& operator=(const TypeA&) = delete;
    TypeA& operator=(TypeA&&) = default;
    ~TypeA() final = default;

    TypeA(size_t idx, const std::shared_ptr<Neutral>& lambda, const std::shared_ptr<Neutral>& k0s, const PartPartResults& res, double nucleon_mass)
        : Sexaquark{idx,                                      //
                    {0.5 * (res.q.pca.X() + res.t.pca.X()),   //
                     0.5 * (res.q.pca.Y() + res.t.pca.Y()),   //
                     0.5 * (res.q.pca.Z() + res.t.pca.Z())},  //
                    lambda->PxPyPzE(res.q.ds) + k0s->PxPyPzE(res.t.ds) - ROOT::Math::PxPyPzEVector{0., 0., 0., nucleon_mass},
                    nucleon_mass},  //
          fLambda{lambda},
          fLambdaRes{res.q},
          fKaonZero{k0s},
          fKaonZeroRes{res.t} {}

    [[nodiscard]] std::shared_ptr<Neutral> Lambda() const { return fLambda; }
    [[nodiscard]] std::shared_ptr<Neutral> KaonZero() const { return fKaonZero; }
    [[nodiscard]] ROOT::Math::XYZPoint LambdaPCA() const { return fLambdaRes.pca; }
    [[nodiscard]] ROOT::Math::XYZPoint KaonZeroPCA() const { return fKaonZeroRes.pca; }
    [[nodiscard]] ROOT::Math::PxPyPzEVector LambdaMom() const { return fLambda->PxPyPzE(fLambdaRes.ds); }
    [[nodiscard]] ROOT::Math::PxPyPzEVector KaonZeroMom() const { return fKaonZero->PxPyPzE(fKaonZeroRes.ds); }
    [[nodiscard]] ROOT::Math::XYZVector LambdaDir() const { return LambdaMom().Vect(); }
    [[nodiscard]] ROOT::Math::XYZVector KaonZeroDir() const { return KaonZeroMom().Vect(); }

    [[nodiscard]] double DecayLengthLa() const { return (SecondaryVertex() - fLambda->DecayVertex()).R(); }
    [[nodiscard]] double DecayLengthK0() const { return (SecondaryVertex() - fKaonZero->DecayVertex()).R(); }
    [[nodiscard]] double DCALaSV() const { return (SecondaryVertex() - LambdaPCA()).R(); }
    [[nodiscard]] double DCAK0SV() const { return (SecondaryVertex() - KaonZeroPCA()).R(); }
    [[nodiscard]] double DCAbtwV0s() const { return (LambdaPCA() - KaonZeroPCA()).R(); }

    [[nodiscard]] double DCALaNegSV() const {
        VtxrResults res{Vertexer::MinimizeDistanceHelixVertex(*fLambda->Neg(), SecondaryVertex())};
        return (SecondaryVertex() - res.pca).R();
    }
    [[nodiscard]] double DCALaPosSV() const {
        VtxrResults res{Vertexer::MinimizeDistanceHelixVertex(*fLambda->Pos(), SecondaryVertex())};
        return (SecondaryVertex() - res.pca).R();
    }
    [[nodiscard]] double DCAK0NegSV() const {
        VtxrResults res{Vertexer::MinimizeDistanceHelixVertex(*fKaonZero->Neg(), SecondaryVertex())};
        return (SecondaryVertex() - res.pca).R();
    }
    [[nodiscard]] double DCAK0PosSV() const {
        VtxrResults res{Vertexer::MinimizeDistanceHelixVertex(*fKaonZero->Pos(), SecondaryVertex())};
        return (SecondaryVertex() - res.pca).R();
    }

   protected:
    std::shared_ptr<Neutral> fLambda;
    VtxrResults fLambdaRes;
    std::shared_ptr<Neutral> fKaonZero;
    VtxrResults fKaonZeroRes;
};

}  // namespace Tree2Secondaries

#endif  // T2S_SECONDARY_CHANNEL_A_HXX
