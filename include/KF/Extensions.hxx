#ifndef T2S_KF_EXTENSIONS_HXX
#define T2S_KF_EXTENSIONS_HXX

#include <cmath>

#include <KFParticle.hxx>
#include <KFParticle_Math.hxx>

#include "Math/Constants.hxx"
#include "Math/Math.hxx"

// Extensions of `KF::Particle` //

namespace Tree2Secondaries::KF {

struct alignas(T2S_SIMD_ALIGN) Track : ::KF::Particle {
    Track() = delete;
    Track(const ::KF::Particle& p, int idx_track) : ::KF::Particle{p}, idx{idx_track} {}
    Track(const ::KF::Vector<6>& p, const ::KF::SymMatrix<6>& cov, int charge, double mass, int idx_track)
        : ::KF::Particle{p, cov, charge, mass}, idx{idx_track} {}
    Track(const ::KF::Vector<7>& p, const ::KF::SymMatrix<7>& cov, int charge, int idx_track) : ::KF::Particle{p, cov, charge}, idx{idx_track} {}

    int idx{};
};

struct alignas(T2S_SIMD_ALIGN) V0 : ::KF::Particle {
    // constructors //
    V0() = delete;
    // -- to use when fitting V0 //
    V0(int idx_v0, Tree2Secondaries::EParticle v0_hypothesis, const KF::Track& neg, const KF::Track& pos)
        : Neg(neg), Pos(pos), idx{idx_v0}, hypothesis{v0_hypothesis} {}
    // -- to use when unpacking V0 //
    V0(int idx_v0, Tree2Secondaries::EParticle v0_hypothesis, const ::KF::Vector<7>& p, const ::KF::SymMatrix<7>& cov, const KF::Track& neg,
       const KF::Track& pos)
        : ::KF::Particle{p, cov, 0}, Neg(neg), Pos(pos), idx{idx_v0}, hypothesis{v0_hypothesis} {}

    // main //
    void DoFit(float bz, const double* mass = nullptr) {
        AddDaughter(Neg, bz);
        AddDaughter(Pos, bz);
        if (mass != nullptr) AddMassConstraint(*mass);
    }

    // utilities -- only usable after vertex fit //
    [[nodiscard]] ::KF::Vector<3> Neg_PCA_XYZ() const { return GetPCA(0).xyz; };
    [[nodiscard]] ::KF::Vector<3> Pos_PCA_XYZ() const { return GetPCA(1).xyz; };
    [[nodiscard]] ::KF::Vector<3> Neg_PCA_PxPyPz() const { return GetPCA(0).dir; };
    [[nodiscard]] ::KF::Vector<3> Pos_PCA_PxPyPz() const { return GetPCA(1).dir; };
    [[nodiscard]] double ArmenterosQt() const {
        return Tree2Secondaries::Math::ArmenterosQt(Px(), Py(), Pz(),  //
                                                    Neg_PCA_PxPyPz()[0], Neg_PCA_PxPyPz()[1], Neg_PCA_PxPyPz()[2]);
    }
    [[nodiscard]] double ArmenterosAlpha() const {
        return Tree2Secondaries::Math::ArmenterosAlpha(Px(), Py(), Pz(),                                               //
                                                       Neg_PCA_PxPyPz()[0], Neg_PCA_PxPyPz()[1], Neg_PCA_PxPyPz()[2],  //
                                                       Pos_PCA_PxPyPz()[0], Pos_PCA_PxPyPz()[1], Pos_PCA_PxPyPz()[2]);
    }

    // cuts -- only usable after vertex fit //
    [[nodiscard]] double AbsZ() const { return std::abs(Z()); }
    [[nodiscard]] double AbsEta() const { return std::abs(Eta()); }
    [[nodiscard]] double DCA_Daughters() const { return GetDCA(0, 1); }
    [[nodiscard]] double DCA_Neg_V0() const { return GetDCA(0); }
    [[nodiscard]] double DCA_Pos_V0() const { return GetDCA(1); }
    [[nodiscard]] double AbsArmQtOverAlpha() const { return ArmenterosQt() / std::abs(ArmenterosAlpha()); };
    [[nodiscard]] double DCA_Point(float x, float y, float z) const {
        return Tree2Secondaries::Math::FastDCALineVertex({Px(), Py(), Pz()}, {X(), Y(), Z()}, {x, y, z});
    }
    [[nodiscard]] double CPA_Point(float x, float y, float z) const {
        return Tree2Secondaries::Math::CosinePointingAngle({Px(), Py(), Pz()}, {X(), Y(), Z()}, {x, y, z});
    }

    // member vars //
    KF::Track Neg;
    KF::Track Pos;
    int idx{};
    Tree2Secondaries::EParticle hypothesis{};
};

struct alignas(T2S_SIMD_ALIGN) Sexaquark : ::KF::Particle {
    // constructors //
    Sexaquark() = delete;
    explicit Sexaquark(double nucleon_mass) : Nucleon_Mass{nucleon_mass} {};

    // utilities //
    [[nodiscard]] double E_MinusNucleon() const { return E() - Nucleon_Mass; };
    [[nodiscard]] double Mass_MinusNucleon() const {
        double mass2{E_MinusNucleon() * E_MinusNucleon() - P2()};
        if (mass2 > 0.) return std::sqrt(mass2);
        return Tree2Secondaries::Const::DummyDouble;
    };

    // cuts //
    [[nodiscard]] double AbsRapidity_MinusNucleon() const { return std::abs(std::log((E_MinusNucleon() + Pz()) / (E_MinusNucleon() - Pz())) / 2.); };
    [[nodiscard]] double CPA_Point(float x, float y, float z) const {
        return Tree2Secondaries::Math::CosinePointingAngle({Px(), Py(), Pz()}, {X(), Y(), Z()}, {x, y, z});
    }

    // member vars //
    double Nucleon_Mass{};
};

struct alignas(T2S_SIMD_ALIGN) ChannelA : Sexaquark {
    // constructors //
    ChannelA() = delete;
    ChannelA(const KF::V0& v0a, const KF::V0& v0b)  //
        : Sexaquark{Tree2Secondaries::Particle::Mass[EParticle::Neutron]}, V0A{v0a}, V0B{v0b} {};

    // main //
    void DoFit(float bz) {
        AddDaughter(V0A, bz);
        AddDaughter(V0B, bz);
    }

    // utilities //
    [[nodiscard]] ::KF::Vector<3> V0A_PCA_XYZ() const { return GetPCA(0).xyz; };
    [[nodiscard]] ::KF::Vector<3> V0B_PCA_XYZ() const { return GetPCA(1).xyz; };

    // cuts //
    [[nodiscard]] double DecayLength_V0A() const {
        return ::KF::Math::Norm(::KF::Vector<3>{V0A.X() - V0A_PCA_XYZ()[0], V0A.Y() - V0A_PCA_XYZ()[1], V0A.Z() - V0A_PCA_XYZ()[2]});
    };
    [[nodiscard]] double DecayLength_V0B() const {
        return ::KF::Math::Norm(::KF::Vector<3>{V0B.X() - V0B_PCA_XYZ()[0], V0B.Y() - V0B_PCA_XYZ()[1], V0B.Z() - V0B_PCA_XYZ()[2]});
    };
    [[nodiscard]] double DCA_btw_V0s() const { return GetDCA(0, 1); };
    [[nodiscard]] double DCA_V0A_wrt_SV() const { return GetDCA(0); };
    [[nodiscard]] double DCA_V0B_wrt_SV() const { return GetDCA(1); };
    [[nodiscard]] double DCA_V0ANeg_wrt_SV(float bz) const { return Tree2Secondaries::Math::FastDCAHelixVertex(V0A.Neg, X(), Y(), Z(), bz); };
    [[nodiscard]] double DCA_V0APos_wrt_SV(float bz) const { return Tree2Secondaries::Math::FastDCAHelixVertex(V0A.Pos, X(), Y(), Z(), bz); };
    [[nodiscard]] double DCA_V0BNeg_wrt_SV(float bz) const { return Tree2Secondaries::Math::FastDCAHelixVertex(V0B.Neg, X(), Y(), Z(), bz); };
    [[nodiscard]] double DCA_V0BPos_wrt_SV(float bz) const { return Tree2Secondaries::Math::FastDCAHelixVertex(V0B.Pos, X(), Y(), Z(), bz); };

    // daughters //
    KF::V0 V0A;
    KF::V0 V0B;
};

struct alignas(T2S_SIMD_ALIGN) ChannelD : Sexaquark {
    // constructors //
    ChannelD() = delete;
    ChannelD(const KF::V0& v0, const KF::Track& kaon) : Sexaquark{Tree2Secondaries::Particle::Mass[EParticle::Proton]}, V0{v0}, Kaon{kaon} {};

    // main //
    void DoFit(float bz) {
        AddDaughter(V0, bz);
        AddDaughter(Kaon, bz);
    }

    // utilities //
    [[nodiscard]] ::KF::Vector<3> V0_PCA_XYZ() const { return GetPCA(0).xyz; };
    [[nodiscard]] ::KF::Vector<3> Kaon_PCA_XYZ() const { return GetPCA(1).xyz; };
    [[nodiscard]] ::KF::Vector<3> Kaon_PCA_PxPyPz() const { return GetPCA(1).dir; };

    // cuts //
    [[nodiscard]] double DCA_btw_V0_Kaon() const { return GetDCA(0, 1); };
    [[nodiscard]] double DCA_V0_wrt_SV() const { return GetDCA(0); };
    [[nodiscard]] double DCA_Kaon_wrt_SV() const { return GetDCA(1); };
    [[nodiscard]] double DCA_V0Neg_wrt_SV(float bz) const { return Tree2Secondaries::Math::FastDCAHelixVertex(V0.Neg, X(), Y(), Z(), bz); };
    [[nodiscard]] double DCA_V0Pos_wrt_SV(float bz) const { return Tree2Secondaries::Math::FastDCAHelixVertex(V0.Pos, X(), Y(), Z(), bz); };
    [[nodiscard]] double DCA_btw_V0Neg_Kaon(float bz) const { return Tree2Secondaries::Math::FastDCAHelixHelix(V0.Neg, Kaon, bz); };
    [[nodiscard]] double DCA_btw_V0Pos_Kaon(float bz) const { return Tree2Secondaries::Math::FastDCAHelixHelix(V0.Pos, Kaon, bz); };

    // daughters //
    KF::V0 V0;
    KF::Track Kaon;
};

}  // namespace Tree2Secondaries::KF

#endif  // T2S_KF_EXTENSIONS_HXX
