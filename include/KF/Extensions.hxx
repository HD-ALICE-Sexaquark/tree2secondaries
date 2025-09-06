#ifndef T2S_KF_EXTENSIONS_HXX
#define T2S_KF_EXTENSIONS_HXX

#include <cmath>

#include <KFParticle.hxx>
#include <KFParticle_Math.hxx>

#include "Math/Constants.hxx"
#include "Math/Math.hxx"

// Extensions of `KF::Particle` //

namespace KF {

struct alignas(32) Track : Particle {
    // constructors //
    Track() = delete;
    Track(const Particle& p, int idx_track) : Particle{p}, idx{idx_track} {}
    Track(const Vector<6>& p, const SymMatrix<6>& cov, int charge, double mass, int idx_track) : Particle{p, cov, charge, mass}, idx{idx_track} {}
    Track(const Vector<7>& p, const SymMatrix<7>& cov, int charge, int idx_track) : Particle{p, cov, charge}, idx{idx_track} {}

    // member vars //
    int idx{};
};

struct alignas(32) V0 : Particle {
    // constructors //
    V0() = delete;
    // -- to use when fitting V0 //
    V0(int idx_v0, Tree2Secondaries::EParticle v0_hypothesis, const KF::Track& neg, const KF::Track& pos)
        : Neg(neg), Pos(pos), idx{idx_v0}, hypothesis{v0_hypothesis} {}
    // -- to use when unpacking V0 //
    V0(int idx_v0, Tree2Secondaries::EParticle v0_hypothesis, const Vector<7>& p, const SymMatrix<7>& cov, const KF::Track& neg, const KF::Track& pos)
        : Particle{p, cov, 0}, Neg(neg), Pos(pos), idx{idx_v0}, hypothesis{v0_hypothesis} {}

    // main //
    void DoFit(float bz, const double* mass = nullptr) {
        AddDaughter(Neg, bz);
        AddDaughter(Pos, bz);
        if (mass != nullptr) AddMassConstraint(*mass);
    }

    // utilities -- only usable after vertex fit //
    KF::Vector<3> Neg_PCA_XYZ() const { return GetPCA(0).xyz; };
    KF::Vector<3> Pos_PCA_XYZ() const { return GetPCA(1).xyz; };
    KF::Vector<3> Neg_PCA_PxPyPz() const { return GetPCA(0).dir; };
    KF::Vector<3> Pos_PCA_PxPyPz() const { return GetPCA(1).dir; };
    double ArmenterosQt() const {
        return Tree2Secondaries::Math::ArmenterosQt(Px(), Py(), Pz(),  //
                                                    Neg_PCA_PxPyPz()[0], Neg_PCA_PxPyPz()[1], Neg_PCA_PxPyPz()[2]);
    }
    double ArmenterosAlpha() const {
        return Tree2Secondaries::Math::ArmenterosAlpha(Px(), Py(), Pz(),                                               //
                                                       Neg_PCA_PxPyPz()[0], Neg_PCA_PxPyPz()[1], Neg_PCA_PxPyPz()[2],  //
                                                       Pos_PCA_PxPyPz()[0], Pos_PCA_PxPyPz()[1], Pos_PCA_PxPyPz()[2]);
    }

    // cuts -- only usable after vertex fit //
    double AbsZ() const { return std::abs(Z()); }
    double AbsEta() const { return std::abs(Eta()); }
    double DCA_Daughters() const { return GetDCA(0, 1); }
    double DCA_Neg_V0() const { return GetDCA(0); }
    double DCA_Pos_V0() const { return GetDCA(1); }
    double AbsArmQtOverAlpha() const { return ArmenterosQt() / std::abs(ArmenterosAlpha()); };
    double DCA_Point(float x, float y, float z) const {
        return Tree2Secondaries::Math::FastDCALineVertex({Px(), Py(), Pz()}, {X(), Y(), Z()}, {x, y, z});
    }
    double CPA_Point(float x, float y, float z) const {
        return Tree2Secondaries::Math::CosinePointingAngle({Px(), Py(), Pz()}, {X(), Y(), Z()}, {x, y, z});
    }

    // member vars //
    KF::Track Neg;
    KF::Track Pos;
    int idx{};
    Tree2Secondaries::EParticle hypothesis{};
};

struct alignas(32) Sexaquark : Particle {
    // constructors //
    Sexaquark() = delete;
    explicit Sexaquark(double nucleon_mass) : Nucleon_Mass{nucleon_mass} {};

    // utilities //
    double E_MinusNucleon() const { return E() - Nucleon_Mass; };
    double Mass_MinusNucleon() const {
        double mass2{E_MinusNucleon() * E_MinusNucleon() - P2()};
        if (mass2 > 0.) return std::sqrt(mass2);
        return Tree2Secondaries::Const::DummyDouble;
    };

    // cuts //
    double AbsRapidity_MinusNucleon() const { return std::abs(std::log((E_MinusNucleon() + Pz()) / (E_MinusNucleon() - Pz())) / 2.); };
    double CPA_Point(float x, float y, float z) const {
        return Tree2Secondaries::Math::CosinePointingAngle({Px(), Py(), Pz()}, {X(), Y(), Z()}, {x, y, z});
    }

    // member vars //
    double Nucleon_Mass{};
};

struct alignas(32) ChannelA : Sexaquark {
    // constructors //
    ChannelA() = delete;
    ChannelA(const KF::V0& v0a, const KF::V0& v0b)  //
        : Sexaquark{Tree2Secondaries::Particle::Mass[Tree2Secondaries::EParticle::Neutron]}, V0A{v0a}, V0B{v0b} {};

    // main //
    void DoFit(float bz) {
        AddDaughter(V0A, bz);
        AddDaughter(V0B, bz);
    }

    // utilities //
    KF::Vector<3> V0A_PCA_XYZ() const { return GetPCA(0).xyz; };
    KF::Vector<3> V0B_PCA_XYZ() const { return GetPCA(1).xyz; };

    // cuts //
    double DecayLength_V0A() const {
        return KF::Math::Norm(KF::Vector<3>{V0A.X() - V0A_PCA_XYZ()[0], V0A.Y() - V0A_PCA_XYZ()[1], V0A.Z() - V0A_PCA_XYZ()[2]});
    };
    double DecayLength_V0B() const {
        return KF::Math::Norm(KF::Vector<3>{V0B.X() - V0B_PCA_XYZ()[0], V0B.Y() - V0B_PCA_XYZ()[1], V0B.Z() - V0B_PCA_XYZ()[2]});
    };
    double DCA_btw_V0s() const { return GetDCA(0, 1); };
    double DCA_V0A_wrt_SV() const { return GetDCA(0); };
    double DCA_V0B_wrt_SV() const { return GetDCA(1); };
    double DCA_V0ANeg_wrt_SV(float bz) const { return Tree2Secondaries::Math::FastDCAHelixVertex(V0A.Neg, X(), Y(), Z(), bz); };
    double DCA_V0APos_wrt_SV(float bz) const { return Tree2Secondaries::Math::FastDCAHelixVertex(V0A.Pos, X(), Y(), Z(), bz); };
    double DCA_V0BNeg_wrt_SV(float bz) const { return Tree2Secondaries::Math::FastDCAHelixVertex(V0B.Neg, X(), Y(), Z(), bz); };
    double DCA_V0BPos_wrt_SV(float bz) const { return Tree2Secondaries::Math::FastDCAHelixVertex(V0B.Pos, X(), Y(), Z(), bz); };

    // daughters //
    KF::V0 V0A;
    KF::V0 V0B;
};

struct alignas(32) ChannelD : Sexaquark {
    // constructors //
    ChannelD() = delete;
    ChannelD(const KF::V0& v0, const KF::Track& kaon)
        : Sexaquark{Tree2Secondaries::Particle::Mass[Tree2Secondaries::EParticle::Proton]}, V0{v0}, Kaon{kaon} {};

    // main //
    void DoFit(float bz) {
        AddDaughter(V0, bz);
        AddDaughter(Kaon, bz);
    }

    // utilities //
    KF::Vector<3> V0_PCA_XYZ() const { return GetPCA(0).xyz; };
    KF::Vector<3> Kaon_PCA_XYZ() const { return GetPCA(1).xyz; };
    KF::Vector<3> Kaon_PCA_PxPyPz() const { return GetPCA(1).dir; };

    // cuts //
    double DCA_btw_V0_Kaon() const { return GetDCA(0, 1); };
    double DCA_V0_wrt_SV() const { return GetDCA(0); };
    double DCA_Kaon_wrt_SV() const { return GetDCA(1); };
    double DCA_V0Neg_wrt_SV(float bz) const { return Tree2Secondaries::Math::FastDCAHelixVertex(V0.Neg, X(), Y(), Z(), bz); };
    double DCA_V0Pos_wrt_SV(float bz) const { return Tree2Secondaries::Math::FastDCAHelixVertex(V0.Pos, X(), Y(), Z(), bz); };
    double DCA_btw_V0Neg_Kaon(float bz) const { return Tree2Secondaries::Math::FastDCAHelixHelix(V0.Neg, Kaon, bz); };
    double DCA_btw_V0Pos_Kaon(float bz) const { return Tree2Secondaries::Math::FastDCAHelixHelix(V0.Pos, Kaon, bz); };

    // daughters //
    KF::V0 V0;
    KF::Track Kaon;
};

struct alignas(32) ChannelE : Sexaquark {
    // PENDING
};

struct alignas(32) ChannelH : Sexaquark {
    // PENDING
};

}  // namespace KF

#endif  // T2S_KF_EXTENSIONS_HXX
