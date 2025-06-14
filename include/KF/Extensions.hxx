#ifndef T2S_KF_EXTENSIONS_HXX
#define T2S_KF_EXTENSIONS_HXX

#include <array>
#include <cmath>
#include <cstddef>

#include <KFParticle.hxx>
#include <KFParticle_Math.hxx>

#include "Math/Constants.hxx"
#include "Math/Math.hxx"

// Extensions of `KF::Particle` //

namespace KF {

struct alignas(32) Track : Particle {
    // constructors //
    Track() = default;
    Track(const Particle& p, size_t idx_track) : Particle{p}, idx{idx_track} {}
    Track(const Vector<7>& p, const SymMatrix<7>& cov, int charge, size_t idx_track) : Particle{p, cov, charge}, idx{idx_track} {}

    // member vars //
    size_t idx{};
};

struct alignas(32) V0 : Particle {
    // constructors //
    V0() = default;
    V0(const Vector<7>& p, const SymMatrix<7>& cov, int charge, size_t idx_v0, size_t idx_v0_neg, size_t idx_v0_pos, int pdg_code_v0_hyp)
        : Particle{p, cov, charge}, idx{idx_v0}, idx_neg{idx_v0_neg}, idx_pos{idx_v0_pos}, pdg_code_hyp{pdg_code_v0_hyp} {}

    // utilities //
    void SetIndices(const std::array<size_t, 3>& indices) {
        idx = indices[0];
        idx_neg = indices[1];
        idx_pos = indices[2];
    }
    KF::Vector<3> PCA_Neg() const { return GetPCA(0).xyz; };
    KF::Vector<3> PCA_Pos() const { return GetPCA(1).xyz; };
    KF::Vector<3> Mom_Neg() const { return GetPCA(0).dir; };
    KF::Vector<3> Mom_Pos() const { return GetPCA(1).dir; };
    double ArmenterosQt() const {
        return Tree2Secondaries::Math::ArmenterosQt(Px(), Py(), Pz(),  //
                                                    Mom_Neg()[0], Mom_Neg()[1], Mom_Neg()[2]);
    }
    double ArmenterosAlpha() const {
        return Tree2Secondaries::Math::ArmenterosAlpha(Px(), Py(), Pz(),                          //
                                                       Mom_Neg()[0], Mom_Neg()[1], Mom_Neg()[2],  //
                                                       Mom_Pos()[0], Mom_Pos()[1], Mom_Pos()[2]);
    }

    // cuts //
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
    size_t idx{};
    size_t idx_neg{};
    size_t idx_pos{};
    double Neg_Energy{};
    double Pos_Energy{};
    int pdg_code_hyp{};  // hypothesis
};

struct alignas(32) Sexaquark : Particle {
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
    size_t idx{};
    double Nucleon_Mass{};
};

struct alignas(32) ChannelA : Sexaquark {
    // utilities //
    void SetIndices(const std::array<size_t, 7>& indices) {
        idx = indices[0];
        idx_lambda = indices[1];
        idx_lambda_neg = indices[2];
        idx_lambda_pos = indices[3];
        idx_k0s = indices[4];
        idx_k0s_neg = indices[5];
        idx_k0s_pos = indices[6];
    }
    void SetAdditionalInfo_V0A(const KF::Vector<3>& decay_vtx, double energy) {
        V0A_DecayVtx = decay_vtx;
        V0A_Energy = energy;
    }
    void SetAdditionalInfo_V0B(const KF::Vector<3>& decay_vtx, double energy) {
        V0B_DecayVtx = decay_vtx;
        V0B_Energy = energy;
    }
    KF::Vector<3> PCA_V0A() const { return GetPCA(0).xyz; };
    KF::Vector<3> PCA_V0B() const { return GetPCA(1).xyz; };
    KF::Vector<3> Mom_V0A() const { return GetPCA(0).dir; };
    KF::Vector<3> Mom_V0B() const { return GetPCA(1).dir; };

    // cuts //
    double DecayLength_V0A() const {
        KF::Vector<3> diff{};
        for (int i{0}; i < 3; ++i) diff[i] += V0A_DecayVtx[i] - PCA_V0A()[i];
        return KF::Math::SquaredNorm(diff);
    };
    double DecayLength_V0B() const {
        KF::Vector<3> diff{};
        for (int i{0}; i < 3; ++i) diff[i] += V0B_DecayVtx[i] - PCA_V0B()[i];
        return KF::Math::SquaredNorm(diff);
    };
    double DCA_V0s() const { return GetDCA(0, 1); };
    double DCA_V0A() const { return GetDCA(0); };
    double DCA_V0B() const { return GetDCA(1); };
    double DCA_V0A_Neg() const;  // PENDING: not trivial
    double DCA_V0A_Pos() const;  // PENDING: not trivial
    double DCA_V0B_Neg() const;  // PENDING: not trivial
    double DCA_V0B_Pos() const;  // PENDING: not trivial

    // member vars //
    KF::Vector<3> V0A_DecayVtx{};
    KF::Vector<3> V0B_DecayVtx{};
    double V0A_Energy{};
    double V0B_Energy{};
    size_t idx_lambda{};
    size_t idx_lambda_neg{};
    size_t idx_lambda_pos{};
    size_t idx_k0s{};
    size_t idx_k0s_neg{};
    size_t idx_k0s_pos{};
};

struct alignas(32) ChannelD : Sexaquark {
    // utilities //
    void SetIndices(const std::array<size_t, 5>& indices) {
        idx = indices[0];
        idx_lambda = indices[1];
        idx_lambda_neg = indices[2];
        idx_lambda_pos = indices[3];
        idx_kaon = indices[4];
    }
    void SetAdditionalInfo_V0(const KF::Vector<3>& decay_vtx, double energy) {
        V0_DecayVtx = decay_vtx;
        V0_Energy = energy;
    }
    KF::Vector<3> PCA_V0() const { return GetPCA(0).xyz; };
    KF::Vector<3> PCA_Kaon() const { return GetPCA(1).xyz; };
    KF::Vector<3> Mom_V0() const { return GetPCA(0).dir; };
    KF::Vector<3> Mom_Kaon() const { return GetPCA(1).dir; };

    // cuts //
    double DCA_V0_Kaon() const { return GetDCA(0, 1); };
    double DCA_V0() const { return GetDCA(0); };
    double DCA_Kaon() const { return GetDCA(1); };
    double DCA_V0_Neg() const;       // PENDING: not trivial
    double DCA_V0_Pos() const;       // PENDING: not trivial
    double DCA_V0_Neg_Kaon() const;  // PENDING: not trivial
    double DCA_V0_Pos_Kaon() const;  // PENDING: not trivial

    // member vars //
    KF::Vector<3> V0_DecayVtx{};
    double V0_Energy{};
    double Kaon_Energy{};
    size_t idx_lambda{};
    size_t idx_lambda_neg{};
    size_t idx_lambda_pos{};
    size_t idx_kaon{};
};

struct alignas(32) ChannelE : Sexaquark {
    size_t idx_lambda{};
    size_t idx_lambda_neg{};
    size_t idx_lambda_pos{};
    size_t idx_kaon{};
    size_t idx_pim{};
    size_t idx_pip{};
};

struct alignas(32) ChannelH : Sexaquark {
    size_t idx_kaon1{};
    size_t idx_kaon2{};
};

}  // namespace KF

#endif  // T2S_KF_EXTENSIONS_HXX
