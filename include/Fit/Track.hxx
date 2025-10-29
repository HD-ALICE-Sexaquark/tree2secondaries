#pragma once

#include <KFParticle.hxx>
#include <KFParticle_Math.hxx>

#include "Math/Constants.hxx"

namespace Tree2Secondaries::Fit {

struct alignas(T2S_SIMD_ALIGN) Track : KF::Particle {
    Track() = delete;
    Track(const KF::Particle& p, int idx)  //
        : KF::Particle{p}, Index{idx} {}
    Track(const KF::Vector<6>& p, const KF::SymMatrix<6>& cov, int charge, double mass, int idx) : KF::Particle{p, cov, charge, mass}, Index{idx} {}
    Track(const KF::Vector<7>& p, const KF::SymMatrix<7>& cov, int charge, int idx)  //
        : KF::Particle{p, cov, charge}, Index{idx} {}

    int Index{};
};

}  // namespace Tree2Secondaries::Fit
