#pragma once

#include <KFParticle.hxx>
#include <KFParticle_Math.hxx>

#include "Math/Constants.hxx"
#include "View/Reconstructed/ViewTrack.hxx"

namespace Tree2Secondaries::Fit {

struct alignas(T2S_SIMD_ALIGN) Track : KF::Particle {
    Track(const View::Rec::Track& view, double mass)
        : KF::Particle{view.State_NoE(), view.CovMatrix_NoE(), view.Charge(), mass},  //
          View{view} {}

    View::Rec::Track View;
};

}  // namespace Tree2Secondaries::Fit
