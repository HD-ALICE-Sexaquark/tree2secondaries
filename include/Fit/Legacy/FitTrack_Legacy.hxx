#pragma once

#include <KFPTrack.h>
#include <KFParticle.h>

#include "Math/Constants.hxx"
#include "View/Reconstructed/Legacy/ViewTrack_Legacy.hxx"

namespace Tree2Secondaries::Fit {

struct alignas(T2S_SIMD_ALIGN) Track : KFParticle {
    Track(const View::Rec::Track& view, double mass)  //
        : View{view} {
        Create(view.State_NoE().data(), view.CovMatrix_NoE().data(), view.Charge(), static_cast<float>(mass));
    }

    View::Rec::Track View;
};

}  // namespace Tree2Secondaries::Fit
