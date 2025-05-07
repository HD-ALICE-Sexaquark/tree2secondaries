#ifndef T2S_SECONDARY_CHANNEL_E_HXX
#define T2S_SECONDARY_CHANNEL_E_HXX

#include <memory>

#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

#include "Math/VtxrResults.hxx"
#include "Secondary/Charged.hxx"
#include "Secondary/Neutral.hxx"
#include "Secondary/Sexaquark.hxx"

namespace Tree2Secondaries {

class TypeE final : public Sexaquark {
   public:
    TypeE(const TypeE&) = delete;
    TypeE(TypeE&&) noexcept = default;
    TypeE& operator=(const TypeE&) = delete;
    TypeE& operator=(TypeE&&) = default;
    ~TypeE() final = default;

   protected:
};

}  // namespace Tree2Secondaries

#endif  // T2S_SECONDARY_CHANNEL_E_HXX
