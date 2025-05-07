#ifndef T2S_SECONDARY_CHANNEL_H_HXX
#define T2S_SECONDARY_CHANNEL_H_HXX

#include <memory>

#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

#include "Math/VtxrResults.hxx"
#include "Secondary/Charged.hxx"
#include "Secondary/Sexaquark.hxx"

namespace Tree2Secondaries {

class TypeH final : public Sexaquark {
   public:
    TypeH(const TypeH&) = delete;
    TypeH(TypeH&&) noexcept = default;
    TypeH& operator=(const TypeH&) = delete;
    TypeH& operator=(TypeH&&) = default;
    ~TypeH() final = default;

   protected:
};

}  // namespace Tree2Secondaries

#endif  // T2S_SECONDARY_CHANNEL_H_HXX
