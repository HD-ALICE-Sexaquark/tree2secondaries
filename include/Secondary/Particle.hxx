#ifndef T2S_SECONDARY_PARTICLE_HXX
#define T2S_SECONDARY_PARTICLE_HXX

#include <cstdlib>

namespace Tree2Secondaries {

// Abstract base class for all particles.
class Particle {
   public:
    virtual ~Particle() = default;
    [[nodiscard]] virtual size_t Index() const = 0;
};

}  // namespace Tree2Secondaries

#endif  // T2S_SECONDARY_PARTICLE_HXX
