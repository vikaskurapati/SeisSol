#ifndef SEISSOL_STF_H
#define SEISSOL_STF_H

#include <stddef.h>

#include "DynamicRupture/Misc.h"
#include "Initializer/DynamicRupture.h"
#include "Numerical_aux/GaussianNucleationFunction.h"
#include "Numerical_aux/RegularizedYoffe.h"
#include "Kernels/precision.hpp"

namespace seissol {
namespace initializers {
class Layer;
}
} // namespace seissol
namespace seissol {
namespace initializers {
struct DynamicRupture;
}
} // namespace seissol

namespace seissol::dr::friction_law {
class YoffeSTF {
  private:
  real (*onsetTime)[misc::numPaddedPoints];
  real (*tauS)[misc::numPaddedPoints];
  real (*tauR)[misc::numPaddedPoints];

  public:
  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime);

  real evaluate(real currentTime,
                [[maybe_unused]] real timeIncrement,
                size_t ltsFace,
                size_t pointIndex);
};

class GaussianSTF {
  private:
  real (*onsetTime)[misc::numPaddedPoints];
  real (*riseTime)[misc::numPaddedPoints];

  public:
  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime);

  real evaluate(real currentTime, real timeIncrement, size_t ltsFace, size_t pointIndex);
};

} // namespace seissol::dr::friction_law
#endif // SEISSOL_STF_H
