#include "STF.h"

#include "Initializer/DynamicRupture.h"
#include "Initializer/tree/Layer.hpp"
#include "Numerical_aux/GaussianNucleationFunction.h"
#include "Numerical_aux/RegularizedYoffe.h"

namespace seissol::dr::friction_law {
void YoffeSTF::copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                                  seissol::initializers::DynamicRupture* dynRup,
                                  real fullUpdateTime) {
  auto* concreteLts = dynamic_cast<seissol::initializers::LTS_ImposedSlipRatesYoffe*>(dynRup);
  onsetTime = layerData.var(concreteLts->onsetTime);
  tauS = layerData.var(concreteLts->tauS);
  tauR = layerData.var(concreteLts->tauR);
}

real YoffeSTF::evaluate(real currentTime,
                        [[maybe_unused]] real timeIncrement,
                        size_t ltsFace,
                        size_t pointIndex) {
  return regularizedYoffe::regularizedYoffe(currentTime - onsetTime[ltsFace][pointIndex],
                                            tauS[ltsFace][pointIndex],
                                            tauR[ltsFace][pointIndex]);
}

void GaussianSTF::copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                                     seissol::initializers::DynamicRupture* dynRup,
                                     real fullUpdateTime) {
  auto* concreteLts = dynamic_cast<seissol::initializers::LTS_ImposedSlipRatesGaussian*>(dynRup);
  onsetTime = layerData.var(concreteLts->onsetTime);
  riseTime = layerData.var(concreteLts->riseTime);
}

real GaussianSTF::evaluate(real currentTime,
                           real timeIncrement,
                           size_t ltsFace,
                           size_t pointIndex) {
  real smoothStepIncrement = gaussianNucleationFunction::smoothStepIncrement(
      currentTime - onsetTime[ltsFace][pointIndex], timeIncrement, riseTime[ltsFace][pointIndex]);
  return smoothStepIncrement / timeIncrement;
}
} // namespace seissol::dr::friction_law
