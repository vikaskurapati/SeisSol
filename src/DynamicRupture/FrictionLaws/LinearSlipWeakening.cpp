#include "LinearSlipWeakening.h"

#include <math.h>

#include "DynamicRupture/Parameters.h"
#include "init.h"
#include "kernel.h"

namespace seissol::dr::friction_law {

void NoSpecialization::resampleSlipRate(
    real (&resampledSlipRate)[dr::misc::numPaddedPoints],
    real const (&slipRateMagnitude)[dr::misc::numPaddedPoints]) {
  dynamicRupture::kernel::resampleParameter resampleKrnl;
  resampleKrnl.resample = init::resample::Values;
  resampleKrnl.originalQ = slipRateMagnitude;
  resampleKrnl.resampledQ = resampledSlipRate;
  resampleKrnl.execute();
}
void BiMaterialFault::copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                                         seissol::initializers::DynamicRupture* dynRup,
                                         real fullUpdateTime) {
  // maybe change later to const_cast?
  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTS_LinearSlipWeakeningBimaterial*>(dynRup);
  regularisedStrength = layerData.var(concreteLts->regularisedStrength);
}

real BiMaterialFault::strengthHook(real faultStrength,
                                   real localSlipRate,
                                   real deltaT,
                                   unsigned int ltsFace,
                                   unsigned int pointIndex) {
  // modify strength according to Prakash-Clifton
  // see e.g.: Pelties - Verification of an ADER-DG method for complex dynamic rupture problems
  real expterm = std::exp(-(std::max(static_cast<real>(0.0), localSlipRate) + drParameters.vStar) *
                          deltaT / drParameters.prakashLength);
  real newStrength =
      regularisedStrength[ltsFace][pointIndex] * expterm + faultStrength * (1.0 - expterm);
  regularisedStrength[ltsFace][pointIndex] = newStrength;
  return newStrength;
}

} // namespace seissol::dr::friction_law
