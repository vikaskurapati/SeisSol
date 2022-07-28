#ifndef SEISSOL_LINEARSLIPWEAKENINGINITIALIZER_H
#define SEISSOL_LINEARSLIPWEAKENINGINITIALIZER_H

#include <string>
#include <unordered_map>

#include "BaseDRInitializer.h"
#include "Initializer/tree/LTSInternalNode.hpp"
#include "Kernels/precision.hpp"

namespace seissol {
namespace initializers {
class LTSTree;
}
} // namespace seissol
namespace seissol {
namespace initializers {
struct DynamicRupture;
}
} // namespace seissol

namespace seissol::dr::initializers {

/**
 * Derived initializer class for the LinearSlipWeakening friction law
 */
class LinearSlipWeakeningInitializer : public BaseDRInitializer {
  public:
  using BaseDRInitializer::BaseDRInitializer;
  /**
   * Computes initial friction and slip rates
   */
  virtual void initializeFault(seissol::initializers::DynamicRupture* dynRup,
                               seissol::initializers::LTSTree* dynRupTree) override;

  protected:
  /**
   * Adds the additional parameters mu_s, mu_d, d_c, cohesion and if available forced_rupture_time.
   */
  virtual void
      addAdditionalParameters(std::unordered_map<std::string, real*>& parameterToStorageMap,
                              seissol::initializers::DynamicRupture* dynRup,
                              seissol::initializers::LTSInternalNode::leaf_iterator& it) override;
};

/**
 * Derived initializer class for the LinearSlipWeakening friction law with bimaterial regularization
 */
class LinearSlipWeakeningBimaterialInitializer : public LinearSlipWeakeningInitializer {
  public:
  using LinearSlipWeakeningInitializer::LinearSlipWeakeningInitializer;
  /**
   * Computes initial value for the regularized strength
   */
  virtual void initializeFault(seissol::initializers::DynamicRupture* dynRup,
                               seissol::initializers::LTSTree* dynRupTree) override;
};

} // namespace seissol::dr::initializers
#endif // SEISSOL_LINEARSLIPWEAKENINGINITIALIZER_H
