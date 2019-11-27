#ifndef DEVICE_GEN_CODE_DEVICE_KERNEL_H_
#define DEVICE_GEN_CODE_DEVICE_KERNEL_H_
#include <cmath>
#include <limits>
#include "tensor.h"
#include <stdlib.h>
using namespace seissol;
namespace device_gen_code {
  namespace kernel {
    struct volume {
      constexpr static unsigned long const NonZeroFlops = 34839;
      constexpr static unsigned long const HardwareFlops = 123336;

      double const* I{};
      double* Q{};
      tensor::kDivM::Container<double const*> kDivM;
      tensor::star::Container<double const*> star;

      unsigned* I_indices{};
      unsigned* Q_indices{};
      tensor::kDivM::Container<unsigned*> kDivM_indices;
      tensor::star::Container<unsigned*> star_indices;

      size_t num_elements{};

      void execute();
    };
  }
  namespace kernel {
    struct localFlux {
      constexpr static unsigned long const NonZeroFlops[] = {9936, 10080, 31968, 27216};
      constexpr static unsigned long const HardwareFlops[] = {49248, 49248, 49248, 49248};

      double const* AplusT{};
      double const* I{};
      double* Q{};
      tensor::fMrT::Container<double const*> fMrT;
      tensor::rDivM::Container<double const*> rDivM;

      unsigned* AplusT_indices{};
      unsigned* I_indices{};
      unsigned* Q_indices{};
      tensor::fMrT::Container<unsigned*> fMrT_indices;
      tensor::rDivM::Container<unsigned*> rDivM_indices;

      size_t num_elements{};

      struct Prefetch {
        double const* I{};
        double const* Q{};
      };
      Prefetch _prefetch;

      void execute0();
      void execute1();
      void execute2();
      void execute3();
      typedef void (localFlux::* const member_function_ptr)(void);
      constexpr static member_function_ptr ExecutePtrs[] = {&localFlux::execute0, &localFlux::execute1, &localFlux::execute2, &localFlux::execute3};
      constexpr static member_function_ptr findExecute(unsigned i0) {
        return ExecutePtrs[1*i0];
      }
      inline void execute(unsigned i0) {
        (this->*findExecute(i0))();
      }
      constexpr static unsigned long nonZeroFlops(unsigned i0) {
        return NonZeroFlops[1*i0];
      }
      constexpr static unsigned long hardwareFlops(unsigned i0) {
        return HardwareFlops[1*i0];
      }
    };
  }
  namespace kernel {
    struct derivativeTaylorExpansion {
      constexpr static unsigned long const NonZeroFlops[] = {504, 630, 360, 180, 72, 18};
      constexpr static unsigned long const HardwareFlops[] = {0, 0, 0, 0, 0, 0};

      double power = std::numeric_limits<double>::signaling_NaN();
      double* I{};
      tensor::dQ::Container<double const*> dQ;

      unsigned* I_indices{};
      tensor::dQ::Container<unsigned*> dQ_indices;

      size_t num_elements{};

      void execute0();
      void execute1();
      void execute2();
      void execute3();
      void execute4();
      void execute5();
      typedef void (derivativeTaylorExpansion::* const member_function_ptr)(void);
      constexpr static member_function_ptr ExecutePtrs[] = {&derivativeTaylorExpansion::execute0, &derivativeTaylorExpansion::execute1, &derivativeTaylorExpansion::execute2, &derivativeTaylorExpansion::execute3, &derivativeTaylorExpansion::execute4, &derivativeTaylorExpansion::execute5};
      constexpr static member_function_ptr findExecute(unsigned i0) {
        return ExecutePtrs[1*i0];
      }
      inline void execute(unsigned i0) {
        (this->*findExecute(i0))();
      }
      constexpr static unsigned long nonZeroFlops(unsigned i0) {
        return NonZeroFlops[1*i0];
      }
      constexpr static unsigned long hardwareFlops(unsigned i0) {
        return HardwareFlops[1*i0];
      }
    };
  }
  namespace kernel {
    struct derivative {
      constexpr static unsigned long const NonZeroFlops[] = {0, 34524, 13806, 4716, 1260, 216};
      constexpr static unsigned long const HardwareFlops[] = {0, 122472, 45360, 17496, 3672, 2376};

      tensor::dQ::Container<double*> dQ;
      tensor::kDivMT::Container<double const*> kDivMT;
      tensor::star::Container<double const*> star;

      tensor::dQ::Container<unsigned*> dQ_indices;
      tensor::kDivMT::Container<unsigned*> kDivMT_indices;
      tensor::star::Container<unsigned*> star_indices;

      size_t num_elements{};

      void execute1();
      void execute2();
      void execute3();
      void execute4();
      void execute5();
      typedef void (derivative::* const member_function_ptr)(void);
      constexpr static member_function_ptr ExecutePtrs[] = {nullptr, &derivative::execute1, &derivative::execute2, &derivative::execute3, &derivative::execute4, &derivative::execute5};
      constexpr static member_function_ptr findExecute(unsigned i0) {
        return ExecutePtrs[1*i0];
      }
      inline void execute(unsigned i0) {
        (this->*findExecute(i0))();
      }
      constexpr static unsigned long nonZeroFlops(unsigned i0) {
        return NonZeroFlops[1*i0];
      }
      constexpr static unsigned long hardwareFlops(unsigned i0) {
        return HardwareFlops[1*i0];
      }
    };
  }
}
#endif
