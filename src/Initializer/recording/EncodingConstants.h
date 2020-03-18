#ifndef ENCODING_CONSTANTS_H_
#define ENCODING_CONSTANTS_H_

#include "specific_types.h"

namespace seissol {
  namespace initializers {
    namespace recording {

      enum struct VariableID : encode_t {
        dofs = 0,
        idofs,
        star,
        buffers,
        derivatives,
        AplusT,
        AminusT,
        godunov,
        fluxSolver,
        ivelocities, // 6th, 7the and 8th columns of idofs
        displacements,
        NodalStressTensor,
        Pstrains,
        elements_ids,
        Count
      };

      const encode_t ALL_BITS = ~static_cast<encode_t>(0);
      constexpr encode_t encode_any(unsigned Count) { return ~(ALL_BITS << Count); }

      enum struct KernelNames : encode_t {
        time = 1 << 0,
        volume = 1 << 1,
        local_flux = 1 << 2,
        neighbor_flux = 1 << 3,
        displacements = 1 << 4,
        plasticity = 1 << 5,
        Count = 6,
        any = encode_any(Count)
      };

      enum struct ComputationKind : encode_t {
        without_derivatives = 1 << 0,
        with_derivatives = 1 << 1,
        with_lts_derivatives = 1 << 2,
        with_gts_derivatives = 1 << 3,
        with_gts_buffers = 1 << 4,
        with_lts_buffers = 1 << 5,
        Count = 6,
        any = encode_any(Count)
      };

      enum struct FaceKinds : encode_t {
        regular = 1 << 0,
        freeSurface = 1 << 1,
        outflow = 1 << 2,
        dynamicRupture = 1 << 3,
        periodic = 1 << 4,
        Count = 5,
        any = encode_any(Count)
      };

      enum struct FaceId : encode_t { Count = 4, any = ALL_BITS };
      enum struct FaceRelations : encode_t { Count = 48, any = ALL_BITS };
      enum struct DrFaceRelations : encode_t { Count = 16, any = ALL_BITS };

      enum struct ExchangeInfo : encode_t {
        buffers = 1 << 0,
        derivatives = 1 << 1,
        Count = 2,
        any = encode_any(Count)
      };

    } // namespace recording
  }   // namespace initializers
} // namespace seissol
#endif // ENCODING_CONSTANTS_H_