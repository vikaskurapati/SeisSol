#ifndef ENCODING_CONSTANTS_H_
#define ENCODING_CONSTANTS_H_

#include "specific_types.h"


namespace seissol {
    namespace initializers {
        namespace binning {

            enum struct VariableID: encode_t{dofs = 0,
                                             idofs,
                                             start,
                                             buffers,
                                             derivatives,
                                             AplusT,
                                             AminusT,
                                             elements_ids,
                                             Count};


            const encode_t ALL_BITS = ~static_cast<encode_t>(0);
            constexpr encode_t encode_any(unsigned Count) {
                return ~(ALL_BITS << Count);
            }

            enum struct KernelNames: encode_t{time = 1 << 0,
                                              volume = 1 << 1,
                                              local_flux = 1 << 2,
                                              neighbor_flux = 1 << 3,
                                              Count = 4,
                                              any = encode_any(Count)};


            enum struct TimeComputationKind: encode_t{with_derivatives = 1 << 0,
                                                      without_derivatives = 1 << 1,
                                                      with_gts_buffers = 1 << 2,
                                                      with_lts_buffers = 1 << 3,
                                                      Count = 4,
                                                      any = encode_any(Count)};


            enum struct FaceKinds: encode_t{regular = 1 << 0, 
                                            freeSurface = 1 << 1, 
                                            dynamicRupture = 1 << 2,
                                            outflow = 1 << 3,
                                            periodic = 1 << 4,
                                            Count = 5,
                                            any = encode_any(Count)};


            enum struct FaceRelations: encode_t{Count = 48,
                                                any = ALL_BITS};


            enum struct ExchangeInfo: encode_t{buffers = 1 << 0,
                                               derivatives = 1 << 1,
                                               Count = 2,
                                               any = encode_any(Count)};

        }
    }
}
#endif  //ENCODING_CONSTANTS_H_