#ifndef CONDITIONAL_TABLE_H_
#define CONDITIONAL_TABLE_H_

#include <utility>
#include "specific_types.h"

namespace seissol {
    namespace initializers {
        namespace binning {

            struct ConditionalKey {
                ConditionalKey(encode_t i_kernel, 
                               encode_t i_type,
                               encode_t i_relation) : kernel(i_kernel),
                                                      type(i_type),
                                                      face_relation(i_relation) {};
                encode_t kernel;
                encode_t type;
                encode_t face_relation;

                bool operator==(const ConditionalKey &key) const {
                    return ((kernel == key.kernel) && 
                            (type == key.type) && 
                            (face_relation == key.face_relation));
                }
            };

            template <class T>
            inline void _hash_combine(std::size_t& seed, const T& v) {
                std::hash<T> hasher;
                seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
            }


            template <class T>
            class ConditionalHash;

            template<>
            struct ConditionalHash<ConditionalKey>
            {
                std::size_t operator()(ConditionalKey const& key) const 
                {
                    std::size_t result = 0;
                    _hash_combine(result, key.kernel);
                    _hash_combine(result, key.type);
                    _hash_combine(result, key.face_relation);
                    return result;
                }
            };

        }
    }
}

#endif  // CONDITIONAL_TABLE_H_