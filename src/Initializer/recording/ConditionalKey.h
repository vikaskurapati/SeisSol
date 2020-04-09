#ifndef CONDITIONAL_KEY_H_
#define CONDITIONAL_KEY_H_

#include "specific_types.h"
#include <limits>
#include <utility>
#include <functional>

namespace seissol {
  namespace initializers {
    namespace recording {

      struct ConditionalKey {
        ConditionalKey(encode_t Kernel, encode_t Type = std::numeric_limits<encode_t>::max(),
                       encode_t FaceId = std::numeric_limits<encode_t>::max(),
                       encode_t Relation = std::numeric_limits<encode_t>::max())
            : m_Kernel(Kernel), m_Type(Type), m_FaceId(FaceId), m_FaceRelation(Relation){};

        bool operator==(const ConditionalKey &Key) const {
          return ((m_Kernel == Key.m_Kernel) && (m_Type == Key.m_Type) && (m_FaceId == Key.m_FaceId) &&
                  (m_FaceRelation == Key.m_FaceRelation));
        }

        encode_t m_Kernel;
        encode_t m_Type;
        encode_t m_FaceId;
        encode_t m_FaceRelation;
      };

      template <class T> inline void hashCombine(std::size_t &Seed, const T &Value) {
        std::hash<T> Hasher;
        Seed ^= Hasher(Value) + 0x9e3779b9 + (Seed << 6) + (Seed >> 2);
      }

      template <class T> class ConditionalHash;

      template <> struct ConditionalHash<ConditionalKey> {
        std::size_t operator()(ConditionalKey const &Key) const {
          std::size_t Result = 0;
          hashCombine(Result, Key.m_Kernel);
          hashCombine(Result, Key.m_Type);
          hashCombine(Result, Key.m_FaceId);
          hashCombine(Result, Key.m_FaceRelation);
          return Result;
        }
      };
    } // namespace recording
  }   // namespace initializers
} // namespace seissol
#endif // CONDITIONAL_KEY_H_
