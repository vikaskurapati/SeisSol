#ifndef CONDITION_H_
#define CONDITION_H_

#include <assert.h>
#include <type_traits>

#include "EncodingConstants.h"
#include "specific_types.h"

using namespace seissol::initializers::recording;
using namespace seissol;

template <class T> class Condition {
public:
  Condition() = delete;

  Condition(T i_value) : m_Value(static_cast<encode_t>(i_value)) {
    /*
    static_assert(std::is_same<encode_t,
                               typename std::underlying_type<KernelNames>::type>::value,
                  "enum base type must be int");
    */

    static_assert(static_cast<encode_t>(T::Count) > 1,
                  "enum must have a positive value of Count"
                  "and greater than 1");

    m_high_bits_mask = ~((~encode_t(0)) << static_cast<encode_t>(T::Count));
  }

  Condition &operator!() {
    m_Value = m_high_bits_mask & (~m_Value);
    return *this;
  }

  Condition &operator||(const Condition &other) {
    m_Value = m_Value | other.m_Value;
    return *this;
  }

  Condition &negate() { return !(*this); }

  encode_t get_encoding() { return m_Value; }
  encode_t m_high_bits_mask;

private:
  Condition(encode_t Value, encode_t Count) : m_Value(Value), m_Count(Count) {}

  encode_t m_Value;
  encode_t m_Count;
};

/** Fast way to perform "OR" operation
 * NOTE: this method return 'encode_t' data type. You won't be able negate
 * the condition after that. Please, use Condition class to perform
 * more sophisticated logical operations
 *
 * NOTE: The function was designed to handle simple condition encodings
 * Refer to Condition Class if you need much more sophisticated behaviour
 */
template <typename T>
typename std::enable_if<std::is_same<FaceKinds, T>::value || std::is_same<KernelNames, T>::value ||
                            std::is_same<ComputationKind, T>::value || std::is_same<ExchangeInfo, T>::value,
                        encode_t>::type
operator||(const T &Lhs, const T &Rhs) {
  return (static_cast<encode_t>(Lhs) | static_cast<encode_t>(Rhs));
}

template <typename T>
typename std::enable_if<std::is_same<FaceKinds, T>::value || std::is_same<KernelNames, T>::value ||
                            std::is_same<ComputationKind, T>::value || std::is_same<ExchangeInfo, T>::value,
                        encode_t>::type
operator!(const T &Condition) {
  encode_t HighBitsMask = ~((~encode_t(0)) << static_cast<encode_t>(T::Count));
  return HighBitsMask & (~static_cast<encode_t>(Condition));
}

/** Returns the actual value of enum item.
 *
 * the behavior is similar as pointer dereference.
 *
 * NOTE: The function was designed to handle simple condition encodings
 * Refer to Condition Class if you need much more sophisticated behaviour
 */
template <typename T>
constexpr
    typename std::enable_if<std::is_same<FaceKinds, T>::value || std::is_same<KernelNames, T>::value ||
                                std::is_same<FaceId, T>::value || std::is_same<FaceRelations, T>::value ||
                                std::is_same<DrFaceRelations, T>::value || std::is_same<ComputationKind, T>::value ||
                                std::is_same<VariableID, T>::value || std::is_same<ExchangeInfo, T>::value,
                            encode_t>::type
    operator*(const T &Condition) {
  return static_cast<encode_t>(Condition);
}

#endif // CONDITION_H_
