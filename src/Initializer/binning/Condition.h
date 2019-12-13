#ifndef CONDITION_H_
#define CONDITION_H_

#include <assert.h>
#include <type_traits>


#include "EncodingConstants.h"
#include "specific_types.h"

using namespace seissol::initializers::binning;
using namespace seissol;


template <class T>
class Condition {
public:
    Condition() = delete;

    Condition(T i_value) : m_value(static_cast<encode_t>(i_value)) {
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

    Condition& operator!() {
        m_value = m_high_bits_mask & (~m_value);
        return *this;
    }

    Condition& operator||(const Condition &other) {
        m_value = m_value | other.m_value;
        return *this;
    }

    Condition& negate() {
       return !(*this);
    }

    encode_t get_encoding() {return m_value;}
    encode_t m_high_bits_mask;

private:
    Condition(encode_t value, encode_t count) : m_value(value), m_count(count) {}


    encode_t m_value;
    encode_t m_count;
};


/** Fast way to perform "OR" operation
 * NOTE: this method return 'encode_t' data type. You won't be able negate
 * the condition after that. Please, use Condition class to perform 
 * more sophisticated logical operations
 *
 * NOTE: The function was designed to handle simple condition encodings
 * Refer to Condition Class if you need much more sophisticated behaviour
*/
template<typename T>
typename std::enable_if<std::is_same<FaceKinds, T>::value ||
                        std::is_same<KernelNames, T>::value ||
                        std::is_same<ComputationKind, T>:: value ||
                        std::is_same<ExchangeInfo, T>::value, encode_t>::type
operator||(const T& lhs, const T& rhs) {
    return (static_cast<encode_t>(lhs) | static_cast<encode_t>(rhs));
}


template<typename T>
typename std::enable_if<std::is_same<FaceKinds, T>::value ||
                        std::is_same<KernelNames, T>::value ||
                        std::is_same<ComputationKind, T>:: value ||
                        std::is_same<ExchangeInfo, T>::value, encode_t>::type
operator!(const T& condition) {
    encode_t high_bits_mask = ~((~encode_t(0)) << static_cast<encode_t>(T::Count));
    return high_bits_mask & (~static_cast<encode_t>(condition));
}

/** Returns the actual value of enum item.
 * 
 * the behavior is similar as pointer dereference.
 * 
 * NOTE: The function was designed to handle simple condition encodings
 * Refer to Condition Class if you need much more sophisticated behaviour
*/
template<typename T>
constexpr
typename std::enable_if<std::is_same<FaceKinds, T>::value ||
                        std::is_same<KernelNames, T>::value ||
                        std::is_same<FaceId, T>::value ||
                        std::is_same<FaceRelations, T>::value ||
                        std::is_same<DrFaceRelations, T>::value ||
                        std::is_same<ComputationKind, T>:: value ||
                        std::is_same<VariableID, T>::value ||
                        std::is_same<ExchangeInfo, T>::value, encode_t>::type
operator*(const T& condition) {
    return static_cast<encode_t>(condition);
}

#endif  //CONDITION_H_
