#ifndef CONDITIONAL_TABLE_H
#define CONDITIONAL_TABLE_H


#include "Initializer/recording/EncodingConstants.h"
#include "Initializer/recording/Condition.h"
#include "Initializer/recording/ConditionalKey.h"
#include "Initializer/recording/PointersTable.h"

using conditional_table_t = std::unordered_map<ConditionalKey, PointersTable, ConditionalHash<ConditionalKey>>;

#endif //CONDITIONAL_TABLE_H
