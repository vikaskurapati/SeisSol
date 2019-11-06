#ifndef INDEX_TABLE_H_ 
#define INDEX_TABLE_H_

#include <unordered_map>
#include <vector>
#include <utility>

#include "specific_types.h"

namespace seissol {
    namespace initializers {
        namespace binning {
            
            struct IndexTable{
                std::vector<std::pair<index_t, index_t>> element_face_pairs;
                std::unordered_map<encode_t, std::vector<index_t>> variable_indices;
            };
            
        }
    }
}

#endif 