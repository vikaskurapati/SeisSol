#ifndef INDEX_TABLE_H_ 
#define INDEX_TABLE_H_

#include <unordered_map>
#include <vector>
#include <utility>
#include <string>

#include "specific_types.h"
#include <device_utils.h>

namespace seissol {
    namespace initializers {
        namespace binning {

            class BasicIndices {
            public:
                BasicIndices(std::vector<index_t> i_indices) : m_indices(i_indices), m_device_ptr(nullptr) {

                    if (!m_indices.empty()) {
                        m_device_ptr = (index_t *) device_malloc(m_indices.size() * sizeof(index_t));
                        device_copy_to(m_device_ptr, m_indices.data(), m_indices.size() * sizeof(index_t));
                    }
                }

                explicit
                BasicIndices(const BasicIndices& other) : m_indices(other.m_indices), m_device_ptr(nullptr) {
                    if (!m_indices.empty()) {
                        if (other.m_device_ptr != nullptr) {
                            m_device_ptr = (index_t *) device_malloc(other.m_indices.size() * sizeof(index_t));
                            device_copy_between(m_device_ptr, other.m_device_ptr, other.m_indices.size() * sizeof(index_t));
                        }
                    }
                }

                BasicIndices& operator=(const BasicIndices& other) = delete;

                virtual ~BasicIndices() {
                    if (m_device_ptr != nullptr) {
                        device_free(m_device_ptr);
                        m_device_ptr = nullptr;
                    }
                    m_indices.clear();
                }


            // private:
                std::vector<index_t> m_indices{};
                index_t *m_device_ptr{};

            };

            class RelativeIndices : public BasicIndices {
            public:

                explicit
                RelativeIndices(std::vector<index_t> i_indices,
                                index_t i_cell_id) : cell_id(i_cell_id), BasicIndices(i_indices) {}

                explicit
                RelativeIndices(const RelativeIndices& other) : BasicIndices(other), cell_id(other.cell_id) {}

                RelativeIndices& operator=(const RelativeIndices& other) = delete;

                virtual ~RelativeIndices() {}

                /** Defines an index of the element relatively which the indices were recorded.
                * */

            // private:
                index_t cell_id{};
            };


        }
    }
}

namespace seissol {
    namespace initializers {
        namespace binning {
            struct IndexTable {

                IndexTable() : m_is_empty(true) {
                    for (auto& container_ptr: variable_indices)
                        container_ptr = nullptr;
                }

            public:
                std::array<BasicIndices*, *VariableID::Count> variable_indices;

                bool is_empty() { return m_is_empty; };
                void set_not_empty_flag() { m_is_empty = false; }
            private:
                bool m_is_empty;  // TODO: doesn't make sense to keep this flag???
            };
        }
    }
}


#endif 