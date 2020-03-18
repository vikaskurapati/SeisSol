#ifndef LAYER_CONTAINER_H_
#define LAYER_CONTAINER_H_

#include <assert.h>
#include <utility>
#include <vector>
#include <unordered_map>

#include <Kernels/precision.hpp>
#include "Initializer/recording/ConditionalTable.h"

namespace seissol {
    namespace initializers {
        class LayerContainer {

        public:
            void setConditionalTable(conditional_table_t Table) {
              m_CondTable = Table;
            }

            conditional_table_t& getConditionalTable() {
                //assert(!m_CondTable.empty() && "conditional table hasn't been initialized");
                return m_CondTable;
            }

            void freeConditionalTable() {
                for (auto& IndexTable: m_CondTable) {
                    for (auto Pointers: IndexTable.second.m_Container) {
                        if (Pointers != nullptr) {
                            delete Pointers;
                        }
                    }
                }
            }

            conditional_table_t& getTableReferenceToInit() {
                return m_CondTable;
            }

        private:
            conditional_table_t m_CondTable{};
        };
    }
}
#endif  //LAYER_CONTAINER_H_
