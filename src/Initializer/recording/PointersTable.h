#ifndef INDEX_TABLE_H_
#define INDEX_TABLE_H_

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "specific_types.h"
#include <device.h>

namespace seissol {
  namespace initializers {
    namespace recording {

      class DevicePointers {
      public:
        DevicePointers(std::vector<real *> Pointers) : m_Pointers(Pointers), m_DevicePtrs(nullptr) {
          if (!m_Pointers.empty()) {
            m_DevicePtrs = (real **)m_Device.api->allocGlobMem(m_Pointers.size() * sizeof(real *));
            m_Device.api->copyTo(m_DevicePtrs, m_Pointers.data(), m_Pointers.size() * sizeof(real *));
          }
        }

        explicit DevicePointers(const DevicePointers &Other) : m_Pointers(Other.m_Pointers), m_DevicePtrs(nullptr) {
          if (!m_Pointers.empty()) {
            if (Other.m_DevicePtrs != nullptr) {
              m_DevicePtrs = (real **)m_Device.api->allocGlobMem(Other.m_Pointers.size() * sizeof(real *));
              m_Device.api->copyBetween(m_DevicePtrs, Other.m_DevicePtrs, Other.m_Pointers.size() * sizeof(real *));
            }
          }
        }

        DevicePointers &operator=(const DevicePointers &Other) = delete;

        virtual ~DevicePointers() {
          if (m_DevicePtrs != nullptr) {
            m_Device.api->freeMem(m_DevicePtrs);
            m_DevicePtrs = nullptr;
          }
          m_Pointers.clear();
        }

        real **getPointers() { return m_DevicePtrs; }
        index_t getSize() { return m_Pointers.size(); }

      private:
        std::vector<real *> m_Pointers{};
        real **m_DevicePtrs{};
        device::DeviceInstance &m_Device = device::DeviceInstance::getInstance();
      };
    } // namespace recording
  }   // namespace initializers
} // namespace seissol

namespace seissol {
  namespace initializers {
    namespace recording {
      struct PointersTable {

      public:
        PointersTable() : m_IsEmpty(true) {
          for (auto &Pointers : m_Container)
            Pointers = nullptr;
        }

        std::array<DevicePointers *, *VariableID::Count> m_Container;

        bool isEmpty() { return m_IsEmpty; };
        void setNotEmpty() { m_IsEmpty = false; }

      private:
        bool m_IsEmpty; // TODO: doesn't make sense to keep this flag???
      };
    } // namespace recording
  }   // namespace initializers
} // namespace seissol
#endif