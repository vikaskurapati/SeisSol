#ifndef SEISSOL_RECORDER_H
#define SEISSOL_RECORDER_H

#include "ConditionalTable.h"

#include <vector>
#include <Initializer/LTS.h>
#include <Initializer/tree/Layer.hpp>

namespace seissol {
  namespace initializers {
    namespace recording {
      class AbstractRecorder {
      public:
        virtual ~AbstractRecorder() {}

        virtual void record(seissol::initializers::LTS &Handler, seissol::initializers::Layer &Layer) = 0;

      protected:
        void checkKey(conditional_table_t &Table, ConditionalKey &Key) {
          if (Table.find(Key) != Table.end()) {
            throw std::string("ERROR::BINNING:: a table key conflict detected. Problems with hashing.");
          }
        }
      };

      class CompositeRecorder : public AbstractRecorder {
      public:
        ~CompositeRecorder() {
          for (auto Recorder : m_ConcreteRecorders)
            delete Recorder;
        }

        void record(seissol::initializers::LTS &Handler, seissol::initializers::Layer &Layer) override {
          for (auto Recorder : m_ConcreteRecorders) {
            Recorder->record(Handler, Layer);
          }
        }

        void addRecorder(AbstractRecorder *Recorder) { m_ConcreteRecorders.push_back(Recorder); }

        void removeRecorder(long unsigned RecorderIndex) {
          if (RecorderIndex < m_ConcreteRecorders.size()) {
            m_ConcreteRecorders.erase(m_ConcreteRecorders.begin() + RecorderIndex);
          }
        }

      private:
        std::vector<AbstractRecorder *> m_ConcreteRecorders{};
      };

      class LocalIntegrationRecorder : public AbstractRecorder {
      public:
        void record(seissol::initializers::LTS &Handler, seissol::initializers::Layer &Layer) override;
      };

      class NeighbIntegrationRecorder : public AbstractRecorder {
      public:
        void record(seissol::initializers::LTS &Handler, seissol::initializers::Layer &Layer) override;
      };

      class PlasticityRecorder : public AbstractRecorder {
      public:
        void record(seissol::initializers::LTS &Handler, seissol::initializers::Layer &Layer) override;
      };
    } // namespace recording
  }   // namespace initializers
} // namespace seissol

#endif // SEISSOL_RECORDER_H
