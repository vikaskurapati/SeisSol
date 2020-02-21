#ifndef SEISSOL_RECORDER_H
#define SEISSOL_RECORDER_H

#include <vector>

#include "../tree/Layer.hpp"
#include "../LTS.h"

namespace seissol {
  namespace initializers {
    namespace recording {
      class AbstractRecorder {
      public:
        virtual ~AbstractRecorder() {}

        virtual void record(seissol::initializers::LTS &handler, seissol::initializers::Layer &layer) = 0;

      protected:
        void checkKey(conditional_table_t &table, ConditionalKey &key) {
          if (table.find(key) != table.end()) {
            throw std::string("ERROR::BINNING:: a table key conflict detected. Problems with hashing.");
          }
        }
      };

      class CompositeRecorder : public AbstractRecorder {
      public:
        ~CompositeRecorder() {
          for (auto Recorder: ConcreteRecorders)
            delete Recorder;
        }

        void record(seissol::initializers::LTS &handler, seissol::initializers::Layer &layer) override {
          for (auto Recorder: ConcreteRecorders) {
            Recorder->record(handler, layer);
          }
        }

        void addRecorder(AbstractRecorder *Recorder) { ConcreteRecorders.push_back(Recorder); }

        void removeRecorder(long unsigned RecorderIndex) {
          if (RecorderIndex < ConcreteRecorders.size()) {
            ConcreteRecorders.erase(ConcreteRecorders.begin() + RecorderIndex);
          }
        }

      private:
        std::vector<AbstractRecorder *> ConcreteRecorders{};
      };


      class LocalIntegrationRecorder : public AbstractRecorder {
      public:
        void record(seissol::initializers::LTS &handler, seissol::initializers::Layer &layer) override;
      };


      class NeighbIntegrationRecorder : public AbstractRecorder {
      public:
        void record(seissol::initializers::LTS &handler, seissol::initializers::Layer &layer) override;
      };

      class PlasticityRecorder : public AbstractRecorder {
      public:
        void record(seissol::initializers::LTS &handler, seissol::initializers::Layer &layer) override;
      };
    }
  }
}

#endif //SEISSOL_RECORDER_H
