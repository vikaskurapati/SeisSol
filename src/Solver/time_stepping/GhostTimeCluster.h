#ifndef SEISSOL_GHOSTTIMECLUSTER_H
#define SEISSOL_GHOSTTIMECLUSTER_H

#include <list>
#include "Initializer/typedefs.hpp"
#include "AbstractTimeCluster.h"

namespace seissol::time_stepping {

class GhostTimeCluster : public AbstractTimeCluster {
 private:
  const int globalClusterId;
  const int otherGlobalClusterId;
  const MeshStructure* meshStructure;
  std::list<unsigned int> sendQueue;
  std::list<unsigned int> receiveQueue;

  double lastSendTime = -1.0;

  void sendCopyLayer();
  void receiveGhostLayer();

  using CallbackType = void (*)(const MeshStructure* meshStructure, unsigned int region);
  bool testQueue(MPI_Request* requests,
                 std::list<unsigned int>& queue,
                 CallbackType postprocess = nullptr);

  static void postprocessReceive(const MeshStructure* meshStructure, unsigned int region);

  bool testForCopyLayerSends();
  bool testForGhostLayerReceives();

  void start() override;
  void predict() override;
  void correct() override;
  bool mayPredict() override;
  bool mayCorrect() override;
  bool maySync() override;
  void handleAdvancedPredictionTimeMessage(const NeighborCluster& neighborCluster) override;
  void handleAdvancedCorrectionTimeMessage(const NeighborCluster& neighborCluster) override;
  void printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) override;

 public:
  GhostTimeCluster(double maxTimeStepSize,
                   int timeStepRate,
                   int globalTimeClusterId,
                   int otherGlobalTimeClusterId,
                   const MeshStructure* meshStructure
  );
  void reset() override;
  ActResult act() override;

};


}



#endif //SEISSOL_GHOSTTIMECLUSTER_H
