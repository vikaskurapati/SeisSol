#include <Parallel/MPI.h>
#include <Solver/time_stepping/GhostTimeCluster.h>

#include "GhostTimeCluster.h"

#ifdef ACL_DEVICE
#include "device.h"
#endif // ACL_DEVICE


namespace seissol::time_stepping {
void GhostTimeCluster::sendCopyLayer(){
  SCOREP_USER_REGION( "sendCopyLayer", SCOREP_USER_REGION_TYPE_FUNCTION )
  assert(ct.correctionTime > lastSendTime);

#ifdef ACL_DEVICE
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
#endif

  lastSendTime = ct.correctionTime;
  for (unsigned int region = 0; region < meshStructure->numberOfRegions; ++region) {
    if (meshStructure->neighboringClusters[region][1] == static_cast<int>(otherGlobalClusterId)) {
      const auto messageSize = static_cast<int>(meshStructure->copyRegionSizes[region]);
#ifdef ACL_DEVICE
      device.api->copyBetween(meshStructure->deviceCopyRegions[region],
                              meshStructure->copyRegions[region],
                              messageSize * sizeof(real));
      real* sendBuffer = meshStructure->deviceCopyRegions[region];
#else
      real* sendBuffer = meshStructure->copyRegions[region];
#endif

      MPI_Isend(sendBuffer,
                messageSize,
                MPI_C_REAL,
                meshStructure->neighboringClusters[region][0],
                timeData + meshStructure->sendIdentifiers[region],
                seissol::MPI::mpi.comm(),
                meshStructure->sendRequests + region);
      sendQueue.push_back(region);
    }
  }
} void GhostTimeCluster::receiveGhostLayer(){
  SCOREP_USER_REGION( "receiveGhostLayer", SCOREP_USER_REGION_TYPE_FUNCTION )
  assert(ct.predictionTime > lastSendTime);
  for (unsigned int region = 0; region < meshStructure->numberOfRegions; ++region) {
    if (meshStructure->neighboringClusters[region][1] == static_cast<int>(otherGlobalClusterId) ) {
#ifdef ACL_DEVICE
      real* receiveBuffer = meshStructure->deviceGhostRegions[region];
#else
      real* receiveBuffer = meshStructure->ghostRegions[region];
#endif
      MPI_Irecv(receiveBuffer,
                static_cast<int>(meshStructure->ghostRegionSizes[region]),
                MPI_C_REAL,
                meshStructure->neighboringClusters[region][0],
                timeData + meshStructure->receiveIdentifiers[region],
                seissol::MPI::mpi.comm(),
                meshStructure->receiveRequests + region);
      receiveQueue.push_back(region);
    }
  }
}

void GhostTimeCluster::postprocessReceive(const MeshStructure* meshStructure, unsigned int region) {
#ifdef ACL_DEVICE
  device::DeviceInstance& device = device::DeviceInstance::getInstance();

  const auto messageSize = static_cast<int>(meshStructure->ghostRegionSizes[region]);
  device.api->copyBetween(meshStructure->ghostRegions[region],
                          meshStructure->deviceGhostRegions[region],
                          messageSize * sizeof(real));
#endif // ACL_DEVICE
}

bool GhostTimeCluster::testQueue(MPI_Request* requests,
                                 std::list<unsigned int>& regions,
                                 CallbackType postprocess) {
  for (auto region = regions.begin(); region != regions.end();) {
    MPI_Request *regionRequest = &requests[*region];
    int testSuccess = 0;
    MPI_Test(regionRequest, &testSuccess, MPI_STATUS_IGNORE);
    if (testSuccess) {
      if (postprocess) {
        postprocess(meshStructure, *region);
      }
      region = regions.erase(region);
    } else {
      ++region;
    }
  }
  return regions.empty();
}

bool GhostTimeCluster::testForGhostLayerReceives(){
  SCOREP_USER_REGION( "testForGhostLayerReceives", SCOREP_USER_REGION_TYPE_FUNCTION )
  return testQueue(meshStructure->receiveRequests,
                   receiveQueue,
                   GhostTimeCluster::postprocessReceive);
}


bool GhostTimeCluster::testForCopyLayerSends(){
  SCOREP_USER_REGION( "testForCopyLayerSends", SCOREP_USER_REGION_TYPE_FUNCTION )
  return testQueue(meshStructure->sendRequests, sendQueue);
}

ActResult GhostTimeCluster::act() {
  // Always check for receives/send for quicker MPI progression.
  testForGhostLayerReceives();
  testForCopyLayerSends();
  return AbstractTimeCluster::act();
}

void GhostTimeCluster::start() {
  assert(testForGhostLayerReceives());
  receiveGhostLayer();
}

void GhostTimeCluster::predict() {
    // Doesn't do anything
}

void GhostTimeCluster::correct() {
    // Doesn't do anything
}
bool GhostTimeCluster::mayCorrect() {
  return testForCopyLayerSends() && AbstractTimeCluster::mayCorrect();
}
bool GhostTimeCluster::mayPredict() {
  return testForGhostLayerReceives() && AbstractTimeCluster::mayPredict();
}

bool GhostTimeCluster::maySync() {
  return testForGhostLayerReceives() && testForCopyLayerSends() && AbstractTimeCluster::maySync();
}

void GhostTimeCluster::handleAdvancedPredictionTimeMessage(const NeighborCluster&) {
  assert(testForCopyLayerSends());
  sendCopyLayer();
}
void GhostTimeCluster::handleAdvancedCorrectionTimeMessage(const NeighborCluster& neighborCluster) {
  assert(testForGhostLayerReceives());

  auto upcomingCorrectionSteps = ct.stepsSinceLastSync;
  if (state == ActorState::Predicted) {
      upcomingCorrectionSteps = ct.nextCorrectionSteps();
  }

  const bool ignoreMessage = upcomingCorrectionSteps >= ct.stepsUntilSync;

  // If we are already at a sync point, we must not post an additional receive, as otherwise start() posts an additional
  // request!
  // This is also true for the last sync point (i.e. end of simulation), as in this case we do not want to have any
  // hanging request.
  if (!ignoreMessage) {
    receiveGhostLayer();
  }
}

GhostTimeCluster::GhostTimeCluster(double maxTimeStepSize,
                                   int timeStepRate,
                                   int globalTimeClusterId,
                                   int otherGlobalTimeClusterId,
                                   const MeshStructure *meshStructure)
    : AbstractTimeCluster(maxTimeStepSize, timeStepRate),
      globalClusterId(globalTimeClusterId),
      otherGlobalClusterId(otherGlobalTimeClusterId),
      meshStructure(meshStructure) {
}
void GhostTimeCluster::reset() {
  AbstractTimeCluster::reset();
  assert(testForGhostLayerReceives());
  lastSendTime = -1;
}

  void GhostTimeCluster::printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) {
    const auto rank = MPI::mpi.rank();
    logWarning(rank)
        << "Ghost: No update since " << timeSinceLastUpdate.count()
        << "[s] for global cluster " << globalClusterId
        << " with other cluster id " << otherGlobalClusterId
        << " at state " << actorStateToString(state)
        << " mayPredict = " << mayPredict()
        << " mayPredict (steps) = " << AbstractTimeCluster::mayPredict()
        << " mayCorrect = " << mayCorrect()
        << " mayCorrect (steps) = " << AbstractTimeCluster::mayCorrect()
        << " maySync = " << maySync();
    for (auto& neighbor : neighbors) {
      logWarning(rank)
        << "Neighbor with rate = " << neighbor.ct.timeStepRate
        << "PredTime = " << neighbor.ct.predictionTime
        << "CorrTime = " << neighbor.ct.correctionTime
        << "predictionsSinceSync = " << neighbor.ct.predictionsSinceLastSync
        << "correctionsSinceSync = " << neighbor.ct.stepsSinceLastSync;
    }
  }


}
