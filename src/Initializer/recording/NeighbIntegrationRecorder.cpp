#include "Recorder.h"
#include <Kernels/Interface.hpp>
#include <yateto.h>

using namespace device;
using namespace seissol::initializers;
using namespace seissol::initializers::recording;

void NeighbIntegrationRecorder::record(seissol::initializers::LTS &Handler, seissol::initializers::Layer &Layer) {
  kernels::NeighborData::Loader Loader;
  Loader.load(Handler, Layer);
  LayerContainer &Container = Layer.getLayerContainer();
  auto &Table = Container.getTableReferenceToInit();

  std::unordered_map<real *, real *> IDofsAddressRegistery{};
  size_t IDofsAddressCounter = 0;

  // get base pointers for the temp. memory
  real *IDofsScratchMem = static_cast<real *>(Layer.scratch_mem(Handler.idofs_scratch));

  real *(*FaceNeighbors)[4] = Layer.var(Handler.faceNeighbors);

  {
    if (Layer.getNumberOfCells()) {
      std::vector<real *> LtsIDofsPtrs{};
      std::vector<real *> LtsDerivativesPtrs{};

      std::vector<real *> GtsDerivativesPtrs{};
      std::vector<real *> GtsIDofsPtrs{};

      for (unsigned Cell = 0; Cell < Layer.getNumberOfCells(); ++Cell) {
        auto Data = Loader.entry(Cell);

        for (unsigned Face = 0; Face < 4; ++Face) {
          real *NeighbourBuffer = FaceNeighbors[Cell][Face];

          // check whether a neighbour element idofs has not been counted twice
          if ((IDofsAddressRegistery.find(NeighbourBuffer) == IDofsAddressRegistery.end())) {

            // maybe, because of BCs, a pointer can be a nullptr, i.e. skip it
            if (NeighbourBuffer != nullptr) {

              if (Data.cellInformation.faceTypes[Face] != outflow &&
                  Data.cellInformation.faceTypes[Face] != dynamicRupture) {

                bool IsNeighbProvidesDerivatives = ((Data.cellInformation.ltsSetup >> Face) % 2) == 1;

                if (IsNeighbProvidesDerivatives) {
                  real *NextTempIDofsPtr = &IDofsScratchMem[IDofsAddressCounter];

                  bool IsGtsNeigbour = ((Data.cellInformation.ltsSetup >> (Face + 4)) % 2) == 1;
                  if (IsGtsNeigbour) {

                    IDofsAddressRegistery[NeighbourBuffer] = NextTempIDofsPtr;
                    GtsIDofsPtrs.push_back(NextTempIDofsPtr);
                    GtsDerivativesPtrs.push_back(NeighbourBuffer);

                  } else {
                    IDofsAddressRegistery[NeighbourBuffer] = NextTempIDofsPtr;
                    LtsIDofsPtrs.push_back(NextTempIDofsPtr);
                    LtsDerivativesPtrs.push_back(NeighbourBuffer);
                  }
                  IDofsAddressCounter += tensor::I::size();
                } else {
                  IDofsAddressRegistery[NeighbourBuffer] = NeighbourBuffer;
                }
              }
            }
          }
        }
      }

      // step 1: time evaluation of dofs
      if (!GtsIDofsPtrs.empty()) {
        ConditionalKey Key(*KernelNames::neighbor_flux, *ComputationKind::with_gts_derivatives);
        checkKey(Table, Key);
        Table[Key].m_Container[*VariableID::derivatives] = new DevicePointers(GtsDerivativesPtrs);
        Table[Key].m_Container[*VariableID::idofs] = new DevicePointers(GtsIDofsPtrs);
        Table[Key].setNotEmpty();
      }

      if (!LtsIDofsPtrs.empty()) {
        ConditionalKey Key(*KernelNames::neighbor_flux, *ComputationKind::with_lts_derivatives);
        checkKey(Table, Key);
        Table[Key].m_Container[*VariableID::derivatives] = new DevicePointers(LtsDerivativesPtrs);
        Table[Key].m_Container[*VariableID::idofs] = new DevicePointers(LtsIDofsPtrs);
        Table[Key].setNotEmpty();
      }
    }
  }

  std::array<std::vector<real *>[*FaceRelations::Count], *FaceId::Count> RegularPeriodicDofs {};
  std::array<std::vector<real *>[*FaceRelations::Count], *FaceId::Count> RegularPeriodicIDofs {};
  std::array<std::vector<real *>[*FaceRelations::Count], *FaceId::Count> RegularPeriodicAminusT {};

  std::vector<real *> FreeSurfaceDofs[*FaceId::Count];
  std::vector<real *> FreeSurfaceIDofs[*FaceId::Count];
  std::vector<real *> FreeSurfaceAminusT[*FaceId::Count];

  std::array<std::vector<real *>[*DrFaceRelations::Count], *FaceId::Count> DrDofs {};
  std::array<std::vector<real *>[*DrFaceRelations::Count], *FaceId::Count> DrGodunov {};
  std::array<std::vector<real *>[*DrFaceRelations::Count], *FaceId::Count> DrFluxSolver {};

  CellDRMapping(*DrMapping)[4] = Layer.var(Handler.drMapping);

  // step 2: evaluation of neighbour flux integrals
  // auto neighbour_data = Loader.entry(data.cellInformation.faceNeighborIds[face]);
  for (unsigned Cell = 0; Cell < Layer.getNumberOfCells(); ++Cell) {
    auto Data = Loader.entry(Cell);
    for (unsigned int Face = 0; Face < 4; Face++) {

      real *NeighbourBufferPtr = FaceNeighbors[Cell][Face];
      // maybe, because of BCs, a pointer can be a nullptr, i.e. skip it
      if (NeighbourBufferPtr != nullptr) {

        switch (Data.cellInformation.faceTypes[Face]) {
        case regular:
          // Fallthrough intended
        case periodic: {
          // compute face type relation
          unsigned FaceRelation =
              Data.cellInformation.faceRelations[Face][1] + 3 * Data.cellInformation.faceRelations[Face][0] + 12 * Face;

          assert((*FaceRelations::Count) > FaceRelation &&
                 "ERROR::BINNING::NEIGB.::reg./period. incorrect face ralation count has been "
                 "detected");

          RegularPeriodicDofs[Face][FaceRelation].push_back(static_cast<real *>(Data.dofs));
          RegularPeriodicIDofs[Face][FaceRelation].push_back(IDofsAddressRegistery[NeighbourBufferPtr]);
          RegularPeriodicAminusT[Face][FaceRelation].push_back(
              static_cast<real *>(Data.neighboringIntegrationDevice.nAmNm1[Face]));
          break;
        }
        case freeSurface: {
          FreeSurfaceDofs[Face].push_back(static_cast<real *>(Data.dofs));
          FreeSurfaceIDofs[Face].push_back(IDofsAddressRegistery[NeighbourBufferPtr]);
          FreeSurfaceAminusT[Face].push_back(static_cast<real *>(Data.neighboringIntegrationDevice.nAmNm1[Face]));
          break;
        }
        case dynamicRupture: {
          unsigned face_relation = DrMapping[Cell][Face].side + 4 * DrMapping[Cell][Face].faceRelation;
          assert((*DrFaceRelations::Count) > face_relation &&
                 "ERROR::BINNING::NEIGB.::dyn. rupture incorrect face ralation count has been "
                 "detected");

          DrDofs[Face][face_relation].push_back(static_cast<real *>(Data.dofs));
          DrGodunov[Face][face_relation].push_back(DrMapping[Cell][Face].godunov);
          DrFluxSolver[Face][face_relation].push_back(DrMapping[Cell][Face].fluxSolver);

          break;
        }
        case outflow:
          break;
        default: {
          std::cout << "condition::" << Data.cellInformation.faceTypes[Face] << std::endl;
          throw std::string("ERROR::unknown boundary condition");
        }
        }
      }
    }
  }

  for (unsigned int Face = 0; Face < 4; Face++) {
    // regular and periodic
    for (unsigned FaceRelation = 0; FaceRelation < (*FaceRelations::Count); ++FaceRelation) {
      if (!RegularPeriodicDofs[Face][FaceRelation].empty()) {
        ConditionalKey Key(
            *KernelNames::neighbor_flux, (FaceKinds::regular || FaceKinds::periodic), Face, FaceRelation);
        checkKey(Table, Key);

        Table[Key].m_Container[*VariableID::idofs] = new DevicePointers(RegularPeriodicIDofs[Face][FaceRelation]);
        Table[Key].m_Container[*VariableID::dofs] = new DevicePointers(RegularPeriodicDofs[Face][FaceRelation]);
        Table[Key].m_Container[*VariableID::AminusT] = new DevicePointers(RegularPeriodicAminusT[Face][FaceRelation]);
      }
    }

    // free surface
    if (!FreeSurfaceDofs[Face].empty()) {
      ConditionalKey Key(*KernelNames::neighbor_flux, *FaceKinds::freeSurface, Face);
      checkKey(Table, Key);

      Table[Key].m_Container[*VariableID::idofs] = new DevicePointers(FreeSurfaceIDofs[Face]);
      Table[Key].m_Container[*VariableID::dofs] = new DevicePointers(FreeSurfaceDofs[Face]);
      Table[Key].m_Container[*VariableID::AminusT] = new DevicePointers(FreeSurfaceAminusT[Face]);
    }

    // dynamic rupture
    for (unsigned FaceRelation = 0; FaceRelation < (*DrFaceRelations::Count); ++FaceRelation) {
      if (!DrDofs[Face][FaceRelation].empty()) {
        ConditionalKey Key(*KernelNames::neighbor_flux, *FaceKinds::dynamicRupture, Face, FaceRelation);
        checkKey(Table, Key);

        Table[Key].m_Container[*VariableID::dofs] = new DevicePointers(DrDofs[Face][FaceRelation]);
        Table[Key].m_Container[*VariableID::godunov] = new DevicePointers(DrGodunov[Face][FaceRelation]);
        Table[Key].m_Container[*VariableID::fluxSolver] = new DevicePointers(DrFluxSolver[Face][FaceRelation]);
      }
    }
  }
}