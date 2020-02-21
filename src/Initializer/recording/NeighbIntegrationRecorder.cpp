#include "Recorder.h"
#include <Kernels/Interface.hpp>
#include <yateto.h>

using namespace device;
using namespace seissol::initializers;
using namespace seissol::initializers::recording;

void NeighbIntegrationRecorder::record(seissol::initializers::LTS &handler, seissol::initializers::Layer &layer) {
  kernels::NeighborData::Loader loader;
  loader.load(handler, layer);
  LayerContainer &container = layer.getLayerContainer();
  auto &table = container.getTableReferenceToInit();

  std::unordered_map<real*, real*> idofs_address_registery{};
  size_t idofs_address_counter = 0;

  // get base pointers for the temp. memory
  real* idofs_scratch_mem = static_cast<real*>(layer.scratch_mem(handler.idofs_scratch));

  real* (*faceNeighbors)[4] = layer.var(handler.faceNeighbors);

  {
    if (layer.getNumberOfCells()) {
      std::vector<real *> lts_idofs_ptrs{};
      std::vector<real *> lts_derivatives_ptrs{};

      std::vector<real *> gts_derivatives_ptrs{};
      std::vector<real *> gts_idofs_ptrs{};

      for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
        auto data = loader.entry(cell);

        for (unsigned face = 0; face < 4; ++face) {
          real *neighbour_buffer = faceNeighbors[cell][face];

          // check whether a neighbour element idofs has not been counted twice
          if ((idofs_address_registery.find(neighbour_buffer) == idofs_address_registery.end())) {

            // maybe, because of BCs, a pointer can be a nullptr, i.e. skip it
            if (neighbour_buffer != nullptr) {

              if(data.cellInformation.faceTypes[face] != outflow && data.cellInformation.faceTypes[face] != dynamicRupture) {

                bool l_NeighbProvidesDerivatives = ((data.cellInformation.ltsSetup >> face) % 2) == 1;

                if (l_NeighbProvidesDerivatives) {
                  real *next_temp_idofs_ptr = &idofs_scratch_mem[idofs_address_counter];

                  bool l_isGtsNeigbour = ((data.cellInformation.ltsSetup >> (face + 4)) % 2) == 1;
                  if (l_isGtsNeigbour) {

                    idofs_address_registery[neighbour_buffer] = next_temp_idofs_ptr;
                    gts_idofs_ptrs.push_back(next_temp_idofs_ptr);
                    gts_derivatives_ptrs.push_back(neighbour_buffer);

                  } else {
                    idofs_address_registery[neighbour_buffer] = next_temp_idofs_ptr;
                    lts_idofs_ptrs.push_back(next_temp_idofs_ptr);
                    lts_derivatives_ptrs.push_back(neighbour_buffer);
                  }
                  idofs_address_counter += tensor::I::size();
                } else {
                  idofs_address_registery[neighbour_buffer] = neighbour_buffer;
                }
              }
            }
          }
        }
      }

      // step 1: time evaluation of dofs
      if (!gts_idofs_ptrs.empty()) {
        ConditionalKey key(*KernelNames::neighbor_flux,
                           *ComputationKind::with_gts_derivatives);
        checkKey(table, key);
        table[key].container[*VariableID::derivatives] = new DevicePointers(gts_derivatives_ptrs);
        table[key].container[*VariableID::idofs] = new DevicePointers(gts_idofs_ptrs);
        table[key].set_not_empty_flag();
      }

      if (!lts_idofs_ptrs.empty()) {
        ConditionalKey key(*KernelNames::neighbor_flux,
                           *ComputationKind::with_lts_derivatives);
        checkKey(table, key);
        table[key].container[*VariableID::derivatives] = new DevicePointers(lts_derivatives_ptrs);
        table[key].container[*VariableID::idofs] = new DevicePointers(lts_idofs_ptrs);
        table[key].set_not_empty_flag();
      }
    }
  }

  std::array<std::vector<real*>[*FaceRelations::Count], *FaceId::Count> regular_periodic_dofs{};
  std::array<std::vector<real*>[*FaceRelations::Count], *FaceId::Count> regular_periodic_idofs{};
  std::array<std::vector<real*>[*FaceRelations::Count], *FaceId::Count> regular_periodic_AminusT{};


  std::vector<real*> freeSurface_dofs[*FaceId::Count];
  std::vector<real*> freeSurface_idofs[*FaceId::Count];
  std::vector<real*> freeSurface_AminusT[*FaceId::Count];

  std::array<std::vector<real*>[*DrFaceRelations::Count], *FaceId::Count> dr_dofs{};
  std::array<std::vector<real*>[*DrFaceRelations::Count], *FaceId::Count> dr_godunov{};
  std::array<std::vector<real*>[*DrFaceRelations::Count], *FaceId::Count> dr_flux_solver{};

  CellDRMapping (*drMapping)[4] = layer.var(handler.drMapping);

  // step 2: evaluation of neighbour flux integrals
  //auto neighbour_data = loader.entry(data.cellInformation.faceNeighborIds[face]);
  for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
    auto data = loader.entry(cell);
    for(unsigned int face = 0; face < 4; face++) {

      real *neighbour_buffer_ptr = faceNeighbors[cell][face];
      // maybe, because of BCs, a pointer can be a nullptr, i.e. skip it
      if (neighbour_buffer_ptr != nullptr) {

        switch (data.cellInformation.faceTypes[face]) {
          case regular:
            // Fallthrough intended
          case periodic: {
            // compute face type relation
            unsigned face_relation = data.cellInformation.faceRelations[face][1]
                                     + 3 * data.cellInformation.faceRelations[face][0]
                                     + 12 * face;

            assert((*FaceRelations::Count) > face_relation &&
                   "ERROR::BINNING::NEIGB.::reg./period. incorrect face ralation count has been detected");

            regular_periodic_dofs[face][face_relation].push_back(static_cast<real*>(data.dofs));
            regular_periodic_idofs[face][face_relation].push_back(idofs_address_registery[neighbour_buffer_ptr]);
            regular_periodic_AminusT[face][face_relation].push_back(static_cast<real*>(data.neighboringIntegrationDevice.nAmNm1[face]));
            break;
          }
          case freeSurface: {
            freeSurface_dofs[face].push_back(static_cast<real*>(data.dofs));
            freeSurface_idofs[face].push_back(idofs_address_registery[neighbour_buffer_ptr]);
            freeSurface_AminusT[face].push_back(static_cast<real*>(data.neighboringIntegrationDevice.nAmNm1[face]));
            break;
          }
          case dynamicRupture: {
            unsigned face_relation = drMapping[cell][face].side + 4 * drMapping[cell][face].faceRelation;
            assert((*DrFaceRelations::Count) > face_relation &&
                   "ERROR::BINNING::NEIGB.::dyn. rupture incorrect face ralation count has been detected");

            dr_dofs[face][face_relation].push_back(static_cast<real*>(data.dofs));
            dr_godunov[face][face_relation].push_back(drMapping[cell][face].godunov);
            dr_flux_solver[face][face_relation].push_back(drMapping[cell][face].fluxSolver);

            break;
          }
          case outflow:
            break;
          default: {
            std::cout << "condition::" << data.cellInformation.faceTypes[face] << std::endl;
            throw std::string("ERROR::unknown boundary condition");
          }
        }
      }
    }
  }

  for(unsigned int face = 0; face < 4; face++) {
    // regular and periodic
    for (unsigned face_relation = 0; face_relation < (*FaceRelations::Count); ++face_relation) {
      if (!regular_periodic_dofs[face][face_relation].empty()) {
        ConditionalKey key(*KernelNames::neighbor_flux,
                           (FaceKinds::regular || FaceKinds::periodic),
                           face,
                           face_relation);
        checkKey(table, key);

        table[key].container[*VariableID::idofs] = new DevicePointers(regular_periodic_idofs[face][face_relation]);
        table[key].container[*VariableID::dofs] = new DevicePointers(regular_periodic_dofs[face][face_relation]);
        table[key].container[*VariableID::AminusT] = new DevicePointers(regular_periodic_AminusT[face][face_relation]);
      }
    }

    // free surface
    if (!freeSurface_dofs[face].empty()) {
      ConditionalKey key(*KernelNames::neighbor_flux,
                         *FaceKinds::freeSurface,
                         face);
      checkKey(table, key);

      table[key].container[*VariableID::idofs] = new DevicePointers(freeSurface_idofs[face]);
      table[key].container[*VariableID::dofs] = new DevicePointers(freeSurface_dofs[face]);
      table[key].container[*VariableID::AminusT] = new DevicePointers(freeSurface_AminusT[face]);
    }

    // dynamic rupture
    for (unsigned face_relation = 0; face_relation < (*DrFaceRelations::Count); ++face_relation) {
      if (!dr_dofs[face][face_relation].empty()) {
        ConditionalKey key(*KernelNames::neighbor_flux,
                           *FaceKinds::dynamicRupture,
                           face,
                           face_relation);
        checkKey(table, key);

        table[key].container[*VariableID::dofs] = new DevicePointers(dr_dofs[face][face_relation]);
        table[key].container[*VariableID::godunov] = new DevicePointers(dr_godunov[face][face_relation]);
        table[key].container[*VariableID::fluxSolver] = new DevicePointers(dr_flux_solver[face][face_relation]);
      }
    }
  }
}