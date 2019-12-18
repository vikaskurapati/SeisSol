#ifdef ACL_DEVICE

#include <iostream>
#include <utility>
#include <vector>
#include <array>
#include <unordered_map>
#include <iterator>
#include <string>
#include <limits>

#include "algorithm.h"
#include "EncodingConstants.h"
#include "Condition.h"
#include "ConditionalTable.h"
#include "PointersTable.h"
#include "specific_types.h"

#include <Kernels/Interface.hpp>
#include <yateto.h>

using namespace seissol::initializers::binning;
using namespace seissol;

namespace seissol {
  namespace initializers {
    namespace binning {

      void check_key(conditional_table_t &table, ConditionalKey &key) {
        if (table.find(key) != table.end()) {
          throw std::string("ERROR::BINNING:: a table key conflict detected. Problems with hashing.");
        }
      }

      void local_integral(seissol::initializers::LTS &handler,
                          seissol::initializers::Layer &layer);

      void neighbour_integral(seissol::initializers::LTS &handler,
                              seissol::initializers::Layer &layer);
    }
  }
}

void seissol::initializers::binning::test(seissol::initializers::LTS &handler, seissol::initializers::Layer &layer) {
  local_integral(handler, layer);
  neighbour_integral(handler, layer);
}


void seissol::initializers::binning::local_integral(seissol::initializers::LTS &handler, seissol::initializers::Layer &layer) {

  kernels::LocalData::Loader loader;
  loader.load(handler, layer);
  LayerContainer &container = layer.getLayerContainer();
  auto &table = container.get_table_reference_to_init();

  // allocate counters and the registry
  std::unordered_map<index_t, real*> idofs_address_registry{};

  size_t idofs_address_counter = 0;
  size_t derivatives_address_counter = 0;

  // get base pointers for the temp. memory
  real* idofs_scratch_mem = static_cast<real*>(layer.scratch_mem(handler.idofs_scratch));
  real* derivatives_scratch_mem = static_cast<real*>(layer.scratch_mem(handler.derivatives_scratch));

  // *************************************** time and volume integrals ***************************************
  {
    if (layer.getNumberOfCells()) {

        std::vector<real*> dofs_ptrs{};
        std::vector<real*> star_ptrs{};
        std::vector<real*> idofs_ptrs{};
        std::vector<real*> dQ_ptrs{};

        std::vector<real*> lts_buffers{};
        std::vector<real*> idofs_for_lts_buffers{};

        real ** derivatives = layer.var(handler.derivatives);
        real ** buffers = layer.var(handler.buffers);

        for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
          auto data = loader.entry(cell);

          // dofs
          dofs_ptrs.push_back(static_cast<real *>(data.dofs));

          // idofs
          real *next_idof_ptr = &idofs_scratch_mem[idofs_address_counter];
          bool l_BuffersProvided = ((data.cellInformation.ltsSetup >> 8) % 2) == 1;
          bool l_LtsBuffers = ((data.cellInformation.ltsSetup >> 10) % 2) == 1;

          if (l_BuffersProvided) {
            if (l_LtsBuffers) {
              // lts buffers may require either accumulation or overriding (in case of reset command)
              idofs_ptrs.push_back(next_idof_ptr);

              idofs_for_lts_buffers.push_back(next_idof_ptr);
              lts_buffers.push_back(buffers[cell]);

              idofs_address_registry[cell] = next_idof_ptr;
              idofs_address_counter += tensor::I::size();
            }
            else {
              // gts buffers have to be always overridden
              idofs_ptrs.push_back(buffers[cell]);
              idofs_address_registry[cell] = buffers[cell];
            }
          }
          else {
            idofs_ptrs.push_back(next_idof_ptr);
            idofs_address_registry[cell] = next_idof_ptr;
            idofs_address_counter += tensor::I::size();
          }

          // stars
          star_ptrs.push_back(static_cast<real *>(data.localIntegrationDevice.starMatrices[0]));

          // derivatives
          bool l_DerivativesProvided = ((data.cellInformation.ltsSetup >> 9) % 2) == 1;
          if (l_DerivativesProvided) {
            dQ_ptrs.push_back(derivatives[cell]);

          } else {
            dQ_ptrs.push_back(&derivatives_scratch_mem[derivatives_address_counter]);
            derivatives_address_counter += yateto::computeFamilySize<tensor::dQ>();
          }
        }

      ConditionalKey key(KernelNames::time || KernelNames::volume);
      check_key(table, key);

      table[key].container[*VariableID::dofs] = new DevicePointers(dofs_ptrs);
      table[key].container[*VariableID::star] = new DevicePointers(star_ptrs);
      table[key].container[*VariableID::idofs] = new DevicePointers(idofs_ptrs);
      table[key].container[*VariableID::derivatives] = new DevicePointers(dQ_ptrs);

      table[key].set_not_empty_flag();

      if (!idofs_for_lts_buffers.empty()) {
        ConditionalKey key(*KernelNames::time, *ComputationKind::with_lts_buffers);

        table[key].container[*VariableID::buffers] = new DevicePointers(lts_buffers);
        table[key].container[*VariableID::idofs] = new DevicePointers(idofs_for_lts_buffers);
        table[key].set_not_empty_flag();
      }
    }
  }

  // *************************************** local flux integral ***************************************
  {
    for (unsigned face = 0; face < 4; ++face) {

      std::vector<real*> idofs_ptrs{};
      std::vector<real*> dofs_ptrs{};
      std::vector<real*> AplusT_ptrs{};

      for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {

        auto data = loader.entry(cell);

        // no element local contribution in the case of dynamic rupture boundary conditions
        if(data.cellInformation.faceTypes[face] != dynamicRupture) {
          idofs_ptrs.push_back(idofs_address_registry[cell]);
          dofs_ptrs.push_back(static_cast<real*>(data.dofs));
          AplusT_ptrs.push_back(static_cast<real*>(data.localIntegrationDevice.nApNm1[face]));
        }
      }

      // NOTE: we can check any container, but we must check that a set is not empty!
      if (!dofs_ptrs.empty()) {
        ConditionalKey key(*KernelNames::local_flux, !FaceKinds::dynamicRupture, face);
        check_key(table, key);
        table[key].container[*VariableID::idofs] = new DevicePointers(idofs_ptrs);
        table[key].container[*VariableID::dofs] = new DevicePointers(dofs_ptrs);
        table[key].container[*VariableID::AplusT] = new DevicePointers(AplusT_ptrs);

        table[key].set_not_empty_flag();
      }
    }
  }
}


void seissol::initializers::binning::neighbour_integral(seissol::initializers::LTS &handler,
                                                        seissol::initializers::Layer &layer) {

  kernels::NeighborData::Loader loader;
  loader.load(handler, layer);
  LayerContainer &container = layer.getLayerContainer();
  auto &table = container.get_table_reference_to_init();

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
        check_key(table, key);
        table[key].container[*VariableID::derivatives] = new DevicePointers(gts_derivatives_ptrs);
        table[key].container[*VariableID::idofs] = new DevicePointers(gts_idofs_ptrs);
        table[key].set_not_empty_flag();
      }

      if (!lts_idofs_ptrs.empty()) {
        ConditionalKey key(*KernelNames::neighbor_flux,
                           *ComputationKind::with_lts_derivatives);
        check_key(table, key);
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
        check_key(table, key);

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
      check_key(table, key);

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
        check_key(table, key);

        table[key].container[*VariableID::dofs] = new DevicePointers(dr_dofs[face][face_relation]);
        table[key].container[*VariableID::godunov] = new DevicePointers(dr_godunov[face][face_relation]);
        table[key].container[*VariableID::fluxSolver] = new DevicePointers(dr_flux_solver[face][face_relation]);
      }
    }
  }
}
#endif