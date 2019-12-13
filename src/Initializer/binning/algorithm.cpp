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
#include "IndexTable.h"
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

  // allocate counters for the main analysis
  size_t idofs_address_counter = 0;        // NOTE: reset at the very end
  std::unordered_map<index_t, index_t> idfos_address_registry{};
  size_t derivatives_address_counter = 0;  // NOTE: reset at the very end
  int distance = 0;

  // *************************************** time integral ***************************************

  // --------------------------------------- primary analysis::start ---------------------------------------
  // perform a primary analysis and find cells which have their own derivatives,
  // and cells which require to provide scratch memory for their derivatives computations
  std::vector<index_t> cells_without_derivatives{};
  std::vector<index_t> cells_with_derivatives{};

  for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
    auto data = loader.entry(cell);
    bool l_DerivativesProvided = ((data.cellInformation.ltsSetup >> 9) % 2) == 1;

    if (l_DerivativesProvided)
      cells_with_derivatives.push_back(cell);
    else
      cells_without_derivatives.push_back(cell);
  }

  // compute the "check-sum"
  index_t bin_check_sum = cells_with_derivatives.size() + cells_without_derivatives.size();
  assert(bin_check_sum == layer.getNumberOfCells() && "ERROR::BINNING::TIME_INTEGRATION::a bin check sum is different.");
  // if we pass the check above, it means that we prepare all cells for puting into appropriate bins
  // --------------------------------------- primary analysis::end ---------------------------------------

  // bin cells without their own derivatives buffers
  if (!cells_without_derivatives.empty()) {
    std::vector<index_t> dofs_indices{};
    std::vector<index_t> start_indices{};
    std::vector<index_t> idofs_indices{};  // relative to the tmp buffer, namely: idof part
    std::vector<index_t> dQ_indices{}; //relative to the tmp buffer, namely: derivative part

    // debugging
    std::vector<index_t> elements_ids{};

    index_t base_cell_id = cells_without_derivatives[0];
    auto base_data = loader.entry(base_cell_id);

    for (auto cell: cells_without_derivatives) {
      auto current_data = loader.entry(cell);

      // dofs
      distance = std::distance(static_cast<real *>(base_data.dofs), static_cast<real *>(current_data.dofs));
      assert((distance >= 0) && "ERROR::BINNING::TIME_INTEGRATION:: negative distance between elements detected");
      dofs_indices.push_back(distance);

      // idofs
      idofs_indices.push_back(idofs_address_counter);
      idfos_address_registry[cell] = idofs_address_counter;
      idofs_address_counter += tensor::I::size();

      // stars
      distance = std::distance(static_cast<real *>(base_data.localIntegration.starMatrices[0]),
                               static_cast<real *>(current_data.localIntegration.starMatrices[0]));

      assert((distance >= 0) && "ERROR::BINNING::TIME_INTEGRATION:: negative distance between elements detected");
      start_indices.push_back(distance);

      // dQ; NOTE: the first derivatives is dofs itself; recording the address of the second derivative
      dQ_indices.push_back(derivatives_address_counter);

      // NOTE: the size is bigger than needed, namely: we take into account the first derivative
      derivatives_address_counter += yateto::computeFamilySize<tensor::dQ>();

      // debugging
      elements_ids.push_back(cell);
    }


    // compute a key for this particular condition
    ConditionalKey key(*KernelNames::time, *ComputationKind::without_derivatives);
    check_key(table, key);

    // record indices to the table
    table[key].variable_indices[*VariableID::dofs] = new RelativeIndices(dofs_indices, base_cell_id);
    table[key].variable_indices[*VariableID::start] = new RelativeIndices(start_indices, base_cell_id);
    table[key].variable_indices[*VariableID::idofs] = new BasicIndices(idofs_indices);
    table[key].variable_indices[*VariableID::derivatives] = new BasicIndices(dQ_indices);

    //debugging
    table[key].variable_indices[*VariableID::elements_ids] = new BasicIndices(elements_ids);

    table[key].set_not_empty_flag();
  }


  //------------------------------- cells with derivatives ------------------------------------
  // bin cells with their own derivatives buffers
  if (!cells_with_derivatives.empty()) {

    std::vector<index_t> dofs_indices{};
    std::vector<index_t> start_indices{};
    std::vector<index_t> idofs_indices{};  // relative to the tmp buffer, namely: idof part
    std::vector<index_t> dQ_indices{};

    // debugging
    std::vector<index_t> elements_ids{};

    index_t base_cell_id = cells_with_derivatives[0];
    auto base_data = loader.entry(base_cell_id);

    real ** derivatives = layer.var(handler.derivatives);
    real *base_derivative_address = layer.var(handler.derivatives)[base_cell_id];

    for (auto cell: cells_with_derivatives) {
      auto current_data = loader.entry(cell);

      // dofs
      distance = std::distance(static_cast<real *>(base_data.dofs), static_cast<real *>(current_data.dofs));
      assert((distance >= 0) && "ERROR::BINNING::TIME_INTEGRATION:: negative distance between elements detected");
      dofs_indices.push_back(distance);

      // idofs
      idofs_indices.push_back(idofs_address_counter);
      idfos_address_registry[cell] = idofs_address_counter;
      idofs_address_counter += tensor::I::size();

      // stars
      distance = std::distance(static_cast<real *>(base_data.localIntegration.starMatrices[0]),
                               static_cast<real *>(current_data.localIntegration.starMatrices[0]));

      assert((distance >= 0) && "ERROR::BINNING::TIME_INTEGRATION:: negative distance between elements detected");
      start_indices.push_back(distance);

      // dQ
      distance = std::distance(base_derivative_address, derivatives[cell]);
      assert((distance >= 0) && "ERROR::BINNING::TIME_INTEGRATION:: negative distance between elements detected");
      dQ_indices.push_back(distance);

      // debugging
      elements_ids.push_back(cell);
    }

    ConditionalKey key(*KernelNames::time, *ComputationKind::with_derivatives);
    check_key(table, key);

    table[key].variable_indices[*VariableID::dofs] = new RelativeIndices(dofs_indices, base_cell_id);
    table[key].variable_indices[*VariableID::start] = new RelativeIndices(start_indices, base_cell_id);
    table[key].variable_indices[*VariableID::idofs] = new BasicIndices(idofs_indices);
    table[key].variable_indices[*VariableID::derivatives] = new RelativeIndices(dQ_indices, base_cell_id);

    //debugging
    table[key].variable_indices[*VariableID::elements_ids] = new BasicIndices(elements_ids);

    table[key].set_not_empty_flag();
  }


  // *************************************** volume integral ***************************************
  // bin cells for the volume integral
  // NOTE: some clusters can have no elements inside
  if (layer.getNumberOfCells()) {
    std::vector<index_t> idofs_indices{};
    std::vector<index_t> dofs_indices{};
    std::vector<index_t> start_indices{};


    // debugging
    std::vector<index_t> elements_ids{};

    index_t base_cell_id = 0;
    auto base_data = loader.entry(base_cell_id);

    for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
      auto current_data = loader.entry(cell);

      // idofs
      idofs_indices.push_back(idfos_address_registry[cell]);

      // dofs
      distance = std::distance(static_cast<real *>(base_data.dofs),
                               static_cast<real *>(current_data.dofs));
      assert((distance >= 0) && "ERROR::BINNING::VOLUME_INTEGRATION:: negative distance between elements detected");
      dofs_indices.push_back(distance);

      // stars
      distance = std::distance(static_cast<real *>(base_data.localIntegration.starMatrices[0]),
                               static_cast<real *>(current_data.localIntegration.starMatrices[0]));

      assert((distance >= 0) && "ERROR::BINNING::VOLUME_INTEGRATION:: negative distance between elements detected");
      start_indices.push_back(distance);

      // debugging
      elements_ids.push_back(cell);
    }

    ConditionalKey key(*KernelNames::volume);
    check_key(table, key);
    table[key].variable_indices[*VariableID::dofs] = new RelativeIndices(dofs_indices, base_cell_id);
    table[key].variable_indices[*VariableID::start] = new RelativeIndices(start_indices, base_cell_id);
    table[key].variable_indices[*VariableID::idofs] = new BasicIndices(idofs_indices);

    //debugging
    table[key].variable_indices[*VariableID::elements_ids] = new BasicIndices(elements_ids);

    table[key].set_not_empty_flag();
  }

  // *************************************** local flux integral ***************************************

  {
    for (unsigned face = 0; face < 4; ++face) {

      std::vector<index_t> idofs_indices{};
      std::vector<index_t> dofs_indices{};
      std::vector<index_t> AplusT_indices{};

      // debugging
      std::vector<index_t> elements_ids{};

      index_t base_cell_id = 0;
      auto base_data = loader.entry(base_cell_id);

      for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {

        auto current_data = loader.entry(cell);

        // no element local contribution in the case of dynamic rupture boundary conditions
        if(current_data.cellInformation.faceTypes[face] != dynamicRupture) {

          // idofs
          idofs_indices.push_back(idfos_address_registry[cell]);

          // dofs
          distance = std::distance(static_cast<real *>(base_data.dofs),
                                   static_cast<real *>(current_data.dofs));
          assert((distance >= 0) && "ERROR::BINNING::LOCAL_FLUX:: negative distance between elements detected");
          dofs_indices.push_back(distance);

          // AplusT
          distance = std::distance(static_cast<real *>(base_data.localIntegration.nApNm1[face]),
                                   static_cast<real *>(current_data.localIntegration.nApNm1[face]));
          assert((distance >= 0) && "ERROR::BINNING::LOCAL_FLUX:: negative distance between elements detected");
          AplusT_indices.push_back(distance);

          elements_ids.push_back(cell);
        }
      }

      // NOTE: we can check any container, but we must check that a set is not empty!
      if (!dofs_indices.empty()) {
        ConditionalKey key(*KernelNames::local_flux, !FaceKinds::dynamicRupture, face);
        check_key(table, key);
        table[key].variable_indices[*VariableID::idofs] = new BasicIndices(idofs_indices);
        table[key].variable_indices[*VariableID::dofs] = new RelativeIndices(dofs_indices, base_cell_id);
        table[key].variable_indices[*VariableID::AplusT] = new RelativeIndices(AplusT_indices, base_cell_id);

        //debugging
        table[key].variable_indices[*VariableID::elements_ids] = new BasicIndices(elements_ids);

        table[key].set_not_empty_flag();
      }
    }
  }

  // *************************************** end of the local integral ***************************************
  // --------------------------------------- primary analysis::start ---------------------------------------

  // perform a primary analysis and find which cells have gts_buffers,
  // and which cells have lts buffers
  std::vector<index_t> cells_with_lts_buffers{};
  std::vector<index_t> cells_with_gts_buffers{};

  for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
    auto data = loader.entry(cell);
    bool l_BuffersProvided = ((data.cellInformation.ltsSetup >> 8) % 2) == 1;
    bool l_LtsBuffers = ((data.cellInformation.ltsSetup >> 10) % 2) == 1;

    if (l_BuffersProvided) {
      if (l_LtsBuffers)
        cells_with_lts_buffers.push_back(cell);
      else
        cells_with_gts_buffers.push_back(cell);
    }
  }
  // --------------------------------------- primary analysis::end ---------------------------------------

  //  bin cells with lts buffers: must be either overridden or accumulated, namely: run-time behaviour
  if (!cells_with_lts_buffers.empty()) {

    std::vector<index_t> buffers_indices{};
    std::vector<index_t> idofs_indices{};  // relative to the tmp buffer, namely: idof part

    index_t base_cell_id = cells_with_lts_buffers[0];
    real ** buffers = layer.var(handler.buffers);
    real *base_buffer_address = layer.var(handler.buffers)[base_cell_id];

    for (auto cell: cells_with_lts_buffers) {
      auto current_data = loader.entry(cell);

      distance = std::distance(base_buffer_address, buffers[cell]);
      assert((distance >= 0) && "ERROR::BINNING::TIME_INTEGRATION:: negative distance between elements detected");
      buffers_indices.push_back(distance);

      idofs_indices.push_back(idfos_address_registry[cell]);
    }

    ConditionalKey key(*KernelNames::time, *ComputationKind::with_lts_buffers);
    check_key(table, key);

    table[key].variable_indices[*VariableID::buffers] = new RelativeIndices(buffers_indices, base_cell_id);
    table[key].variable_indices[*VariableID::idofs] = new BasicIndices(idofs_indices);
    table[key].set_not_empty_flag();
  }

  // bin cells with gts buffers: must be only overridden: static behaviour, i.e. depends on a mesh
  if (!cells_with_gts_buffers.empty()) {

    std::vector<index_t> buffers_indices{};
    std::vector<index_t> idofs_indices{};  // relative to the tmp buffer, namely: idof part

    index_t base_cell_id = cells_with_gts_buffers[0];
    real ** buffers = layer.var(handler.buffers);
    real *base_buffer_address = layer.var(handler.buffers)[base_cell_id];

    for (auto cell: cells_with_gts_buffers) {
      auto current_data = loader.entry(cell);

      distance = std::distance(base_buffer_address, buffers[cell]);
      assert((distance >= 0) && "ERROR::BINNING::TIME_INTEGRATION:: negative distance between elements detected");
      buffers_indices.push_back(distance);

      idofs_indices.push_back(idfos_address_registry[cell]);
    }

    ConditionalKey key(*KernelNames::time, *ComputationKind::with_gts_buffers);
    check_key(table, key);

    table[key].variable_indices[*VariableID::buffers] = new RelativeIndices(buffers_indices, base_cell_id);
    table[key].variable_indices[*VariableID::idofs] = new BasicIndices(idofs_indices);
    table[key].set_not_empty_flag();
  }
}


void seissol::initializers::binning::neighbour_integral(seissol::initializers::LTS &handler,
                                                        seissol::initializers::Layer &layer) {

  kernels::NeighborData::Loader loader;
  loader.load(handler, layer);
  LayerContainer &container = layer.getLayerContainer();
  auto &table = container.get_table_reference_to_init();

  // allocate counters for the main analysis
  size_t idofs_address_counter = 0;        // NOTE: reset at the very end
  int distance = 0;


  std::pair<std::vector<index_t>, std::vector<index_t>> cells_with_gts_derivatives{};
  std::pair<std::vector<index_t>, std::vector<index_t>> cells_with_lts_derivatives{};
  std::pair<std::vector<index_t>, std::vector<index_t>> cells_with_buffers{};
  std::unordered_map<real*, index_t> idofs_address_registery{};

  real* (*faceNeighbors)[4] = layer.var(handler.faceNeighbors);

  index_t base_cell_id = 0;
  real* base_buffer = layer.var(handler.buffers)[base_cell_id];

  // this step is supposed to uniquely arrange time integrated dofs
  // find cells which are either GTS-derivatives, or LTS-derivatives, or LTS
  for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
    auto data = loader.entry(cell);
    for (unsigned face = 0; face < 4; ++face) {

      // retrive ptr. to derivatives/buffer of a neighbour cell
      real* neighbour_buffer = faceNeighbors[cell][face];

      // check whether a neighbour element idofs has not been counted twice
      if ((idofs_address_registery.find(neighbour_buffer) == idofs_address_registery.end())) {

        // check whether we a neighbour element is valid
        if (neighbour_buffer != nullptr) {

          // record idofs counter to the registery
          idofs_address_registery[neighbour_buffer] = idofs_address_counter;

          distance = std::distance(base_buffer, neighbour_buffer);
          bool l_DerivativesProvided = ((data.cellInformation.ltsSetup >> face) % 2) == 1;

          if (l_DerivativesProvided) {
            bool l_isGtsNeigbour = ((data.cellInformation.ltsSetup >> (face + 4)) % 2) == 1;
            if (l_isGtsNeigbour) {
              // these cells will integrate their derivatives according to the GTS scheme
              // and put the idofs to the corresponding place
              cells_with_gts_derivatives.first.push_back(distance);
              cells_with_gts_derivatives.second.push_back(idofs_address_counter);
            } else {
              // these cells will integrate their derivatives according to the LTS scheme
              // and put the idofs to the corresponding place
              cells_with_lts_derivatives.first.push_back(distance);
              cells_with_lts_derivatives.second.push_back(idofs_address_counter);
            }
          }
          else {
            // these cells will copy the content of their buffers to idofs array
            // b/c their buffers contains integrated dofs
            cells_with_buffers.first.push_back(distance);
            cells_with_buffers.second.push_back(idofs_address_counter);
          }

          // update the address for the next cell
          idofs_address_counter += tensor::I::size();
        }
        else {
          idofs_address_registery[neighbour_buffer] = std::numeric_limits<index_t>::max();
        }
      }
    }
  }

  // step 1: time evaluation of dofs
  {
    if (!cells_with_gts_derivatives.first.empty()) {
      ConditionalKey key(*KernelNames::neighbor_flux,
                         *ComputationKind::with_gts_derivatives);
      check_key(table, key);
      table[key].variable_indices[*VariableID::buffers] = new RelativeIndices(cells_with_gts_derivatives.first, base_cell_id);
      table[key].variable_indices[*VariableID::idofs] = new BasicIndices(cells_with_gts_derivatives.second);
      table[key].set_not_empty_flag();
    }

    if (!cells_with_lts_derivatives.first.empty()) {
      ConditionalKey key(*KernelNames::neighbor_flux,
                         *ComputationKind::with_lts_derivatives);
      check_key(table, key);
      table[key].variable_indices[*VariableID::buffers] = new RelativeIndices(cells_with_lts_derivatives.first, base_cell_id);
      table[key].variable_indices[*VariableID::idofs] = new BasicIndices(cells_with_lts_derivatives.second);
      table[key].set_not_empty_flag();
    }

    if (!cells_with_buffers.first.empty()) {

      ConditionalKey key(*KernelNames::neighbor_flux,
                         (ComputationKind::with_gts_buffers || ComputationKind::with_lts_buffers));
      check_key(table, key);
      table[key].variable_indices[*VariableID::buffers] = new RelativeIndices(cells_with_buffers.first, base_cell_id);
      table[key].variable_indices[*VariableID::idofs] = new BasicIndices(cells_with_buffers.second);
      table[key].set_not_empty_flag();
    }
  }

  // step 2: computations of the neighbour integral
  std::array<std::vector<index_t>[*FaceRelations::Count], *FaceId::Count> regular_periodic_bins{};
  std::vector<index_t> freeSurface_bin[*FaceId::Count];
  std::array<std::vector<index_t>[*DrFaceRelations::Count], *FaceId::Count> dynamic_rupture_bins{};
  CellDRMapping (*drMapping)[4] = layer.var(handler.drMapping);

  //auto neighbour_data = loader.entry(data.cellInformation.faceNeighborIds[face]);
  for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
    auto data = loader.entry(cell);
    for(unsigned int face = 0; face < 4; face++) {

      switch (data.cellInformation.faceTypes[face]) {
        case regular:
          // Fallthrough intended
        case periodic:
        {
          // compute face type relation
          unsigned face_relation = data.cellInformation.faceRelations[face][1]
                                   + 3 * data.cellInformation.faceRelations[face][0]
                                   + 12 * face;

          assert((*FaceRelations::Count) > face_relation && "ERROR::BINNING::NEIGB.::reg./period. incorrect face ralation count has been detected");
          regular_periodic_bins[face][face_relation].push_back(cell);
          break;
        }
        case freeSurface:
        {
          freeSurface_bin[face].push_back(cell);
          break;
        }
        case dynamicRupture:
        {
          unsigned face_relation = drMapping[cell][face].side + 4 * drMapping[cell][face].faceRelation;
          assert((*DrFaceRelations::Count) > face_relation && "ERROR::BINNING::NEIGB.::dyn. rupture incorrect face ralation count has been detected");
          dynamic_rupture_bins[face][face_relation].push_back(cell);
          break;
        }
        case outflow: break;
        default:
        {
          std::cout << "condition::" << data.cellInformation.faceTypes[face] << std::endl;
          throw std::string("ERROR::unknown boundary condition");
        }
      }
    }
  }


  // recodring entries to the table for neighbour flux integration
  for (unsigned face = 0; face < 4; ++face) {

    // REGULAR and PERIODIC cells
    {
      for (unsigned face_relation = 0; face_relation < (*FaceRelations::Count); ++face_relation) {
        if (!regular_periodic_bins[face][face_relation].empty()) {

          std::vector<index_t> idofs_indices{};
          std::vector<index_t> dofs_indices{};
          std::vector<index_t> AminusT_indices{};

          index_t base_cell_id = regular_periodic_bins[face][face_relation][0];
          auto base_data = loader.entry(base_cell_id);
          for (auto cell: regular_periodic_bins[face][face_relation]) {
            real *neighbour_buffer_ptr = faceNeighbors[cell][face];

            if (neighbour_buffer_ptr != nullptr) {
              auto current_data = loader.entry(cell);
              // idofs
              idofs_indices.push_back(idofs_address_registery[neighbour_buffer_ptr]);

              // dofs
              distance = std::distance(static_cast<real *>(base_data.dofs),
                                       static_cast<real *>(current_data.dofs));
              assert((distance >= 0) && "ERROR::BINNING::NEIGHBOUR_FLUX:: negative distance between elements detected");
              dofs_indices.push_back(distance);

              // AminusT
              distance = std::distance(static_cast<real *>(base_data.neighboringIntegration.nAmNm1[face]),
                                       static_cast<real *>(current_data.neighboringIntegration.nAmNm1[face]));
              assert((distance >= 0) && "ERROR::BINNING::NEIGHBOUR_FLUX:: negative distance between elements detected");
              AminusT_indices.push_back(distance);
            }
          }

          if (!dofs_indices.empty()) {
            //std::cout << "reg::" << (FaceKinds::regular || FaceKinds::periodic) << "::" << face << "::" << face_relation << std::endl;

            ConditionalKey key(*KernelNames::neighbor_flux,
                               (FaceKinds::regular || FaceKinds::periodic),
                               face,
                               face_relation);
            check_key(table, key);

            table[key].variable_indices[*VariableID::idofs] = new BasicIndices(idofs_indices);
            table[key].variable_indices[*VariableID::dofs] = new RelativeIndices(dofs_indices, base_cell_id);
            table[key].variable_indices[*VariableID::AminusT] = new RelativeIndices(AminusT_indices, base_cell_id);
          }
        }
      }
    }

    // FREE SURFACE
    {
      if (!freeSurface_bin[face].empty()) {

        std::vector<index_t> idofs_indices{};
        std::vector<index_t> dofs_indices{};
        std::vector<index_t> AminusT_indices{};

        index_t base_cell_id = freeSurface_bin[face][0];
        auto base_data = loader.entry(base_cell_id);
        for (auto cell: freeSurface_bin[face]) {
          real *neighbour_buffer_ptr = faceNeighbors[cell][face];

          if (neighbour_buffer_ptr != nullptr) {
            auto current_data = loader.entry(cell);
            // idofs
            idofs_indices.push_back(idofs_address_registery[neighbour_buffer_ptr]);

            // dofs
            distance = std::distance(static_cast<real *>(base_data.dofs),
                                     static_cast<real *>(current_data.dofs));
            assert((distance >= 0) && "ERROR::BINNING::NEIGHBOUR_FLUX:: negative distance between elements detected");
            dofs_indices.push_back(distance);

            // AminusT
            distance = std::distance(static_cast<real *>(base_data.neighboringIntegration.nAmNm1[face]),
                                     static_cast<real *>(current_data.neighboringIntegration.nAmNm1[face]));
            assert((distance >= 0) && "ERROR::BINNING::NEIGHBOUR_FLUX:: negative distance between elements detected");
            AminusT_indices.push_back(distance);
          }
        }

        if (!dofs_indices.empty()) {
          ConditionalKey key(*KernelNames::neighbor_flux,
                             *FaceKinds::freeSurface,
                             face);
          check_key(table, key);

          table[key].variable_indices[*VariableID::idofs] = new BasicIndices(idofs_indices);
          table[key].variable_indices[*VariableID::dofs] = new RelativeIndices(dofs_indices, base_cell_id);
          table[key].variable_indices[*VariableID::AminusT] = new RelativeIndices(AminusT_indices, base_cell_id);
        }
      }
    }

    // Dynamic Rupture
    {
      for (unsigned face_relation = 0; face_relation < (*DrFaceRelations::Count); ++face_relation) {
        if (!dynamic_rupture_bins[face][face_relation].empty()) {

          std::vector<index_t> flux_solver_indices{};
          std::vector<index_t> godunov_indices{};
          std::vector<index_t> dofs_indices{};

          index_t base_cell_id = dynamic_rupture_bins[face][face_relation][0];
          auto base_data = loader.entry(base_cell_id);
          for (auto cell: dynamic_rupture_bins[face][face_relation]) {
            auto current_data = loader.entry(cell);

            // dofs
            distance = std::distance(static_cast<real *>(base_data.dofs),
                                     static_cast<real *>(current_data.dofs));
            assert((distance >= 0) && "ERROR::BINNING::NEIGHBOUR_FLUX:: negative distance between elements detected");
            dofs_indices.push_back(distance);

            // godunov solver
            distance = std::distance(drMapping[base_cell_id][face].godunov, drMapping[cell][face].godunov);
            godunov_indices.push_back(distance);

            // flux solver
            distance = std::distance(drMapping[base_cell_id][face].fluxSolver, drMapping[cell][face].fluxSolver);
            flux_solver_indices.push_back(distance);
          }

          if (!dofs_indices.empty()) {
            ConditionalKey key(*KernelNames::neighbor_flux,
                               *FaceKinds::dynamicRupture,
                               face,
                               face_relation);
            check_key(table, key);

            table[key].variable_indices[*VariableID::dofs] = new RelativeIndices(dofs_indices, base_cell_id);
            table[key].variable_indices[*VariableID::godunov] = new RelativeIndices(godunov_indices, base_cell_id);
            table[key].variable_indices[*VariableID::fluxSolver] = new RelativeIndices(flux_solver_indices,
                                                                                       base_cell_id);
          }
        }
      }
    }
  }
}
#endif