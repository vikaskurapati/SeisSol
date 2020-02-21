#include "Recorder.h"
#include <Kernels/Interface.hpp>
#include <yateto.h>

#include "EncodingConstants.h"
#include "Condition.h"
#include "ConditionalTable.h"
#include "PointersTable.h"
#include "specific_types.h"

using namespace device;
using namespace seissol::initializers;
using namespace seissol::initializers::recording;

void LocalIntegrationRecorder::record(seissol::initializers::LTS &handler, seissol::initializers::Layer &layer) {

  kernels::LocalData::Loader loader;
  loader.load(handler, layer);
  LayerContainer &container = layer.getLayerContainer();
  auto &table = container.getTableReferenceToInit();

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
      checkKey(table, key);

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
        checkKey(table, key);
        table[key].container[*VariableID::idofs] = new DevicePointers(idofs_ptrs);
        table[key].container[*VariableID::dofs] = new DevicePointers(dofs_ptrs);
        table[key].container[*VariableID::AplusT] = new DevicePointers(AplusT_ptrs);

        table[key].set_not_empty_flag();
      }
    }
  }

  // *************************************** displacements ***************************************
  {
    real** displacements = layer.var(handler.displacements);
    std::vector<real*> ivelocities_ptrs{};
    std::vector<real*> displacements_ptrs{};

    // NOTE: velocity components are between 6th and 8th columns
    const unsigned OffsetToVelocities = 6 * seissol::tensor::I::Shape[0];
    for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
      if (displacements[cell] != nullptr) {
        real *ivelocity = &idofs_address_registry[cell][OffsetToVelocities];
        ivelocities_ptrs.push_back(ivelocity);
        displacements_ptrs.push_back(displacements[cell]);
      }
    }
    if (!displacements_ptrs.empty()) {
      ConditionalKey key(*KernelNames::displacements);
      checkKey(table, key);
      table[key].container[*VariableID::ivelocities] = new DevicePointers(ivelocities_ptrs);
      table[key].container[*VariableID::displacements] = new DevicePointers(displacements_ptrs);

      table[key].set_not_empty_flag();
    }
  }
}