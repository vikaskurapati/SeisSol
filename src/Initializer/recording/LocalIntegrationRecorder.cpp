#include "Recorder.h"
#include <Kernels/Interface.hpp>
#include <yateto.h>

#include "Condition.h"
#include "ConditionalTable.h"
#include "EncodingConstants.h"
#include "PointersTable.h"
#include "specific_types.h"

using namespace device;
using namespace seissol::initializers;
using namespace seissol::initializers::recording;

void LocalIntegrationRecorder::record(seissol::initializers::LTS &Handler, seissol::initializers::Layer &Layer) {

  kernels::LocalData::Loader Loader;
  Loader.load(Handler, Layer);
  LayerContainer &Container = Layer.getLayerContainer();
  auto &Table = Container.getTableReferenceToInit();

  // allocate counters and the registry
  std::unordered_map<index_t, real *> IDofsAddressRegistry{};

  size_t IDofsAddressCounter = 0;
  size_t DerivativesAddressCounter = 0;

  // get base pointers for the temp. memory
  real *IDofsScratchMem = static_cast<real *>(Layer.scratch_mem(Handler.idofs_scratch));
  real *DerivativesScratchMem = static_cast<real *>(Layer.scratch_mem(Handler.derivatives_scratch));

  // ********************** time and volume integrals **********************
  {
    if (Layer.getNumberOfCells()) {

      std::vector<real *> DofsPtrs{};
      std::vector<real *> StarPtrs{};
      std::vector<real *> IDofsPtrs{};
      std::vector<real *> DQPtrs{};

      std::vector<real *> LtsBuffers{};
      std::vector<real *> IDofsForLtsBuffers{};

      real **Derivatives = Layer.var(Handler.derivatives);
      real **Buffers = Layer.var(Handler.buffers);

      for (unsigned Cell = 0; Cell < Layer.getNumberOfCells(); ++Cell) {
        auto Data = Loader.entry(Cell);

        // dofs
        DofsPtrs.push_back(static_cast<real *>(Data.dofs));

        // idofs
        real *NextIDofPtr = &IDofsScratchMem[IDofsAddressCounter];
        bool IsBuffersProvided = ((Data.cellInformation.ltsSetup >> 8) % 2) == 1;
        bool IsLtsBuffers = ((Data.cellInformation.ltsSetup >> 10) % 2) == 1;

        if (IsBuffersProvided) {
          if (IsLtsBuffers) {
            // lts buffers may require either accumulation or overriding (in case of reset command)
            IDofsPtrs.push_back(NextIDofPtr);

            IDofsForLtsBuffers.push_back(NextIDofPtr);
            LtsBuffers.push_back(Buffers[Cell]);

            IDofsAddressRegistry[Cell] = NextIDofPtr;
            IDofsAddressCounter += tensor::I::size();
          } else {
            // gts buffers have to be always overridden
            IDofsPtrs.push_back(Buffers[Cell]);
            IDofsAddressRegistry[Cell] = Buffers[Cell];
          }
        } else {
          IDofsPtrs.push_back(NextIDofPtr);
          IDofsAddressRegistry[Cell] = NextIDofPtr;
          IDofsAddressCounter += tensor::I::size();
        }

        // stars
        StarPtrs.push_back(static_cast<real *>(Data.localIntegrationDevice.starMatrices[0]));

        // derivatives
        bool IsDerivativesProvided = ((Data.cellInformation.ltsSetup >> 9) % 2) == 1;
        if (IsDerivativesProvided) {
          DQPtrs.push_back(Derivatives[Cell]);

        } else {
          DQPtrs.push_back(&DerivativesScratchMem[DerivativesAddressCounter]);
          DerivativesAddressCounter += yateto::computeFamilySize<tensor::dQ>();
        }
      }

      ConditionalKey Key(KernelNames::time || KernelNames::volume);
      checkKey(Table, Key);

      Table[Key].m_Container[*VariableID::dofs] = new DevicePointers(DofsPtrs);
      Table[Key].m_Container[*VariableID::star] = new DevicePointers(StarPtrs);
      Table[Key].m_Container[*VariableID::idofs] = new DevicePointers(IDofsPtrs);
      Table[Key].m_Container[*VariableID::derivatives] = new DevicePointers(DQPtrs);

      Table[Key].setNotEmpty();

      if (!IDofsForLtsBuffers.empty()) {
        ConditionalKey Key(*KernelNames::time, *ComputationKind::with_lts_buffers);

        Table[Key].m_Container[*VariableID::buffers] = new DevicePointers(LtsBuffers);
        Table[Key].m_Container[*VariableID::idofs] = new DevicePointers(IDofsForLtsBuffers);
        Table[Key].setNotEmpty();
      }
    }
  }

  // ********************** local flux integral **********************
  {
    for (unsigned Face = 0; Face < 4; ++Face) {

      std::vector<real *> IDofsPtrs{};
      std::vector<real *> DofsPtrs{};
      std::vector<real *> AplusTPtrs{};

      for (unsigned Cell = 0; Cell < Layer.getNumberOfCells(); ++Cell) {

        auto Data = Loader.entry(Cell);

        // no element local contribution in the case of dynamic rupture boundary conditions
        if (Data.cellInformation.faceTypes[Face] != dynamicRupture) {
          IDofsPtrs.push_back(IDofsAddressRegistry[Cell]);
          DofsPtrs.push_back(static_cast<real *>(Data.dofs));
          AplusTPtrs.push_back(static_cast<real *>(Data.localIntegrationDevice.nApNm1[Face]));
        }
      }

      // NOTE: we can check any container, but we must check that a set is not empty!
      if (!DofsPtrs.empty()) {
        ConditionalKey Key(*KernelNames::local_flux, !FaceKinds::dynamicRupture, Face);
        checkKey(Table, Key);
        Table[Key].m_Container[*VariableID::idofs] = new DevicePointers(IDofsPtrs);
        Table[Key].m_Container[*VariableID::dofs] = new DevicePointers(DofsPtrs);
        Table[Key].m_Container[*VariableID::AplusT] = new DevicePointers(AplusTPtrs);

        Table[Key].setNotEmpty();
      }
    }
  }

  // ********************** displacements **********************
  {
    real **Displacements = Layer.var(Handler.displacements);
    std::vector<real *> IvelocitiesPtrs{};
    std::vector<real *> DisplacementsPtrs{};

    // NOTE: velocity components are between 6th and 8th columns
    const unsigned OFFSET_TO_VELOCITIES = 6 * seissol::tensor::I::Shape[0];
    for (unsigned Cell = 0; Cell < Layer.getNumberOfCells(); ++Cell) {
      if (Displacements[Cell] != nullptr) {
        real *IVelocity = &IDofsAddressRegistry[Cell][OFFSET_TO_VELOCITIES];
        IvelocitiesPtrs.push_back(IVelocity);
        DisplacementsPtrs.push_back(Displacements[Cell]);
      }
    }
    if (!DisplacementsPtrs.empty()) {
      ConditionalKey Key(*KernelNames::displacements);
      checkKey(Table, Key);
      Table[Key].m_Container[*VariableID::ivelocities] = new DevicePointers(IvelocitiesPtrs);
      Table[Key].m_Container[*VariableID::displacements] = new DevicePointers(DisplacementsPtrs);

      Table[Key].setNotEmpty();
    }
  }
}