#include "Recorder.h"
#include <Kernels/Interface.hpp>
#include <yateto.h>

using namespace device;
using namespace seissol::initializers;
using namespace seissol::initializers::recording;

void PlasticityRecorder::record(seissol::initializers::LTS &Handler, seissol::initializers::Layer &Layer) {
  kernels::LocalData::Loader Loader;
  Loader.load(Handler, Layer);
  LayerContainer &Container = Layer.getLayerContainer();
  auto &Table = Container.getTableReferenceToInit();

  // allocate counters and the registry
  std::unordered_map<index_t, real *> IDofsAddressRegistry{};

  real(*Pstrains)[7] = Layer.var(Handler.pstrain);
  size_t NodalStressTensorCounter = 0;
  real *ScratchMem = static_cast<real *>(Layer.scratch_mem(Handler.idofs_scratch));
  if (Layer.getNumberOfCells()) {
    std::vector<real *> DofsPtrs{};
    std::vector<real *> QStressNodalPtrs{};
    std::vector<real *> PstransPtrs{};

    for (unsigned Cell = 0; Cell < Layer.getNumberOfCells(); ++Cell) {
      auto Data = Loader.entry(Cell);
      DofsPtrs.push_back(static_cast<real *>(Data.dofs));
      QStressNodalPtrs.push_back(&ScratchMem[NodalStressTensorCounter]);
      NodalStressTensorCounter += tensor::QStressNodal::size();
      PstransPtrs.push_back(static_cast<real *>(Pstrains[Cell]));
    }
    if (!DofsPtrs.empty()) {
      ConditionalKey Key(*KernelNames::plasticity);
      checkKey(Table, Key);
      Table[Key].m_Container[*VariableID::dofs] = new DevicePointers(DofsPtrs);
      Table[Key].m_Container[*VariableID::NodalStressTensor] = new DevicePointers(QStressNodalPtrs);
      Table[Key].m_Container[*VariableID::Pstrains] = new DevicePointers(PstransPtrs);
      Table[Key].setNotEmpty();
    }
  }
}
