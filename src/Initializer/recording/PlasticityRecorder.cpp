#include "Recorder.h"
#include <Kernels/Interface.hpp>
#include <yateto.h>

using namespace device;
using namespace seissol::initializers;
using namespace seissol::initializers::recording;

void PlasticityRecorder::record(seissol::initializers::LTS &handler, seissol::initializers::Layer &layer) {
  kernels::LocalData::Loader loader;
  loader.load(handler, layer);
  LayerContainer &container = layer.getLayerContainer();
  auto &table = container.getTableReferenceToInit();

  // allocate counters and the registry
  std::unordered_map<index_t, real*> idofs_address_registry{};

  size_t NodalStressTensorCounter = 0;
  real* ScratchMem = static_cast<real*>(layer.scratch_mem(handler.idofs_scratch));
  if (layer.getNumberOfCells()) {
    std::vector<real*> DofsPtrs{};
    std::vector<real*> QStressNodalPtrs{};

    for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
      auto data = loader.entry(cell);
      DofsPtrs.push_back(static_cast<real *>(data.dofs));
      QStressNodalPtrs.push_back(&ScratchMem[NodalStressTensorCounter]);
      NodalStressTensorCounter += NodalStressTensorCounter * tensor::QStressNodal::size();
    }
    if (!DofsPtrs.empty()) {
      ConditionalKey key(*KernelNames::plasticity);
      checkKey(table, key);
      table[key].container[*VariableID::dofs] = new DevicePointers(DofsPtrs);
      table[key].container[*VariableID::NodalStressTensor] = new DevicePointers(QStressNodalPtrs);
      table[key].set_not_empty_flag();
    }
  }
}
