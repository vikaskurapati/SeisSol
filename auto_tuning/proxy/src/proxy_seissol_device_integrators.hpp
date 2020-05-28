/*
Copyright (c) 2015, Intel Corporation

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <generated_code/tensor.h>
#include <device.h>
#include <stdio.h>

namespace tensor = seissol::tensor;
namespace kernels = seissol::kernels;

namespace proxy {
  namespace device {

    void computeAderIntegration() {
      const DeviceInstance &Device = DeviceInstance::getInstance();
      auto& layer = m_ltsTree->child(0).child<Interior>();

      kernels::LocalData::Loader loader;
      loader.load(m_lts, layer);
      kernels::LocalTmp Tmp;

      seissol::initializers::LayerContainer &Container = layer.getLayerContainer();
      conditional_table_t &Table = Container.getConditionalTable();

      m_timeKernel.computeAderWithinWorkItem(static_cast<double>(m_timeStepWidthSimulation), Tmp, Table);
      Device.api->synchDevice();
    }

    void computeLocalWithoutAderIntegration() {
      const DeviceInstance &Device = DeviceInstance::getInstance();
      auto& layer = m_ltsTree->child(0).child<Interior>();

      kernels::LocalData::Loader loader;
      loader.load(m_lts, layer);
      kernels::LocalTmp Tmp;

      seissol::initializers::LayerContainer &Container = layer.getLayerContainer();
      conditional_table_t &Table = Container.getConditionalTable();

      m_localKernel.computeIntegralWithinWorkItem(Table, Tmp);
      Device.api->synchDevice();
    }

    void computeLocalIntegration() {
      const DeviceInstance &Device = DeviceInstance::getInstance();
      auto& layer = m_ltsTree->child(0).child<Interior>();

      kernels::LocalData::Loader loader;
      loader.load(m_lts, layer);
      kernels::LocalTmp Tmp;

      seissol::initializers::LayerContainer &Container = layer.getLayerContainer();
      conditional_table_t &Table = Container.getConditionalTable();

      m_timeKernel.computeAderWithinWorkItem(static_cast<double>(m_timeStepWidthSimulation), Tmp, Table);
      m_localKernel.computeIntegralWithinWorkItem(Table, Tmp);
      Device.api->synchDevice();
    }

    void computeNeighboringIntegration() {
      const DeviceInstance &Device = DeviceInstance::getInstance();
      auto& layer = m_ltsTree->child(0).child<Interior>();

      kernels::NeighborData::Loader loader;
      loader.load(m_lts, layer);

      seissol::initializers::LayerContainer &Container = layer.getLayerContainer();
      conditional_table_t &Table = Container.getConditionalTable();

      seissol::kernels::TimeCommon::computeIntegralsForWorkItem(m_timeKernel,
                                                                0.0,
                                                                static_cast<double>(m_timeStepWidthSimulation),
                                                                Table);
      m_neighborKernel.computeNeighborsIntegralWithinWorkItem(Table);
      Device.api->synchDevice();
    }

    void computeDynRupGodunovState() {
    }
  }  // namespace device
}  // namespace proxy
