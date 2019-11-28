/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2013-2014, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Local kernel of SeisSol.
 **/

#include "Kernels/Local.h"

#ifndef NDEBUG
#pragma message "compiling local kernel with assertions"
#endif

#include <yateto.h>


#include <cassert>
#include <stdint.h>
#include <cstring>

#include <Kernels/common.hpp>
GENERATE_HAS_MEMBER(ET)
GENERATE_HAS_MEMBER(sourceMatrix)


void seissol::kernels::Local::setGlobalData(GlobalData const* global) {
#ifndef NDEBUG
  for (unsigned stiffness = 0; stiffness < 3; ++stiffness) {
    assert( ((uintptr_t)global->stiffnessMatrices(stiffness)) % ALIGNMENT == 0 );
  }
  for (unsigned flux = 0; flux < 4; ++flux) {
    assert( ((uintptr_t)global->localChangeOfBasisMatricesTransposed(flux)) % ALIGNMENT == 0 );
    assert( ((uintptr_t)global->changeOfBasisMatrices(flux)) % ALIGNMENT == 0 );
  }
#endif

  m_volumeKernelPrototype.kDivM = global->stiffnessMatrices;
  m_localFluxKernelPrototype.rDivM = global->changeOfBasisMatrices;
  m_localFluxKernelPrototype.fMrT = global->localChangeOfBasisMatricesTransposed;
}


#ifdef ACL_DEVICE
void seissol::kernels::Local::setGlobalDataOnDevice(GlobalData const* global) {

  m_deviceVolumeKernelPrototype.kDivM = global->stiffnessMatrices;
  m_deviceLocalFluxKernelPrototype.rDivM = global->changeOfBasisMatrices;
  m_deviceLocalFluxKernelPrototype.fMrT = global->localChangeOfBasisMatricesTransposed;
}
#endif


void seissol::kernels::Local::computeIntegral(  real       i_timeIntegratedDegreesOfFreedom[tensor::I::size()],
                                                LocalData& data,
                                                LocalTmp& )
{
  // assert alignments
#ifndef NDEBUG
  assert( ((uintptr_t)i_timeIntegratedDegreesOfFreedom) % ALIGNMENT == 0 );
  assert( ((uintptr_t)data.dofs)              % ALIGNMENT == 0 );
#endif

  kernel::volume volKrnl = m_volumeKernelPrototype;
  volKrnl.Q = data.dofs;
  volKrnl.I = i_timeIntegratedDegreesOfFreedom;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
    volKrnl.star(i) = data.localIntegration.starMatrices[i];
  }

  // Optional source term
  set_ET(volKrnl, get_ptr_sourceMatrix<seissol::model::LocalData>(data.localIntegration.specific));

  kernel::localFlux lfKrnl = m_localFluxKernelPrototype;
  lfKrnl.Q = data.dofs;
  lfKrnl.I = i_timeIntegratedDegreesOfFreedom;
  lfKrnl._prefetch.I = i_timeIntegratedDegreesOfFreedom + tensor::I::size();
  lfKrnl._prefetch.Q = data.dofs + tensor::Q::size();
  
  volKrnl.execute();

  for( unsigned int face = 0; face < 4; face++ ) {
    // no element local contribution in the case of dynamic rupture boundary conditions
    if( data.cellInformation.faceTypes[face] != dynamicRupture ) {
      lfKrnl.AplusT = data.localIntegration.nApNm1[face];
      lfKrnl.execute(face);
    }
  }
}
#ifdef ACL_DEVICE
void seissol::kernels::Local::computeIntegralWithinWorkItem(real* i_timeIntegratedScratchMem,
                                                            kernels::LocalData::Loader &loader,
                                                            conditional_table_t &table,
                                                            LocalTmp& ) {
  // Volume integral
  ConditionalKey key(*KernelNames::volume);
  device_gen_code::kernel::volume volKrnl = m_deviceVolumeKernelPrototype;

  if (table.find(key) != table.end()) {
    IndexTable &index_table = table[key];
    unsigned base_cell_id = dynamic_cast<RelativeIndices*>(index_table.variable_indices[*VariableID::dofs])->cell_id;
    unsigned num_cells = index_table.variable_indices[*VariableID::dofs]->m_indices.size();

    volKrnl.num_elements = num_cells;

    auto data = loader.entry(base_cell_id);

    volKrnl.Q = data.dofs;
    //volKrnl.Q = data.debugging_dofs;
    volKrnl.Q_indices = index_table.variable_indices[*VariableID::dofs]->m_device_ptr;

    volKrnl.I = i_timeIntegratedScratchMem;
    volKrnl.I_indices = index_table.variable_indices[*VariableID::idofs]->m_device_ptr;

    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
      volKrnl.star(i) = data.localIntegrationDevice.starMatrices[i];
      volKrnl.star_indices(i) = index_table.variable_indices[*VariableID::start]->m_device_ptr;
    }

    volKrnl.execute();
  }


  // Local Flux Integral
  device_gen_code::kernel::localFlux lfKrnl = m_deviceLocalFluxKernelPrototype;
  for (unsigned face = 0; face < 4; ++face) {
    key = ConditionalKey(*KernelNames::local_flux, !FaceKinds::dynamicRupture, face);

    if (table.find(key) != table.end()) {
      IndexTable &index_table = table[key];

      unsigned base_cell_id = dynamic_cast<RelativeIndices*>(index_table.variable_indices[*VariableID::dofs])->cell_id;
      unsigned num_cells = index_table.variable_indices[*VariableID::dofs]->m_indices.size();

      auto data = loader.entry(base_cell_id);

      lfKrnl.num_elements = num_cells;

      lfKrnl.Q = data.dofs;
      //lfKrnl.Q = data.debugging_dofs;
      lfKrnl.Q_indices = index_table.variable_indices[*VariableID::dofs]->m_device_ptr;

      lfKrnl.I = i_timeIntegratedScratchMem;
      lfKrnl.I_indices = index_table.variable_indices[*VariableID::idofs]->m_device_ptr;

      lfKrnl.AplusT = data.localIntegrationDevice.nApNm1[face];
      lfKrnl.AplusT_indices = index_table.variable_indices[*VariableID::AplusT]->m_device_ptr;

      lfKrnl.execute(face);
    }
  }

}
#endif

void seissol::kernels::Local::flopsIntegral(  enum faceType const i_faceTypes[4],
                                              unsigned int        &o_nonZeroFlops,
                                              unsigned int        &o_hardwareFlops )
{
  o_nonZeroFlops = seissol::kernel::volume::NonZeroFlops;
  o_hardwareFlops = seissol::kernel::volume::HardwareFlops;

  for( unsigned int face = 0; face < 4; ++face ) {
    if( i_faceTypes[face] != dynamicRupture ) {
      o_nonZeroFlops  += seissol::kernel::localFlux::nonZeroFlops(face);
      o_hardwareFlops += seissol::kernel::localFlux::hardwareFlops(face);
    }
  }
}

unsigned seissol::kernels::Local::bytesIntegral()
{
  unsigned reals = 0;

  // star matrices load
  reals += yateto::computeFamilySize<tensor::star>();
  // flux solvers
  reals += 4 * tensor::AplusT::size();

  // DOFs write
  reals += tensor::Q::size();
  
  return reals * sizeof(real);
}
