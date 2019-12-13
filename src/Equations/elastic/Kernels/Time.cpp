/******************************************************************************
** Copyright (c) 2014-2015, Intel Corporation                                **
** All rights reserved.                                                      **
**                                                                           **
** Redistribution and use in source and binary forms, with or without        **
** modification, are permitted provided that the following conditions        **
** are met:                                                                  **
** 1. Redistributions of source code must retain the above copyright         **
**    notice, this list of conditions and the following disclaimer.          **
** 2. Redistributions in binary form must reproduce the above copyright      **
**    notice, this list of conditions and the following disclaimer in the    **
**    documentation and/or other materials provided with the distribution.   **
** 3. Neither the name of the copyright holder nor the names of its          **
**    contributors may be used to endorse or promote products derived        **
**    from this software without specific prior written permission.          **
**                                                                           **
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       **
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT         **
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR     **
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT      **
** HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,    **
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED  **
** TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    **
** PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    **
** LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      **
** NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        **
** SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              **
******************************************************************************/
/* Alexander Heinecke (Intel Corp.)
******************************************************************************/
/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2013-2015, SeisSol Group
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
 * Time kernel of SeisSol.
 **/

#include "Kernels/TimeBase.h"
#include "Kernels/Time.h"

#ifndef NDEBUG
#pragma message "compiling time kernel with assertions"
extern long long libxsmm_num_total_flops;
#endif

#include <Kernels/common.hpp>
#include <Kernels/denseMatrixOps.hpp>

#include <cstring>
#include <cassert>
#include <stdint.h>
#include <omp.h>

#include <yateto.h>

GENERATE_HAS_MEMBER(ET)
GENERATE_HAS_MEMBER(sourceMatrix)

seissol::kernels::TimeBase::TimeBase() {
  m_derivativesOffsets[0] = 0;
  for (int order = 0; order < CONVERGENCE_ORDER; ++order) {
    if (order > 0) {
      m_derivativesOffsets[order] = tensor::dQ::size(order-1) + m_derivativesOffsets[order-1];
    }
  }
}

void seissol::kernels::Time::setGlobalData(GlobalData const* global) {
  assert( ((uintptr_t)global->stiffnessMatricesTransposed(0)) % ALIGNMENT == 0 );
  assert( ((uintptr_t)global->stiffnessMatricesTransposed(1)) % ALIGNMENT == 0 );
  assert( ((uintptr_t)global->stiffnessMatricesTransposed(2)) % ALIGNMENT == 0 );

  m_krnlPrototype.kDivMT = global->stiffnessMatricesTransposed;
}

#ifdef ACL_DEVICE
void seissol::kernels::Time::setGlobalDataOnDevice(GlobalData const* global) {
  assert( ((uintptr_t)global->stiffnessMatricesTransposed(0)) % ALIGNMENT == 0 );
  assert( ((uintptr_t)global->stiffnessMatricesTransposed(1)) % ALIGNMENT == 0 );
  assert( ((uintptr_t)global->stiffnessMatricesTransposed(2)) % ALIGNMENT == 0 );

  m_DeviceKrnlPrototype.kDivMT = global->stiffnessMatricesTransposed;
}
#endif

void seissol::kernels::Time::computeAder( double                      i_timeStepWidth,
                                          LocalData&                  data,
                                          LocalTmp&                   tmp,
                                          real                        o_timeIntegrated[tensor::I::size()],
                                          real*                       o_timeDerivatives )
{
  /*
   * assert alignments.
   */
  assert( ((uintptr_t)data.dofs)              % ALIGNMENT == 0 );
  assert( ((uintptr_t)o_timeIntegrated )      % ALIGNMENT == 0 );
  assert( ((uintptr_t)o_timeDerivatives)      % ALIGNMENT == 0 || o_timeDerivatives == NULL );

  /*
   * compute ADER scheme.
   */
  // temporary result
  real temporaryBuffer[yateto::computeFamilySize<tensor::dQ>()] __attribute__((aligned(PAGESIZE_STACK)));
  real* derivativesBuffer = (o_timeDerivatives != nullptr) ? o_timeDerivatives : temporaryBuffer;

  kernel::derivative krnl = m_krnlPrototype;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
    krnl.star(i) = data.localIntegration.starMatrices[i];
  }

  // Optional source term
  set_ET(krnl, get_ptr_sourceMatrix<seissol::model::LocalData>(data.localIntegration.specific));

  krnl.dQ(0) = const_cast<real*>(data.dofs);
  for (unsigned i = 1; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    krnl.dQ(i) = derivativesBuffer + m_derivativesOffsets[i];
  }

  kernel::derivativeTaylorExpansion intKrnl;
  intKrnl.I = o_timeIntegrated;
  intKrnl.dQ(0) = data.dofs;
  for (unsigned i = 1; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = derivativesBuffer + m_derivativesOffsets[i];
  }
  
  // powers in the taylor-series expansion
  intKrnl.power = i_timeStepWidth;

  intKrnl.execute0();

  // stream out frist derivative (order 0)
  if (o_timeDerivatives != nullptr) {
    streamstore(tensor::dQ::size(0), data.dofs, o_timeDerivatives);
  }
  
  for (unsigned der = 1; der < CONVERGENCE_ORDER; ++der) {
    krnl.execute(der);

    // update scalar for this derivative
    intKrnl.power *= i_timeStepWidth / real(der+1);    
    intKrnl.execute(der);
  }
}

#ifdef ACL_DEVICE
#include <device_utils.h>
void seissol::kernels::Time::computeAderWithinWorkItem(double i_timeStepWidth,
                                                       LocalTmp& tmp,
                                                       kernels::LocalData::Loader &loader,
                                                       conditional_table_t &table,
                                                       real *o_timeIntegratedScratchMem,
                                                       real *o_timeDerivativesScratchMem,
                                                       real **derivatives) {

  // NOTE: wod - without derivatives; wd - with derivatives
  // NOTE: base prefix stands for the first instance of a sequence
  //       Example 1: base_cell_id - the first cell in a sequence from which indices were recorded
  //       Example 2: base_star_indices - the fist star matrix, start(0), for which the indices matrices were recorded
  //                                      indices for start(1) and start(2) are computed in run-time

  //kernel::derivative derivativesKrnl = m_krnlPrototype;
  device_gen_code::kernel::derivative derivativesKrnl = m_DeviceKrnlPrototype;
  device_gen_code::kernel::derivativeTaylorExpansion intKrnl;

  DeviceTemporaryMemoryMenager &tmp_mem_manager = DeviceTemporaryMemoryMenager::get_instance();

  // compute cells which do NOT require to have their derivatives
  ConditionalKey key(*KernelNames::time, *ComputationKind::without_derivatives);
  if(table.find(key) != table.end()) {

    IndexTable &index_table = table[key];
    unsigned base_cell_id = dynamic_cast<RelativeIndices*>(index_table.variable_indices[*VariableID::dofs])->cell_id;
    unsigned num_cells = index_table.variable_indices[*VariableID::dofs]->m_indices.size();

    derivativesKrnl.num_elements = num_cells;
    intKrnl.num_elements = num_cells;

    auto data = loader.entry(base_cell_id);

    // attach idofs
    intKrnl.I = o_timeIntegratedScratchMem;
    intKrnl.I_indices = index_table.variable_indices[*VariableID::idofs]->m_device_ptr;

    // attach start matrices
    // NOTE: all start matrices are equally shifted, i.e. the indices can be reused for all tree star matrices
    //       only the base data pointer must be shifted
    //unsigned *start_indices = index_table.variable_indices[*VariableID::start]->m_device_ptr;
    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
      derivativesKrnl.star(i) = data.localIntegrationDevice.starMatrices[i];
      derivativesKrnl.star_indices(i) = index_table.variable_indices[*VariableID::start]->m_device_ptr;
    }


    // attach derivatives
    // NOTE: the very first derivative is going to come from dofs themselves

    // zero derivative
    derivativesKrnl.dQ(0) = data.dofs;
    intKrnl.dQ(0) = data.dofs;

    derivativesKrnl.dQ_indices(0) = index_table.variable_indices[*VariableID::dofs]->m_device_ptr;
    intKrnl.dQ_indices(0) = index_table.variable_indices[*VariableID::dofs]->m_device_ptr;

    // other derivatives: from 1 to #order
    unsigned offset = 0;
    for (unsigned i = 1; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
      derivativesKrnl.dQ(i) = &o_timeDerivativesScratchMem[offset];
      intKrnl.dQ(i) = &o_timeDerivativesScratchMem[offset];

      derivativesKrnl.dQ_indices(i) = index_table.variable_indices[*VariableID::derivatives]->m_device_ptr;
      intKrnl.dQ_indices(i) = index_table.variable_indices[*VariableID::derivatives]->m_device_ptr;

      offset += tensor::dQ::size(i);
    }

    intKrnl.power = i_timeStepWidth;
    intKrnl.execute0();

    for (unsigned der = 1; der < CONVERGENCE_ORDER; ++der) {
      derivativesKrnl.execute(der);

      // update scalar for this derivative
      intKrnl.power *= i_timeStepWidth / real(der+1);
      intKrnl.execute(der);
    }
  }


  // compute cells which HAVE their own derivatives
  key = ConditionalKey(*KernelNames::time, *ComputationKind::with_derivatives);
  if(table.find(key) != table.end()) {
    IndexTable &index_table = table[key];
    unsigned base_cell_id = dynamic_cast<RelativeIndices *>(index_table.variable_indices[*VariableID::dofs])->cell_id;
    unsigned num_cells = index_table.variable_indices[*VariableID::dofs]->m_indices.size();

    derivativesKrnl.num_elements = num_cells;
    intKrnl.num_elements = num_cells;

    auto data = loader.entry(base_cell_id);

    // attach idofs
    intKrnl.I = o_timeIntegratedScratchMem;
    intKrnl.I_indices = index_table.variable_indices[*VariableID::idofs]->m_device_ptr;

    // attach start matrices
    // NOTE: all start matrices are equally shifted, i.e. the indices can be reused for all tree star matrices
    //       only the base data pointer must be shifted
    //unsigned *start_indices = index_table.variable_indices[*VariableID::start]->m_device_ptr;
    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
      derivativesKrnl.star(i) = data.localIntegrationDevice.starMatrices[i];
      derivativesKrnl.star_indices(i) = index_table.variable_indices[*VariableID::start]->m_device_ptr;
    }

    // stream dofs to the zero derivative
    device_stream_data(data.dofs,
                       index_table.variable_indices[*VariableID::dofs]->m_device_ptr,
                       derivatives[base_cell_id],
                       index_table.variable_indices[*VariableID::derivatives]->m_device_ptr,
                       tensor::Q::size(),
                       num_cells);

    unsigned offset = 0;
    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {

      derivativesKrnl.dQ(i) = &derivatives[base_cell_id][offset];
      derivativesKrnl.dQ_indices(i) = index_table.variable_indices[*VariableID::derivatives]->m_device_ptr;

      intKrnl.dQ(i) = &derivatives[base_cell_id][offset];
      intKrnl.dQ_indices(i) = index_table.variable_indices[*VariableID::derivatives]->m_device_ptr;

      offset += tensor::dQ::size(i);
    }

    intKrnl.power = i_timeStepWidth;
    intKrnl.execute0();

    for (unsigned der = 1; der < CONVERGENCE_ORDER; ++der) {
      derivativesKrnl.execute(der);

      // update scalar for this derivative
      intKrnl.power *= i_timeStepWidth / real(der + 1);
      intKrnl.execute(der);
    }
  }
}
#endif

void seissol::kernels::Time::flopsAder( unsigned int        &o_nonZeroFlops,
                                        unsigned int        &o_hardwareFlops ) {
  // reset flops
  o_nonZeroFlops = 0; o_hardwareFlops =0;

  // initialization
  o_nonZeroFlops  += kernel::derivativeTaylorExpansion::nonZeroFlops(0);
  o_hardwareFlops += kernel::derivativeTaylorExpansion::hardwareFlops(0);

  // interate over derivatives
  for( unsigned l_derivative = 1; l_derivative < CONVERGENCE_ORDER; l_derivative++ ) {
    o_nonZeroFlops  += kernel::derivative::nonZeroFlops(l_derivative);
    o_hardwareFlops += kernel::derivative::hardwareFlops(l_derivative);

    // update of time integrated DOFs
    o_nonZeroFlops  += kernel::derivativeTaylorExpansion::nonZeroFlops(l_derivative);
    o_hardwareFlops += kernel::derivativeTaylorExpansion::hardwareFlops(l_derivative);
  }

}

unsigned seissol::kernels::Time::bytesAder()
{
  unsigned reals = 0;
  
  // DOFs load, tDOFs load, tDOFs write
  reals += tensor::Q::size() + 2 * tensor::I::size();
  // star matrices, source matrix
  reals += yateto::computeFamilySize<tensor::star>();
           
  /// \todo incorporate derivatives

  return reals * sizeof(real);
}

void seissol::kernels::Time::computeIntegral( double                            i_expansionPoint,
                                              double                            i_integrationStart,
                                              double                            i_integrationEnd,
                                              const real*                       i_timeDerivatives,
                                              real                              o_timeIntegrated[tensor::I::size()] )
{
  /*
   * assert alignments.
   */
  assert( ((uintptr_t)i_timeDerivatives)  % ALIGNMENT == 0 );
  assert( ((uintptr_t)o_timeIntegrated)   % ALIGNMENT == 0 );

  // assert that this is a forwared integration in time
  assert( i_integrationStart + (real) 1.E-10 > i_expansionPoint   );
  assert( i_integrationEnd                   > i_integrationStart );

  /*
   * compute time integral.
   */
  // compute lengths of integration intervals
  real l_deltaTLower = i_integrationStart - i_expansionPoint;
  real l_deltaTUpper = i_integrationEnd   - i_expansionPoint;

  // initialization of scalars in the taylor series expansion (0th term)
  real l_firstTerm  = (real) 1;
  real l_secondTerm = (real) 1;
  real l_factorial  = (real) 1;
  
  kernel::derivativeTaylorExpansion intKrnl;
  intKrnl.I = o_timeIntegrated;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = i_timeDerivatives + m_derivativesOffsets[i];
  }
 
  // iterate over time derivatives
  for(int der = 0; der < CONVERGENCE_ORDER; ++der ) {
    l_firstTerm  *= l_deltaTUpper;
    l_secondTerm *= l_deltaTLower;
    l_factorial  *= (real)(der+1);

    intKrnl.power  = l_firstTerm - l_secondTerm;
    intKrnl.power /= l_factorial;

    intKrnl.execute(der);
  }
}

void seissol::kernels::Time::computeTaylorExpansion( real         time,
                                                     real         expansionPoint,
                                                     real const*  timeDerivatives,
                                                     real         timeEvaluated[tensor::Q::size()] ) {
  /*
   * assert alignments.
   */
  assert( ((uintptr_t)timeDerivatives)  % ALIGNMENT == 0 );
  assert( ((uintptr_t)timeEvaluated)    % ALIGNMENT == 0 );

  // assert that this is a forward evaluation in time
  assert( time >= expansionPoint );

  real deltaT = time - expansionPoint;

  static_assert(tensor::I::size() == tensor::Q::size(), "Sizes of tensors I and Q must match");

  kernel::derivativeTaylorExpansion intKrnl;
  intKrnl.I = timeEvaluated;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = timeDerivatives + m_derivativesOffsets[i];
  }
  intKrnl.power = 1.0;
 
  // iterate over time derivatives
  for(int derivative = 0; derivative < CONVERGENCE_ORDER; ++derivative) {
    intKrnl.execute(derivative);
    intKrnl.power *= deltaT / real(derivative+1);
  }
}

void seissol::kernels::Time::flopsTaylorExpansion(long long& nonZeroFlops, long long& hardwareFlops) {
  // reset flops
  nonZeroFlops = 0; hardwareFlops = 0;

  // interate over derivatives
  for (unsigned der = 0; der < CONVERGENCE_ORDER; ++der) {
    nonZeroFlops  += kernel::derivativeTaylorExpansion::nonZeroFlops(der);
    hardwareFlops += kernel::derivativeTaylorExpansion::hardwareFlops(der);
  }
}
