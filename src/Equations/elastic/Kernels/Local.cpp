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

#ifdef  ACL_DEVICE
void seissol::kernels::Local::computeIntegralWithinWorkItem(conditional_table_t &Table,
                                                            LocalTmp& ) {

  // Volume integral
  ConditionalKey Key(KernelNames::time || KernelNames::volume);
  kernel::gpu_volume VolKrnl = m_deviceVolumeKernelPrototype;
  kernel::gpu_localFlux LocalFluxKrnl = m_deviceLocalFluxKernelPrototype;

  constexpr unsigned MAX_TMP_MEM = (VolKrnl.TmpMaxMemRequiredInBytes > LocalFluxKrnl.TmpMaxMemRequiredInBytes) \
                                   ? VolKrnl.TmpMaxMemRequiredInBytes : LocalFluxKrnl.TmpMaxMemRequiredInBytes;
  real* TmpMem = nullptr;

  if (Table.find(Key) != Table.end()) {
    PointersTable &Entry = Table[Key];

    unsigned MaxNumElements = (Entry.m_Container[*VariableID::dofs])->getSize();
    VolKrnl.numElements = MaxNumElements;

    // volume kernel always contains more elements than any local one
    TmpMem = (real*)(m_Device.api->getStackMemory(MAX_TMP_MEM * MaxNumElements));

    VolKrnl.Q = (Entry.m_Container[*VariableID::dofs])->getPointers();
    VolKrnl.I = const_cast<const real **>((Entry.m_Container[*VariableID::idofs])->getPointers());

    unsigned StarOffset = 0;
    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
      VolKrnl.star(i) = const_cast<const real **>((Entry.m_Container[*VariableID::star])->getPointers());
      VolKrnl.extraOffset_star(i) = StarOffset;
      StarOffset += tensor::star::size(i);
    }
    VolKrnl.linearAllocator.initialize(TmpMem);
    VolKrnl.execute();
  }

  // Local Flux Integral
  for (unsigned Face = 0; Face < 4; ++Face) {
    Key = ConditionalKey(*KernelNames::local_flux, !FaceKinds::dynamicRupture, Face);

    if (Table.find(Key) != Table.end()) {
      PointersTable &Entry = Table[Key];
      LocalFluxKrnl.numElements = Entry.m_Container[*VariableID::dofs]->getSize();
      LocalFluxKrnl.Q = (Entry.m_Container[*VariableID::dofs])->getPointers();
      LocalFluxKrnl.I = const_cast<const real **>((Entry.m_Container[*VariableID::idofs])->getPointers());
      LocalFluxKrnl.AplusT = const_cast<const real **>(Entry.m_Container[*VariableID::AplusT]->getPointers());
      LocalFluxKrnl.linearAllocator.initialize(TmpMem);
      LocalFluxKrnl.execute(Face);
    }
  }
  if (TmpMem != nullptr) {
    m_Device.api->popStackMemory();
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
