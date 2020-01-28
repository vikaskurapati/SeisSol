/******************************************************************************
** Copyright (c) 2015, Intel Corporation                                     **
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
 * @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2013, SeisSol Group
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
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Aligned memory allocation.
 **/

#ifdef ACL_DEVICE
#include "device_utils.h"
#endif

#include <algorithm>
#include "MemoryAllocator.h"

#include <utils/logger.h>

void* seissol::memory::allocate(size_t i_size, size_t i_alignment, enum Memkind i_memkind)
{
    void* l_ptrBuffer;
    bool error = false;

    /* handle zero allocation */
    if ( i_size == 0 ) {
      //logWarning() << "allocation of size 0 requested, returning NULL; (alignment: " << i_alignment << ", memkind: " << i_memkind << ").";
      l_ptrBuffer = NULL;
      return l_ptrBuffer;
    }

#if defined(USE_MEMKIND) || defined(ACL_DEVICE)
    if( i_memkind == 0 ) {
#endif
      if (i_alignment % (sizeof(void*)) != 0) {
        l_ptrBuffer = malloc( i_size );
        error = (l_ptrBuffer == NULL);
      } else {
        error = (posix_memalign( &l_ptrBuffer, i_alignment, i_size ) != 0);
      }
#ifdef USE_MEMKIND
    }
    else if (i_memkind == HighBandwidth) {
      if (i_alignment % (sizeof(void*)) != 0) {
        l_ptrBuffer = hbw_malloc( i_size );
        error = (l_ptrBuffer == NULL);

      } else {
        error = (hbw_posix_memalign( &l_ptrBuffer, i_alignment, i_size ) != 0);
      }
    }
#endif

#ifdef ACL_DEVICE
    }
    else if (i_memkind == DeviceGlobalMemory){
      l_ptrBuffer = device_malloc(i_size);
    }
    else if (i_memkind == DeviceUnifiedMemory){
      l_ptrBuffer = device_malloc_unified(i_size);
    }
    else if (i_memkind == PinnedMemory) {
      l_ptrBuffer = device_malloc_pinned(i_size);
    }
#endif

#if defined(USE_MEMKIND) || defined(ACL_DEVICE)
    else {
        logError() << "unknown memkind type used (" << i_memkind << "). Please, refer to the documentation";
    }
#endif

    if (error) {
        logError() << "The malloc failed (bytes: " << i_size << ", alignment: " << i_alignment << ", memkind: " << i_memkind << ").";
    }

    return l_ptrBuffer;
}

void seissol::memory::free(void* i_pointer, enum Memkind i_memkind) {
#if defined(USE_MEMKIND) || defined(ACL_DEVICE)
    if (i_memkind == Standard) {
#endif
        ::free(i_pointer);
#ifdef USE_MEMKIND
        }
        else if (i_memkind == HighBandwidth) {
          hbw_free(i_pointer);
        }
#endif

#ifdef ACL_DEVICE
    }
    else if ((i_memkind == DeviceGlobalMemory) || (i_memkind == DeviceUnifiedMemory)) {
        device_free(i_pointer);
    }
    else if (i_memkind == PinnedMemory) {
        device_free_pinned(i_pointer);
    }
#endif

#if defined(USE_MEMKIND) || defined(ACL_DEVICE)
    else {
        logError() << "unknown memkind type used (" << i_memkind << "). Please, refer to the documentation";
    }
#endif
}

void seissol::memory::printMemoryAlignment( std::vector< std::vector<unsigned long long> > i_memoryAlignment ) {
  logDebug() << "printing memory alignment per struct";
  for( unsigned long long l_i = 0; l_i < i_memoryAlignment.size(); l_i++ ) {
    logDebug() << i_memoryAlignment[l_i][0] << ", " << i_memoryAlignment[l_i][1];
  }
}


seissol::memory::ManagedAllocator::~ManagedAllocator()
{
  for (AddressVector::const_iterator it = m_dataMemoryAddresses.begin(); it != m_dataMemoryAddresses.end(); ++it) {
      if (it->first == seissol::memory::DeviceGlobalMemory) {
          std::cout << "ERROR::device global mem. in a destructor" << std::endl;
          throw;
      }


      if (it->first == seissol::memory::DeviceUnifiedMemory) {
          std::cout << "ERROR::uniformed mem. in a destructor" << std::endl;
          throw;
      }

      seissol::memory::free(it->second, it->first);
  }

  // reset memory vectors
  m_dataMemoryAddresses.clear();
}

void* seissol::memory::ManagedAllocator::allocateMemory( size_t i_size, size_t i_alignment, enum Memkind i_memkind )
{
  void* l_ptrBuffer = seissol::memory::allocate(i_size, i_alignment, i_memkind);
  m_dataMemoryAddresses.push_back( Address(i_memkind, l_ptrBuffer) );
  return l_ptrBuffer;
}


void seissol::memory::ManagedAllocator::deallocateMemory(void *i_ptr, enum Memkind i_memkind) {
    auto address = std::make_pair(i_memkind, i_ptr);
    auto it = std::find(m_dataMemoryAddresses.begin(), m_dataMemoryAddresses.end(), address);
    if (it != std::end(m_dataMemoryAddresses)) {
        seissol::memory::free(it->second, it->first);
        m_dataMemoryAddresses.erase(it);
    }
}
