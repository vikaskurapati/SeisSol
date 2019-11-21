/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2017, SeisSol Group
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
 **/

#include "GlobalData.h"
#include <generated_code/init.h>
#include <yateto.h>

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace init = seissol::init;

void seissol::initializers::initializeGlobalData(GlobalData& globalData, memory::ManagedAllocator& memoryAllocator, enum seissol::memory::Memkind memkind)
{
  // We ensure that global matrices always start at an aligned memory address,
  // such that mixed cases with aligned and non-aligned global matrices do also work.  

  unsigned globalMatrixMemSize = 0;
  globalMatrixMemSize += yateto::computeFamilySize<init::kDivM>(yateto::alignedReals<real>(ALIGNMENT));
  globalMatrixMemSize += yateto::computeFamilySize<init::kDivMT>(yateto::alignedReals<real>(ALIGNMENT));
  globalMatrixMemSize += yateto::computeFamilySize<init::rDivM>(yateto::alignedReals<real>(ALIGNMENT));
  globalMatrixMemSize += yateto::computeFamilySize<init::rT>(yateto::alignedReals<real>(ALIGNMENT));
  globalMatrixMemSize += yateto::computeFamilySize<init::fMrT>(yateto::alignedReals<real>(ALIGNMENT));
  globalMatrixMemSize += yateto::computeFamilySize<init::fP>(yateto::alignedReals<real>(ALIGNMENT));
  globalMatrixMemSize += yateto::alignedUpper(tensor::evalAtQP::size(),  yateto::alignedReals<real>(ALIGNMENT));
  globalMatrixMemSize += yateto::alignedUpper(tensor::projectQP::size(), yateto::alignedReals<real>(ALIGNMENT));
  
  real* globalMatrixMem = static_cast<real*>(memoryAllocator.allocateMemory( globalMatrixMemSize * sizeof(real), PAGESIZE_HEAP, memkind ));

  real* globalMatrixMemPtr = globalMatrixMem;
  yateto::CopyManager<real> copy_manager;
  copy_manager.copyFamilyToMemAndSetPtr<init::kDivMT>(globalMatrixMemPtr, globalData.stiffnessMatricesTransposed, ALIGNMENT);
  copy_manager.copyFamilyToMemAndSetPtr<init::kDivM>(globalMatrixMemPtr, globalData.stiffnessMatrices, ALIGNMENT);
  copy_manager.copyFamilyToMemAndSetPtr<init::rDivM>(globalMatrixMemPtr, globalData.changeOfBasisMatrices, ALIGNMENT);
  copy_manager.copyFamilyToMemAndSetPtr<init::rT>(globalMatrixMemPtr, globalData.neighbourChangeOfBasisMatricesTransposed, ALIGNMENT);
  copy_manager.copyFamilyToMemAndSetPtr<init::fMrT>(globalMatrixMemPtr, globalData.localChangeOfBasisMatricesTransposed, ALIGNMENT);
  copy_manager.copyFamilyToMemAndSetPtr<init::fP>(globalMatrixMemPtr, globalData.neighbourFluxMatrices, ALIGNMENT);
  copy_manager.copyTensorToMemAndSetPtr<init::evalAtQP>(globalMatrixMemPtr, globalData.evalAtQPMatrix, ALIGNMENT);
  copy_manager.copyTensorToMemAndSetPtr<init::projectQP>(globalMatrixMemPtr, globalData.projectQPMatrix, ALIGNMENT);

  assert(globalMatrixMemPtr == globalMatrixMem + globalMatrixMemSize);

  // @TODO Integrate this step into the code generator
  for (unsigned transposedStiffness = 0; transposedStiffness < 3; ++transposedStiffness) {
    real* matrix = const_cast<real*>(globalData.stiffnessMatricesTransposed(transposedStiffness));
    for (unsigned i = 0; i < init::kDivMT::size(transposedStiffness); ++i) {
      matrix[i] *= -1.0;
    }
  }

  // Dynamic Rupture global matrices
  unsigned drGlobalMatrixMemSize = 0;
  drGlobalMatrixMemSize += yateto::computeFamilySize<init::V3mTo2nTWDivM>(yateto::alignedReals<real>(ALIGNMENT));
  drGlobalMatrixMemSize += yateto::computeFamilySize<init::V3mTo2n>(yateto::alignedReals<real>(ALIGNMENT));
  
  real* drGlobalMatrixMem = static_cast<real*>(memoryAllocator.allocateMemory( drGlobalMatrixMemSize  * sizeof(real), PAGESIZE_HEAP, memkind ));
  
  real* drGlobalMatrixMemPtr = drGlobalMatrixMem;
  copy_manager.copyFamilyToMemAndSetPtr<init::V3mTo2nTWDivM>(drGlobalMatrixMemPtr, globalData.nodalFluxMatrices, ALIGNMENT);
  copy_manager.copyFamilyToMemAndSetPtr<init::V3mTo2n>(drGlobalMatrixMemPtr, globalData.faceToNodalMatrices, ALIGNMENT);

  assert(drGlobalMatrixMemPtr == drGlobalMatrixMem + drGlobalMatrixMemSize);

  // Plasticity global matrices
  unsigned plasticityGlobalMatrixMemSize = 0;
  plasticityGlobalMatrixMemSize += yateto::alignedUpper(tensor::v::size(),    yateto::alignedReals<real>(ALIGNMENT));
  plasticityGlobalMatrixMemSize += yateto::alignedUpper(tensor::vInv::size(), yateto::alignedReals<real>(ALIGNMENT));

  real* plasticityGlobalMatrixMem = static_cast<real*>(memoryAllocator.allocateMemory( plasticityGlobalMatrixMemSize * sizeof(real), PAGESIZE_HEAP, memkind ));
  
  real* plasticityGlobalMatrixMemPtr = plasticityGlobalMatrixMem;
  copy_manager.copyTensorToMemAndSetPtr<init::v>(plasticityGlobalMatrixMemPtr, globalData.vandermondeMatrix, ALIGNMENT);
  copy_manager.copyTensorToMemAndSetPtr<init::vInv>(plasticityGlobalMatrixMemPtr, globalData.vandermondeMatrixInverse, ALIGNMENT);


  assert(plasticityGlobalMatrixMemPtr == plasticityGlobalMatrixMem + plasticityGlobalMatrixMemSize);
  
  // thread-local LTS integration buffers  
  int l_numberOfThreads = 1;
#ifdef _OPENMP
  l_numberOfThreads = omp_get_max_threads();
#endif
  real* integrationBufferLTS = (real*) memoryAllocator.allocateMemory( l_numberOfThreads*(4*tensor::I::size())*sizeof(real), PAGESIZE_STACK, memkind ) ;

  // initialize w.r.t. NUMA
#ifdef _OPENMP
  #pragma omp parallel
  {
    size_t l_threadOffset = omp_get_thread_num()*(4*tensor::I::size());
#else
    size_t l_threadOffset = 0;
#endif
    for ( unsigned int l_dof = 0; l_dof < (4*tensor::I::size()); l_dof++ ) {
      integrationBufferLTS[l_dof + l_threadOffset] = (real)0.0;
    }
#ifdef _OPENMP
  }
#endif
  
  globalData.integrationBufferLTS = integrationBufferLTS;
}


#ifdef ACL_DEVICE
#include <device_utils.h>
// Allocates and inits global data structures on a device (gpu)
void seissol::initializers::initializeGlobalDataOnDevice(GlobalDataOnDevice& globalData,
                                                         memory::ManagedAllocator& memoryAllocator) {

    seissol::memory::Memkind memkind = seissol::memory::Memkind::DeviceGlobalMemory;
    // We ensure that global matrices always start at an aligned memory address,
    // such that mixed cases with aligned and non-aligned global matrices do also work.

    // compute a memory size needed to store global tensors and matrices
    // NOTE: the memory size includes padding for efficient vectorization
    unsigned globalMatrixMemSize = 0;
    const unsigned DEVICE_ALIGNMENT_IN_BYTES = 256; // TODO: move it out of here (it must be known during compilation)
    globalMatrixMemSize += yateto::computeFamilySize<init::kDivM>(yateto::alignedReals<real>(DEVICE_ALIGNMENT_IN_BYTES));
    globalMatrixMemSize += yateto::computeFamilySize<init::kDivMT>(yateto::alignedReals<real>(DEVICE_ALIGNMENT_IN_BYTES));
    globalMatrixMemSize += yateto::computeFamilySize<init::rDivM>(yateto::alignedReals<real>(DEVICE_ALIGNMENT_IN_BYTES));
    globalMatrixMemSize += yateto::computeFamilySize<init::rT>(yateto::alignedReals<real>(DEVICE_ALIGNMENT_IN_BYTES));
    globalMatrixMemSize += yateto::computeFamilySize<init::fMrT>(yateto::alignedReals<real>(DEVICE_ALIGNMENT_IN_BYTES));
    globalMatrixMemSize += yateto::computeFamilySize<init::fP>(yateto::alignedReals<real>(DEVICE_ALIGNMENT_IN_BYTES));
    globalMatrixMemSize += yateto::alignedUpper(tensor::evalAtQP::size(),  yateto::alignedReals<real>(DEVICE_ALIGNMENT_IN_BYTES));
    globalMatrixMemSize += yateto::alignedUpper(tensor::projectQP::size(), yateto::alignedReals<real>(DEVICE_ALIGNMENT_IN_BYTES));

    // allocate memory
    const unsigned MATRIX_ALIGNMENT = 1;
    real* globalMatrixMem = static_cast<real*>(memoryAllocator.allocateMemory(globalMatrixMemSize * sizeof(real),
                                                                              MATRIX_ALIGNMENT,
                                                                              memkind));
    globalData.address_registry.push_back(globalMatrixMem);

    // copy data defined in source files generated by yateto to the alligned memeory
    real* globalMatrixMemPtr = globalMatrixMem;
    yateto::DeviceCopyManager<real> copy_manager;

    copy_manager.copyFamilyToMemAndSetPtr<init::kDivM>(globalMatrixMemPtr,
                                                       globalData.stiffnessMatrices,
                                                       DEVICE_ALIGNMENT_IN_BYTES);

    copy_manager.copyFamilyToMemAndSetPtr<init::kDivMT>(globalMatrixMemPtr,
                                                        globalData.stiffnessMatricesTransposed,
                                                        DEVICE_ALIGNMENT_IN_BYTES);


    copy_manager.copyFamilyToMemAndSetPtr<init::rDivM>(globalMatrixMemPtr,
                                                       globalData.changeOfBasisMatrices,
                                                       DEVICE_ALIGNMENT_IN_BYTES);

    copy_manager.copyFamilyToMemAndSetPtr<init::rT>(globalMatrixMemPtr,
                                                    globalData.neighbourChangeOfBasisMatricesTransposed,
                                                    DEVICE_ALIGNMENT_IN_BYTES);

    copy_manager.copyFamilyToMemAndSetPtr<init::fMrT>(globalMatrixMemPtr,
                                                      globalData.localChangeOfBasisMatricesTransposed,
                                                      DEVICE_ALIGNMENT_IN_BYTES);

    copy_manager.copyFamilyToMemAndSetPtr<init::fP>(globalMatrixMemPtr,
                                                    globalData.neighbourFluxMatrices,
                                                    DEVICE_ALIGNMENT_IN_BYTES);

    copy_manager.copyTensorToMemAndSetPtr<init::evalAtQP>(globalMatrixMemPtr,
                                                          globalData.evalAtQPMatrix,
                                                          DEVICE_ALIGNMENT_IN_BYTES);

    copy_manager.copyTensorToMemAndSetPtr<init::projectQP>(globalMatrixMemPtr,
                                                           globalData.projectQPMatrix,
                                                           DEVICE_ALIGNMENT_IN_BYTES);

    assert(globalMatrixMemPtr == globalMatrixMem + globalMatrixMemSize);


    // TODO: implement it on gpu
    // multiply all stifness matricess by -1
    // @TODO Integrate this step into the code generator
    for (unsigned transposedStiffness = 0; transposedStiffness < 3; ++transposedStiffness) {
        device_scale_array(-1.0,
                           init::kDivMT::size(transposedStiffness),
                           const_cast<real*>(globalData.stiffnessMatricesTransposed(transposedStiffness)));
    }
}

void seissol::initializers::compareGlobalData(const GlobalData &host, const GlobalDataOnDevice &device) {
    for (unsigned i = 0; i < 3; ++i) {
        std::string array_name = "kDivM(" + std::to_string(i)+ ")";
        device_compare_with_host_array(host.stiffnessMatrices(i),
                                       device.stiffnessMatrices(i),
                                       seissol::init::kDivM::size(i),
                                       array_name.c_str());
    }


    for (unsigned i = 0; i < 3; ++i) {
        std::string array_name = "kDivMT(" + std::to_string(i)+ ")";
        device_compare_with_host_array(host.stiffnessMatricesTransposed(i),
                                       device.stiffnessMatricesTransposed(i),
                                       seissol::init::kDivMT::size(i),
                                       array_name.c_str());
    }

    for (unsigned i = 0; i < 4; ++i) {
        std::string array_name = "rDivM(" + std::to_string(i)+ ")";
        device_compare_with_host_array(host.changeOfBasisMatrices(i),
                                       device.changeOfBasisMatrices(i),
                                       seissol::init::rDivM::size(i),
                                       array_name.c_str());
    }

    for (unsigned i = 0; i < 4; ++i) {
        std::string array_name = "rT(" + std::to_string(i)+ ")";
        device_compare_with_host_array(host.neighbourChangeOfBasisMatricesTransposed(i),
                                       device.neighbourChangeOfBasisMatricesTransposed(i),
                                       seissol::init::rT::size(i),
                                       array_name.c_str());
    }

    for (unsigned i = 0; i < 4; ++i) {
        std::string array_name = "fMrT(" + std::to_string(i)+ ")";
        device_compare_with_host_array(host.localChangeOfBasisMatricesTransposed(i),
                                       device.localChangeOfBasisMatricesTransposed(i),
                                       seissol::init::fMrT::size(i),
                                       array_name.c_str());
    }

    for (unsigned i = 0; i < 3; ++i) {
        std::string array_name = "fP(" + std::to_string(i)+ ")";
        device_compare_with_host_array(host.neighbourFluxMatrices(i),
                                       device.neighbourFluxMatrices(i),
                                       seissol::init::fP::size(i),
                                       array_name.c_str());
    }

    device_compare_with_host_array(host.evalAtQPMatrix,
                                   device.evalAtQPMatrix,
                                   seissol::init::evalAtQP::size(),
                                   "evalAtQP");

    device_compare_with_host_array(host.projectQPMatrix,
                                   device.projectQPMatrix,
                                   seissol::init::projectQP::size(),
                                   "projectQP");
}

/*

// TODO: comment this function
void seissol::initializers::prepareDeviceData(seissol::initializers::LTSTree &tree,
                                              memory::ManagedAllocator& memoryAllocator) {

    unsigned max_layer_cells = 0;
    // iterate though all layer of a tree
    for (LTSTree::leaf_iterator it = tree.beginLeaf(); it != tree.endLeaf(); ++it) {
        // compute the number of elements in the biggest cluster
        if (it->getNumberOfCells() > max_layer_cells)
            max_layer_cells = it->getNumberOfCells();
    }

    // init a linspace on the device
    // NOTE: the linspace is used by temporary variables within the generated code
    seissol::tensor::device_linspace = (unsigned*)(memoryAllocator.allocateMemory(max_layer_cells * sizeof(unsigned),
                                                                                  1,
                                                                                  seissol::memory::Memkind::DeviceGlobalMemory));
    device_init_linspace(max_layer_cells, seissol::tensor::device_linspace);

    seissol::tensor::device_ones = (unsigned*)(memoryAllocator.allocateMemory(max_layer_cells * sizeof(unsigned),
                                                                              1,
                                                                              seissol::memory::Memkind::DeviceGlobalMemory));
    device_init_array(1, max_layer_cells, seissol::tensor::device_ones);

    seissol::tensor::device_zeros = (unsigned*)(memoryAllocator.allocateMemory(max_layer_cells * sizeof(unsigned),
                                                                               1,
                                                                               seissol::memory::Memkind::DeviceGlobalMemory));
    device_init_array(0, max_layer_cells, seissol::tensor::device_zeros);

}
*/
#endif