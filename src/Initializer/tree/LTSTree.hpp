/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2016, SeisSol Group
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
 * Tree for managing lts data.
 **/
 
#ifndef INITIALIZER_TREE_LTSTREE_HPP_
#define INITIALIZER_TREE_LTSTREE_HPP_

#include "LTSInternalNode.hpp"
#include "TimeCluster.hpp"

#include <Initializer/MemoryAllocator.h>

namespace seissol {
  namespace initializers {
    class LTSTree;
  }
}

class seissol::initializers::LTSTree : public seissol::initializers::LTSInternalNode {
private:
  void** m_vars;
  void** m_buckets;
  std::vector<MemoryInfo> varInfo;
  std::vector<MemoryInfo> bucketInfo;
  seissol::memory::ManagedAllocator m_allocator;

  std::vector<size_t> m_variableSizes{};  /*!< sizes of variables within the entire tree in bytes */
  std::vector<size_t> m_bucketSizes{};    /*!< sizes of buckets within the entire tree in bytes */

#ifdef ACL_DEVICE
  std::vector<MemoryInfo> scratchPadInfo{};
  std::vector<size_t> scratchPadSizes{};  /*!< sizes of variables within the entire tree in bytes */
  void** m_scratchPads;
  std::vector<int> m_StreamIds{};
#endif  // ACL_DEVICE

public:
  LTSTree() : m_vars(NULL), m_buckets(NULL) {}
  
  ~LTSTree() { delete[] m_vars; delete[] m_buckets; }
  
  void setNumberOfTimeClusters(unsigned numberOfTimeCluster) {
    setChildren<TimeCluster>(numberOfTimeCluster);
  }
  
  void fixate() {
    setPostOrderPointers();
    for (LTSTree::leaf_iterator it = beginLeaf(); it != endLeaf(); ++it) {
      it->allocatePointerArrays(varInfo.size(), bucketInfo.size());

#ifdef ACL_DEVICE
      it->allocateScratchPadArrays(scratchPadInfo.size());
#endif
    }
  }
  
  inline TimeCluster& child(unsigned index) {
    return *static_cast<TimeCluster*>(m_children[index]);
  }
  
  inline TimeCluster const& child(unsigned index) const {
    return *static_cast<TimeCluster*>(m_children[index]);
  }

  template<typename T>
  T* var(Variable<T> const& handle) {
    assert(handle.index != std::numeric_limits<unsigned>::max());
    assert(m_vars != NULL/* && m_vars[handle.index] != NULL*/);
    return static_cast<T*>(m_vars[handle.index]);
  }

  MemoryInfo const& info(unsigned index) const {
    return varInfo[index];
  }
  
  inline unsigned getNumberOfVariables() const {
    return varInfo.size();
  }
  
  template<typename T>
  void addVar(Variable<T>& handle, LayerMask mask, size_t alignment, seissol::memory::Memkind memkind) {
    handle.index = varInfo.size();
    handle.mask = mask;
    MemoryInfo m;
    m.bytes = sizeof(T)*handle.count;
    m.alignment = alignment;
    m.mask = mask;
    m.memkind = memkind;
    varInfo.push_back(m);
  }
  
  void addBucket(Bucket& handle, size_t alignment, seissol::memory::Memkind memkind) {
    handle.index = bucketInfo.size();
    MemoryInfo m;
    m.alignment = alignment;
    m.memkind = memkind;
    bucketInfo.push_back(m);
  }

#ifdef ACL_DEVICE
  void addScratchPad(ScratchPadMemory& handle, size_t alignment, seissol::memory::Memkind memkind) {
      handle.index = scratchPadInfo.size();
      MemoryInfo m;
      m.alignment = alignment;
      m.memkind = memkind;
      scratchPadInfo.push_back(m);
  }
#endif
  
  void allocateVariables() {
    m_vars = new void*[varInfo.size()];
    m_variableSizes.resize(varInfo.size(), 0);

    for (LTSTree::leaf_iterator it = beginLeaf(); it != endLeaf(); ++it) {
      it->addVariableSizes(varInfo, m_variableSizes);
    }

    for (unsigned var = 0; var < varInfo.size(); ++var) {
      m_vars[var] = m_allocator.allocateMemory(m_variableSizes[var], varInfo[var].alignment, varInfo[var].memkind);
    }

    std::vector<size_t> offsets(varInfo.size(), 0);
    for (LTSTree::leaf_iterator it = beginLeaf(); it != endLeaf(); ++it) {
      it->setMemoryRegionsForVariables(varInfo, m_vars, offsets);
      it->addVariableSizes(varInfo, offsets);
    }
  }

  void allocateBuckets() {
    m_buckets = new void*[bucketInfo.size()];
    m_bucketSizes.resize(bucketInfo.size(), 0);
    
    for (LTSTree::leaf_iterator it = beginLeaf(); it != endLeaf(); ++it) {
      it->addBucketSizes(m_bucketSizes);
    }
    
    for (unsigned bucket = 0; bucket < bucketInfo.size(); ++bucket) {
      m_buckets[bucket] = m_allocator.allocateMemory(m_bucketSizes[bucket], bucketInfo[bucket].alignment, bucketInfo[bucket].memkind);
    }

    std::vector<size_t> offsets(bucketInfo.size(), 0);
    for (LTSTree::leaf_iterator it = beginLeaf(); it != endLeaf(); ++it) {
      it->setMemoryRegionsForBuckets(m_buckets, offsets);
      it->addBucketSizes(offsets);
    }
  }

#ifdef ACL_DEVICE
  // TODO: document (ravil)
  void allocateScratchPads() {
    m_scratchPads = new void*[scratchPadInfo.size()];
    scratchPadSizes.resize(scratchPadInfo.size(), 0);

    for (LTSTree::leaf_iterator it = beginLeaf(); it != endLeaf(); ++it) {
      it->findMaxScratchPadSizes(scratchPadSizes);
    }

    for (unsigned id = 0; id < scratchPadSizes.size(); ++id) {
      assert((scratchPadSizes[id] > 0) && "ERROR: scratch mem. size is equal to zero");
      m_scratchPads[id] = m_allocator.allocateMemory(scratchPadSizes[id],
                                                     scratchPadInfo[id].alignment,
                                                     scratchPadInfo[id].memkind);
    }

    for (LTSTree::leaf_iterator it = beginLeaf(); it != endLeaf(); ++it) {
      it->setMemoryRegionsForScratchPads(m_scratchPads, scratchPadSizes.size());
    }
  }
#endif
  
  void touchVariables() {
    for (LTSTree::leaf_iterator it = beginLeaf(); it != endLeaf(); ++it) {
      it->touchVariables(varInfo);
    }
  }

  const std::vector<size_t>& getVariableSizes() {
    return m_variableSizes;
  }

  const std::vector<size_t>& getBucketSizes() {
    return m_bucketSizes;
  }

#ifdef ACL_DEVICE
    /** Frees variables allocated on a device(s).
     *
     * NOTE: a tree is initialized with a static object. Deallocation of resources is going to happen
     * at the very end of the program. At that time, any device drive is going to be detached. Deallocation under this
     * condition will lead to a segmentation fault. Thus, the user must explicitly deallocate all variables allocated
     * on a device(s) in advance
     * */
    void freeDeviceVariablesExplicitly() {

        // iterate through all variables and find those which were allocated on a device(s)
        for (unsigned i = 0; i < varInfo.size(); ++i) {
            if ((varInfo[i].memkind == seissol::memory::DeviceGlobalMemory)
                || (varInfo[i].memkind == seissol::memory::DeviceUnifiedMemory)) {

                m_allocator.deallocateMemory(m_vars[i], varInfo[i].memkind);
                m_vars[i] = nullptr;
            }
        }

        // iterate through all buckets and find those which were allocated on a device(s)
        for (unsigned i = 0; i < bucketInfo.size(); ++i) {
            if ((bucketInfo[i].memkind == seissol::memory::DeviceGlobalMemory)
                || (bucketInfo[i].memkind == seissol::memory::DeviceUnifiedMemory)) {

                m_allocator.deallocateMemory(m_buckets[i], bucketInfo[i].memkind);
                m_buckets[i] = nullptr;
            }
        }

      for (unsigned i = 0; i < scratchPadInfo.size(); ++i) {
        if ((scratchPadInfo[i].memkind == seissol::memory::DeviceGlobalMemory)
            || (scratchPadInfo[i].memkind == seissol::memory::DeviceUnifiedMemory)) {

          m_allocator.deallocateMemory(m_scratchPads[i], scratchPadInfo[i].memkind);
          m_scratchPads[i] = nullptr;
        }
      }
    }

    void freeLeavesContainersExplicitly() {
        for (LTSTree::leaf_iterator it = beginLeaf(); it != endLeaf(); ++it) {
          it->getLayerContainer().freeConditionalTable();
        }
    }
#endif
};

#endif
