/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved.
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
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
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
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the following paper when you publish your work:
 *
 * [1] Drisko et al., J. Open Source Softw. 9, 7004 (2024).
 *
 * Good starting points for code and simulation methodology are:
 *
 * [2] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [3] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [4] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [7] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 * [9] Drisko & Gezelter, J. Chem. Theory Comput. 20, 4986-4997 (2024).
 */

#include "brains/BlockSnapshotManager.hpp"

#include <algorithm>
#include <functional>
#include <iostream>
#include <string>

#include "brains/SimInfo.hpp"
#include "io/DumpReader.hpp"

namespace OpenMD {

  BlockSnapshotManager::BlockSnapshotManager(
      SimInfo* info, const std::string& filename, int atomStorageLayout,
      int rigidBodyStorageLayout, int cutoffGroupStorageLayout,
      long long int memSize, int blockCapacity) :
      SnapshotManager(atomStorageLayout, rigidBodyStorageLayout,
                      cutoffGroupStorageLayout),
      blockCapacity_(blockCapacity), memSize_(memSize),
      activeBlocks_(blockCapacity_, -1), activeRefCount_(blockCapacity_, 0) {
    nAtoms_        = info->getNGlobalAtoms();
    nRigidBodies_  = info->getNGlobalRigidBodies();
    nCutoffGroups_ = info->getNCutoffGroups();
    usePBC_        = info->getSimParams()->getUsePeriodicBoundaryConditions();

    // eliminate suspect calls to figure out free memory:
    // RealType physMem = physmem_total();
    // RealType rssMem = residentMem();
    // RealType avaliablePhysMem = physMem - rssMem;

    int bytesPerAtom = DataStorage::getBytesPerStuntDouble(atomStorageLayout);
    int bytesPerRigidBody =
        DataStorage::getBytesPerStuntDouble(rigidBodyStorageLayout);
    int bytesPerCutoffGroup =
        DataStorage::getBytesPerStuntDouble(DataStorage::dslPosition);

    int bytesPerFrameData = Snapshot::getFrameDataSize();
    int bytesPerFrame =
        nRigidBodies_ * bytesPerRigidBody + nAtoms_ * bytesPerAtom +
        nCutoffGroups_ * bytesPerCutoffGroup + bytesPerFrameData;

    // total number of frames that can fit in memory
    // RealType frameCapacity = avaliablePhysMem / bytesPerFrame;
    RealType frameCapacity = (RealType)memSize_ / (RealType)bytesPerFrame;

    // number of frames in each block given the need to hold multiple blocks
    // in memory at the same time:
    nSnapshotPerBlock_ = int(frameCapacity) / blockCapacity_;
    if (nSnapshotPerBlock_ <= 0) {
      std::cerr << "not enough memory to hold two configs!\n";
    }
    reader_     = new DumpReader(info, filename);
    nframes_    = reader_->getNFrames();
    int nblocks = nframes_ / nSnapshotPerBlock_;
    if (nframes_ % int(nSnapshotPerBlock_) != 0) { ++nblocks; }

    for (int i = 0; i < nblocks; ++i) {
      blocks_.push_back(
          SnapshotBlock(i * nSnapshotPerBlock_, (i + 1) * nSnapshotPerBlock_));
    }
    // the last block may not have nSnapshotPerBlock frames, we need
    // to consider this special situation
    blocks_.back().second = nframes_;

    snapshots_.insert(snapshots_.begin(), nframes_,
                      static_cast<Snapshot*>(NULL));

    std::cout << "-----------------------------------------------------\n";
    std::cout << "BlockSnapshotManager memory report:\n";
    std::cout << "\n";
    // std::cout << "  Physical Memory available:\t" << (unsigned long)physMem
    // << " bytes"
    // <<std::endl;
    // std::cout << "     Resident Memory in use:\t" << (unsigned long)rssMem <<
    // " bytes"
    // <<std::endl; std::cout << "Memory available for OpenMD:\t" << (unsigned
    // long)avaliablePhysMem << " bytes" <<std::endl;
    std::cout << "Memory requested for OpenMD:\t" << (unsigned long)memSize_
              << " bytes" << std::endl;
    std::cout << "        Bytes per FrameData:\t"
              << (unsigned long)bytesPerFrameData << std::endl;
    std::cout << "             Bytes per Atom:\t" << (unsigned long)bytesPerAtom
              << std::endl;
    std::cout << "        Bytes per RigidBody:\t"
              << (unsigned long)bytesPerRigidBody << std::endl;
    std::cout << "     Bytes per Cutoff Group:\t"
              << (unsigned long)bytesPerCutoffGroup << std::endl;
    std::cout << "            Bytes per Frame:\t"
              << (unsigned long)bytesPerFrame << std::endl;
    std::cout << "             Frame Capacity:\t"
              << (unsigned long)frameCapacity << std::endl;
    std::cout << "       Frames in trajectory:\t" << (unsigned long)nframes_
              << std::endl;
    std::cout << "        Snapshots per Block:\t"
              << (unsigned long)nSnapshotPerBlock_ << std::endl;
    std::cout << "     Total number of Blocks:\t" << (unsigned long)nblocks
              << std::endl;
    std::cout << "-----------------------------------------------------"
              << std::endl;
  }

  BlockSnapshotManager::~BlockSnapshotManager() {
    currentSnapshot_  = NULL;
    previousSnapshot_ = NULL;

    delete reader_;

    std::vector<int>::iterator i;
    for (i = activeBlocks_.begin(); i != activeBlocks_.end(); ++i) {
      if (*i != -1) { unloadBlock(*i); }
    }
  }

  Snapshot* BlockSnapshotManager::getSnapshot(int id) {
    currentSnapshot_ = snapshots_[id];
    return snapshots_[id];
  }

  int BlockSnapshotManager::getNActiveBlocks() {
    return std::count_if(activeBlocks_.begin(), activeBlocks_.end(),
                         [](int val) { return val != -1; });
  }

  bool BlockSnapshotManager::loadBlock(int block) {
    std::vector<int>::iterator i = findActiveBlock(block);
    bool loadSuccess(false);
    if (i != activeBlocks_.end()) {
      // If the block is already in memory, just increase the
      // reference count:
      ++activeRefCount_[i - activeBlocks_.begin()];
      loadSuccess = true;
    } else if (getNActiveBlocks() < blockCapacity_) {
      // If the number of active blocks is less than the block
      // capacity, just load the block:
      internalLoad(block);
      loadSuccess = true;
    } else if (hasZeroRefBlock()) {
      // If we have already reached the block capacity, we need to
      // unload a block with 0 references:
      int zeroRefBlock = getFirstZeroRefBlock();
      assert(zeroRefBlock != -1);
      internalUnload(zeroRefBlock);
      internalLoad(block);
    } else {
      // We have reached capacity and all blocks in memory are have
      // non-zero references:
      loadSuccess = false;
    }
    return loadSuccess;
  }

  bool BlockSnapshotManager::unloadBlock(int block) {
    bool unloadSuccess;
    std::vector<int>::iterator i = findActiveBlock(block);

    if (i != activeBlocks_.end()) {
      --activeRefCount_[i - activeBlocks_.begin()];
      if (activeRefCount_[i - activeBlocks_.begin()] < 0) {
        // in case, unloadBlock called multiple times
        activeRefCount_[i - activeBlocks_.begin()] = 0;
      }

      if (activeRefCount_[i - activeBlocks_.begin()] == 0) {
        internalUnload(block);
      }

      unloadSuccess = true;
    } else {
      unloadSuccess = false;
    }

    return unloadSuccess;
  }

  void BlockSnapshotManager::internalLoad(int block) {
    for (int i = blocks_[block].first; i < blocks_[block].second; ++i) {
      snapshots_[i] = loadFrame(i);
    }

    std::vector<int>::iterator j;
    j = std::find(activeBlocks_.begin(), activeBlocks_.end(), -1);
    assert(j != activeBlocks_.end());
    *j = block;
    ++activeRefCount_[j - activeBlocks_.begin()];
  }

  void BlockSnapshotManager::internalUnload(int block) {
    for (int i = blocks_[block].first; i < blocks_[block].second; ++i) {
      delete snapshots_[i];
      snapshots_[i] = NULL;
    }
    std::vector<int>::iterator j;
    j = std::find(activeBlocks_.begin(), activeBlocks_.end(), block);
    assert(j != activeBlocks_.end());
    *j = -1;
  }

  bool BlockSnapshotManager::hasZeroRefBlock() {
    return std::find(activeRefCount_.begin(), activeRefCount_.end(), 0) !=
                   activeRefCount_.end() ?
               true :
               false;
  }

  int BlockSnapshotManager::getFirstZeroRefBlock() {
    std::vector<int>::iterator i =
        std::find(activeRefCount_.begin(), activeRefCount_.end(), 0);
    return i != activeRefCount_.end() ?
               activeBlocks_[i - activeRefCount_.begin()] :
               -1;
  }

  std::vector<int> BlockSnapshotManager::getActiveBlocks() {
    std::vector<int> result;
    std::copy_if(
        activeBlocks_.begin(), activeBlocks_.end(), std::back_inserter(result),
        std::bind(std::not_equal_to<int>(), std::placeholders::_1, -1));
    return result;
  }

  Snapshot* BlockSnapshotManager::loadFrame(int frame) {
    Snapshot* snapshot = new Snapshot(
        nAtoms_, nRigidBodies_, nCutoffGroups_, getAtomStorageLayout(),
        getRigidBodyStorageLayout(), getCutoffGroupStorageLayout(), usePBC_);
    snapshot->setID(frame);
    snapshot->clearDerivedProperties();

    currentSnapshot_ = snapshot;
    reader_->readFrame(frame);

    return snapshot;
  }

  int BlockSnapshotManager::getNFrames() { return reader_->getNFrames(); }

  void BlockSnapshotManager::needCOMprops(bool ncp) {
    reader_->setNeedCOMprops(ncp);
  }

}  // namespace OpenMD
