/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * This software is provided "AS IS," without a warranty of any
 * kind. All express or implied conditions, representations and
 * warranties, including any implied warranty of merchantability,
 * fitness for a particular purpose or non-infringement, are hereby
 * excluded.  The University of Notre Dame and its licensors shall not
 * be liable for any damages suffered by licensee as a result of
 * using, modifying or distributing the software or its
 * derivatives. In no event will the University of Notre Dame or its
 * licensors be liable for any lost revenue, profit or data, or for
 * direct, indirect, special, consequential, incidental or punitive
 * damages, however caused and regardless of the theory of liability,
 * arising out of the use of or inability to use software, even if the
 * University of Notre Dame has been advised of the possibility of
 * such damages.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

#include <algorithm>
#include <functional>

#include "brains/BlockSnapshotManager.hpp"
//#include "utils/residentMem.h"
//#include "utils/physmem.h"
#include "utils/Algorithm.hpp"
#include "brains/SimInfo.hpp"
#include "io/DumpReader.hpp"

using namespace std;

namespace OpenMD {
  BlockSnapshotManager::BlockSnapshotManager(SimInfo* info, 
                                             const std::string& filename, 
					     int storageLayout,
                                             long long int memSize,
                                             int blockCapacity) 
    : SnapshotManager(storageLayout),
      blockCapacity_(blockCapacity), memSize_(memSize), 
      activeBlocks_(blockCapacity_, -1), 
      activeRefCount_(blockCapacity_, 0) {
    
    nAtoms_ = info->getNGlobalAtoms();
    nRigidBodies_ = info->getNGlobalRigidBodies();
    nCutoffGroups_ = info->getNCutoffGroups();
    usePBC_ = info->getSimParams()->getUsePeriodicBoundaryConditions();

    // eliminate suspect calls to figure out free memory:
    // RealType physMem = physmem_total();
    // RealType rssMem = residentMem();
    // RealType avaliablePhysMem = physMem - rssMem;
    
    int bytesPerStuntDouble = DataStorage::getBytesPerStuntDouble(storageLayout);
    int bytesPerCutoffGroup = DataStorage::getBytesPerStuntDouble(DataStorage::dslPosition);

    int bytesPerFrameData = Snapshot::getFrameDataSize();
    int bytesPerFrame = (nRigidBodies_ + nAtoms_) * bytesPerStuntDouble 
      + nCutoffGroups_ * bytesPerCutoffGroup 
      + bytesPerFrameData;

    // total number of frames that can fit in memory
    //RealType frameCapacity = avaliablePhysMem / bytesPerFrame;
    RealType frameCapacity = (RealType) memSize_ / (RealType) bytesPerFrame;

    // number of frames in each block given the need to hold multiple blocks 
    // in memory at the same time:
    nSnapshotPerBlock_ = int(frameCapacity) / blockCapacity_;
    if (nSnapshotPerBlock_ <= 0) {
      std::cerr << "not enough memory to hold two configs!" << std::endl;
    }
    reader_ = new DumpReader(info, filename);
    nframes_ = reader_->getNFrames();
    int nblocks = nframes_ / nSnapshotPerBlock_;
    if (nframes_ % int(nSnapshotPerBlock_) != 0) {
      ++nblocks;
    }  
    
    for (int i = 0; i < nblocks; ++i) {
      blocks_.push_back(SnapshotBlock(i*nSnapshotPerBlock_,
                                      (i+1)*nSnapshotPerBlock_));    
    }
    //the last block may not have nSnapshotPerBlock frames, we need
    //to consider this special situation
    blocks_.back().second = nframes_;

    snapshots_.insert(snapshots_.begin(), nframes_,
                      static_cast<Snapshot*>(NULL));   

    std::cout << "-----------------------------------------------------"
              << std::endl;
    std::cout << "BlockSnapshotManager memory report:" << std::endl;
    std::cout << "\n";
    // std::cout << "  Physical Memory available:\t" << (unsigned long)physMem <<  " bytes" <<std::endl;
    //std::cout << "     Resident Memory in use:\t" << (unsigned long)rssMem << " bytes" <<std::endl;
    //std::cout << "Memory available for OpenMD:\t" << (unsigned long)avaliablePhysMem << " bytes" <<std::endl;
    std::cout << "Memory requested for OpenMD:\t" 
              << (unsigned long)memSize_ << " bytes" << std::endl;
    std::cout << "        Bytes per FrameData:\t" 
              << (unsigned long)bytesPerFrameData << std::endl;
    std::cout << "      Bytes per StuntDouble:\t" 
              << (unsigned long)bytesPerStuntDouble << std::endl;
    std::cout << "     Bytes per Cutoff Group:\t" 
              << (unsigned long)bytesPerCutoffGroup << std::endl;
    std::cout << "            Bytes per Frame:\t" 
              << (unsigned long)bytesPerFrame << std::endl;
    std::cout << "             Frame Capacity:\t"
              << (unsigned long)frameCapacity << std::endl;
    std::cout << "       Frames in trajectory:\t" 
              << (unsigned long)nframes_ << std::endl;
    std::cout << "        Snapshots per Block:\t" 
              << (unsigned long)nSnapshotPerBlock_ << std::endl;
    std::cout << "     Total number of Blocks:\t" 
              << (unsigned long)nblocks << std::endl;
    std::cout << "-----------------------------------------------------"
              << std::endl;
    
  }


  BlockSnapshotManager::~BlockSnapshotManager() {
    currentSnapshot_ = NULL;
    previousSnapshot_ = NULL;
    
    delete reader_;

    std::vector<int>::iterator i;
    for (i = activeBlocks_.begin(); i != activeBlocks_.end(); ++i) {
      if (*i != -1) {
	unloadBlock(*i);
      }
    }
  }

  Snapshot* BlockSnapshotManager::getSnapshot(int id) { 
    currentSnapshot_ = snapshots_[id]; 
    return snapshots_[id]; 
  }

  int BlockSnapshotManager::getNActiveBlocks() {
#ifdef __RWSTD   
    int count = 0;
    std::count_if(activeBlocks_.begin(), activeBlocks_.end(),
				  std::bind(std::not_equal_to<int>(), placeholders::_1, -1), count);
    return count;
#else
    return std::count_if(activeBlocks_.begin(), activeBlocks_.end(),
						 std::bind(std::not_equal_to<int>(), placeholders::_1, -1));
#endif
  }



  bool BlockSnapshotManager::loadBlock(int block) {
    std::vector<int>::iterator i = findActiveBlock(block);
    bool loadSuccess(false);
    if (i != activeBlocks_.end()) {
      // If the block is already in memory, just increase the
      // reference count:
      ++activeRefCount_[i - activeBlocks_.begin()];
      loadSuccess = true;
    } else if (getNActiveBlocks() < blockCapacity_){
      // If the number of active blocks is less than the block
      // capacity, just load the block:
      internalLoad(block);
      loadSuccess = true;
    } else if ( hasZeroRefBlock() ) {
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
    
    if (i != activeBlocks_.end()){
      --activeRefCount_[i - activeBlocks_.begin()];
      if (activeRefCount_[i - activeBlocks_.begin()] < 0) {
	//in case, unloadBlock called multiple times
	activeRefCount_[i - activeBlocks_.begin()]  = 0;
      }

      if (activeRefCount_[i-activeBlocks_.begin()] == 0) {
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

  bool BlockSnapshotManager::hasZeroRefBlock(){
    return std::find(activeRefCount_.begin(), activeRefCount_.end(), 0)
      != activeRefCount_.end() ?  true : false;
  }
  
  int BlockSnapshotManager::getFirstZeroRefBlock(){
    std::vector<int>::iterator i = std::find(activeRefCount_.begin(),
                                             activeRefCount_.end(), 0);
    return i != activeRefCount_.end() ?
      activeBlocks_[i - activeRefCount_.begin()] : -1;
  }

  std::vector<int> BlockSnapshotManager::getActiveBlocks() {
    std::vector<int> result;
    OpenMD::copy_if(activeBlocks_.begin(), activeBlocks_.end(),
                    std::back_inserter(result), 
                    std::bind(std::not_equal_to<int>(), placeholders::_1, -1));
    return result;    
  }

  Snapshot* BlockSnapshotManager::loadFrame(int frame){
    Snapshot* snapshot = new Snapshot(nAtoms_, nRigidBodies_, nCutoffGroups_, 
                                      getStorageLayout(), usePBC_);
    snapshot->setID(frame);
    snapshot->clearDerivedProperties();
    
    currentSnapshot_ = snapshot;   
    reader_->readFrame(frame);

    return snapshot;
  }

  int BlockSnapshotManager::getNFrames() {
    return reader_->getNFrames();
  }

  void BlockSnapshotManager::needCOMprops(bool ncp) {
    reader_->setNeedCOMprops(ncp);
  }

}
