/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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
 */
#include <algorithm>
#include "brains/BlockSnapshotManager.hpp"
#include "utils/physmem.h"
#include "brains/SimInfo.hpp"
#include "io/DumpReader.hpp"

namespace oopse {
BlockSnapshotMananger::BlockSnapshotMananger(SimInfo* info, const std::string& filename, 
    int storageLayout, int blockCapacity) 
    : SnapshotManager(storageLayout), info_(info), blockCapacity_(blockCapacity), activeBlocks_(blockCapacity_, -1) {

    nAtoms_ = info->getNGlobalAtoms();
    nRigidBodies_ = info->getNGlobalRigidBodies();

    double avalPhysMem = physmem_available();
    
    int bytesPerStuntDouble = DataStorage::getBytesPerStuntDouble(storageLayout);

    int bytesPerFrame = (nRigidBodies_ + nAtoms_) * bytesPerStuntDouble;

    int frameCapacity = int (avalPhysMem / bytesPerFrame);
    
    nSnapshotPerBlock_ = frameCapacity /blockCapacity_ ;

    reader_ = new DumpReader(info, filename);
    nframes_ = reader_->getNFrames();

    int nblocks = nframes_ / nSnapshotPerBlock_;
    if (nframes_ % nSnapshotPerBlock_ != 0) {
        ++nblocks;
    }  
    
    for (int i = 0; i < nblocks; ++i) {
        blocks_.push_back(SnapshotBlock(i, (i+1)*nSnapshotPerBlock_);    
    }
    //the last block may not have nSnapshotPerBlock frames, we need to consider this special situation
    blocks_.back.second = nframes_;

    snapshots_.insert(snapshots_.begin(), nframes_, NULL);   
    
}


BlockSnapshotMananger::~BlockSnapshotMananger() {
    currentSnapshot_ = NULL;
    previousSnapshot_ = NULL;
    
    delete reader_;
    std::for_each(activeBlocks_.begin(), activeBlocks_.end(), unloadBlock);
}

int BlockSnapshotManager::getNActiveBlocks() {
    return std::count_if(activeBlocks_.begin(), activeBlocks_.end(), std::not);
}

bool BlockSnapshotManager::isBlockActive(int block) {
    return std::find(activeBlocks_.begin(), activeBlocks_.end(), block) != activeBlocks_.end() ? true : false;
}

bool BlockSnapshotManager::loadBlock(int block) {
    bool loadSuccess;
    if (isBlockActive(block)) {
        loadSuccess = true;
    } else if (getNActiveBlocks() < blockCapacity_){

        for (int i = blocks_[block].first; i < blocks_[block].second; ++i) {
            snapshots_[i] = loadFrame(i);
        }
        
        std::vector<int>::iterator j;
        j = std::find(activeBlocks_.begin(), activeBlocks_.end(), -1);
        assert(j != activeBlocks_.end());
        *j = block;    
        loadSuccess = true;
    }else {
        loadSuccess = false;
    }

    return loadSuccess;
}

bool BlockSnapshotManager::unloadBlock(int block) {
    bool unloadSuccess;
    if (!isBlockActive(block)){
        unloadSuccess = false;
    } else {
        for (int i = blocks_[block].first; i < blocks_[block].second; ++i) {
            delete snapshots_[i];
            snapshots_[i] = NULL;
        }
        std::vector<int>::iterator j;
        j = std::find(activeBlocks_.begin(), activeBlocks_.end(), block);
        assert(j != activeBlocks_.end());
        *j = -1;
    }
}

std::vector<int> BlockSnapshotManager::getActiveBlocks() {
    std::vector<int> result;
    std::copy_if(activeBlocks_.begin(), activeBlocks_.end());
}

Snapshot* BlockSnapshotManager::loadFrame(int frame){
    Snapshot* snapshot = new Snapshot(nAtoms_, nRigidBodies_, getStorageLayout());
    snapshot->setID(frame);
    setCurrentSnapshot(snapshot);   /** @todo fixed me */
    reader_->readFrame(frame);
    return snapshot;
}

int BlockSnapshotManager::getNFrames() {
    return reader_->getNFrames();
}

}
