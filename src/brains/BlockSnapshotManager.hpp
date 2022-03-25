/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
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

#ifndef BRAINS_BLOCKSNAPSHOTMANAGER_HPP
#define BRAINS_BLOCKSNAPSHOTMANAGER_HPP

#include <vector>

#include "brains/SnapshotManager.hpp"

namespace OpenMD {

  class SimInfo;
  class DumpReader;

  using SnapshotBlock = std::pair<int, int>;

  /**
   * @class BlockSnapshotManager
   * @todo document
   */
  class BlockSnapshotManager : public SnapshotManager {
  public:
    BlockSnapshotManager(SimInfo* info, const std::string& filename,
                         int storageLayout, long long int memSize,
                         int blockCapacity = 2);
    ~BlockSnapshotManager();

    virtual Snapshot* getSnapshot(int id);

    /** Returns number of snapshot blocks in this BlockSnapshotManager*/
    int getNBlocks() { return blocks_.size(); }

    SnapshotBlock getSnapshotBlock(int block) { return blocks_.at(block); }

    int getNActiveBlocks();

    void needCOMprops(bool ncp);

    bool isBlockActive(int block) {
      return findActiveBlock(block) != activeBlocks_.end() ? true : false;
    }

    bool loadBlock(int block);

    bool unloadBlock(int block);

    std::vector<int> getActiveBlocks();

    int getBlockCapacity() { return blockCapacity_; }

    int getNFrames();

  private:
    std::vector<int>::iterator findActiveBlock(int block) {
      return std::find(activeBlocks_.begin(), activeBlocks_.end(), block);
    }

    bool hasZeroRefBlock();

    int getFirstZeroRefBlock();

    void internalLoad(int block);
    void internalUnload(int block);
    Snapshot* loadFrame(int frame);

    int blockCapacity_;
    long long int memSize_;

    std::vector<Snapshot*> snapshots_;
    std::vector<SnapshotBlock> blocks_;
    std::vector<int> activeBlocks_;
    std::vector<int> activeRefCount_;

    int nAtoms_;
    int nRigidBodies_;
    int nCutoffGroups_;
    bool usePBC_;

    DumpReader* reader_;
    int nframes_;
    int nSnapshotPerBlock_;
  };
}  // namespace OpenMD

#endif
