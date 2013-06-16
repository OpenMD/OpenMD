/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
#ifndef BRAINS_BLOCKSNAPSHOTMANAGER_HPP
#define BRAINS_BLOCKSNAPSHOTMANAGER_HPP
#include <vector>

#include "brains/SnapshotManager.hpp"
namespace OpenMD {

  class SimInfo;
  class DumpReader;

  typedef std::pair<int, int> SnapshotBlock;

  /**
   * @class BlockSnapshotManager
   * @todo document
   */
  class BlockSnapshotManager : public SnapshotManager{

  public:
    BlockSnapshotManager(SimInfo* info, const std::string& filename, int storageLayout, long long int memSize, int blockCapacity = 2);
    ~BlockSnapshotManager();
        
    virtual Snapshot* getSnapshot(int id);


    /** Returns number of snapshot blocks in this BlockSnapshotManager*/
    int getNBlocks() {
      return blocks_.size();
    }

    SnapshotBlock getSnapshotBlock(int block) {
      return blocks_.at(block);
    }
        
    int getNActiveBlocks();

    void needCOMprops(bool ncp);


    bool isBlockActive(int block) {
      return  findActiveBlock(block) != activeBlocks_.end() ? true : false;
    }        

    bool loadBlock(int block);
        
    bool unloadBlock(int block);

    std::vector<int> getActiveBlocks();

    int getBlockCapacity() {
      return blockCapacity_;                
    }

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
        
    SimInfo* info_;
    int blockCapacity_;
    long long int memSize_;

    std::vector<Snapshot*> snapshots_;
    std::vector<SnapshotBlock> blocks_;        
    std::vector<int> activeBlocks_;
    std::vector<int> activeRefCount_;
        
    int nAtoms_;
    int nRigidBodies_;
    int nCutoffGroups_;
    
    DumpReader* reader_;
    int nframes_;
    int nSnapshotPerBlock_;

  };

}

#endif
