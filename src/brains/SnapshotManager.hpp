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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
  
/**
 * @file SnapshotManager.hpp
 * @author tlin
 * @date 10/20/2004
 * @time 23:56am
 * @version 1.0
 */
#ifndef BRAINS_SNAPSHOTMANAGER_HPP
#define BRAINS_SNAPSHOTMANAGER_HPP

#include "brains/Snapshot.hpp"

namespace OpenMD{

  /**
   * @class SnapshotManager SnapshotManager.hpp "brains/SnapshotManager.hpp"
   * @brief SnapshotManager class is an abstract class which maintains 
   * a series of snapshots.
   * 
   * @see SimSnapshotManager
   * @see PropSnapshotManager
   */
  class SnapshotManager {
  public:

    virtual ~SnapshotManager() {
      delete currentSnapshot_;
      delete previousSnapshot_;
    }
            
    virtual bool advance() { return true; }

    virtual Snapshot* getSnapshot(int id) = 0;

    /**
     * Returns the pointer of previous snapshot
     * @return the pointer of previous snapshot
     */
    Snapshot* getPrevSnapshot() {
      return previousSnapshot_;
    }

    /**
     * Returns the pointer of current snapshot
     * @return the pointer of current snapshot
     */            
    Snapshot* getCurrentSnapshot() {
      return currentSnapshot_;
    }

    int getStorageLayout() {
      return storageLayout_;
    }

  private:
    int storageLayout_;

  protected:

    SnapshotManager(int storageLayout) : storageLayout_(storageLayout), currentSnapshot_(NULL), previousSnapshot_(NULL) {
    }
            
    Snapshot* currentSnapshot_;
    Snapshot* previousSnapshot_;
            

  };

}
#endif //BRAINS_SNAPSHOTMANAGER_HPP

