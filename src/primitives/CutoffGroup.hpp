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
 * [4]  Vardeman & Gezelter, in progress (2009).                        
 */
 
#ifndef PRIMITIVES_CUTOFFGROUP_HPP

#define PRIMITIVES_CUTOFFGROUP_HPP

#include "primitives/Atom.hpp"
#include "math/Vector3.hpp"

namespace OpenMD {
  class CutoffGroup {
  public:
    
    CutoffGroup() :  snapshotMan_(NULL) {

      storage_ = &Snapshot::cgData;
      haveTotalMass = false;
      totalMass = 0.0;
    }
    
    /**
     * Sets the Snapshot Manager of this cutoffGroup
     */
    void setSnapshotManager(SnapshotManager* sman) {
      snapshotMan_ = sman;
    }


    void addAtom(Atom *atom) {
      cutoffAtomList.push_back(atom);
    }
    
    Atom *beginAtom(std::vector<Atom *>::iterator & i) {
      i = cutoffAtomList.begin();
      return i != cutoffAtomList.end() ? *i : NULL;
    }
    
    Atom *nextAtom(std::vector<Atom *>::iterator & i) {
      i++;
      return i != cutoffAtomList.end() ? *i : NULL;
    }

    std::vector<Atom*> getAtoms() { return cutoffAtomList; }
    RealType getMass() {
      std::vector<Atom *>::iterator i;
      Atom * atom;
      RealType mass;
      
      if (!haveTotalMass) {
	totalMass = 0;
        
	for(atom = beginAtom(i); atom != NULL; atom = nextAtom(i)) {
	  mass = atom->getMass();
	  totalMass += mass;
	}
        
	haveTotalMass = true;
      }
      
      return totalMass;
    }
    
    void updateCOM() {
      std::vector<Atom *>::iterator i;
      Atom * atom;

      DataStorage&  data = snapshotMan_->getCurrentSnapshot()->*storage_;
     
      totalMass = getMass();
      
      if (cutoffAtomList.size() == 1) {
        std::cerr << "YO!\n";
        std::cerr << "atipos = " << beginAtom(i)->getPos() << "\n";
        std::cerr << "lI = " << localIndex_ << "\n";

        data.position[localIndex_] = beginAtom(i)->getPos();
        std::cerr << "YOYO!\n";
      } else {
        data.position[localIndex_] = V3Zero;
	for(atom = beginAtom(i); atom != NULL; atom = nextAtom(i)) {
          data.position[localIndex_] += atom->getMass() * atom->getPos();
	}        
	data.position[localIndex_] /= totalMass;
      }
    }


    Vector3d getPos() {
      return ((snapshotMan_->getCurrentSnapshot())->*storage_).position[localIndex_];      
    }
    
    int getNumAtom() {
      return cutoffAtomList.size();
    }
    
    int getGlobalIndex() {
      return globalIndex;
    }
    
    void setGlobalIndex(int id) {
      this->globalIndex = id;
    }

    /** 
     * Returns the local index of this cutoffGroup
     * @return the local index of this cutoffGroup
     */
    int getLocalIndex() {
      return localIndex_;
    }

    /**
     * Sets the local index of this cutoffGroup
     * @param index new index to be set
     */        
    void setLocalIndex(int index) {
      localIndex_ = index;
    }
    
  private:
    
    std::vector<Atom *>cutoffAtomList;
    bool haveTotalMass;
    RealType totalMass;
    int globalIndex;    

    int localIndex_;
    DataStoragePointer storage_;
    SnapshotManager* snapshotMan_;

  };  
} //end namespace OpenMD
#endif //PRIMITIVES_CUTOFFGROUP_HPP  
