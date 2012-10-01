/*
 * Copyright (c) 2005, 2010 The University of Notre Dame. All Rights Reserved.
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

#include "selection/DistanceFinder.hpp"
#include "primitives/Molecule.hpp"
namespace OpenMD {

  DistanceFinder::DistanceFinder(SimInfo* info) : info_(info) {

    nStuntDoubles_ = info_->getNGlobalAtoms() + info_->getNGlobalRigidBodies();
    stuntdoubles_.resize(nStuntDoubles_);
    
    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::AtomIterator ai;
    Atom* atom;
    Molecule::RigidBodyIterator rbIter;
    RigidBody* rb;

    
    for (mol = info_->beginMolecule(mi); mol != NULL; 
         mol = info_->nextMolecule(mi)) {
        
      for(atom = mol->beginAtom(ai); atom != NULL; 
          atom = mol->nextAtom(ai)) {
	stuntdoubles_[atom->getGlobalIndex()] = atom;
      }

      for (rb = mol->beginRigidBody(rbIter); rb != NULL; 
           rb = mol->nextRigidBody(rbIter)) {
	stuntdoubles_[rb->getGlobalIndex()] = rb;
      }
    }
  }

  OpenMDBitSet DistanceFinder::find(const OpenMDBitSet& bs, RealType distance) {
    StuntDouble * center;
    Vector3d centerPos;
    Snapshot* currSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
    OpenMDBitSet bsResult(nStuntDoubles_);
    assert(bsResult.size() == bs.size());
   
    for (unsigned int j = 0; j < stuntdoubles_.size(); ++j) {
      if (stuntdoubles_[j]->isRigidBody()) {
        RigidBody* rb = static_cast<RigidBody*>(stuntdoubles_[j]);
        rb->updateAtoms();
      }
    }
   
    // This will fail in parallel because i might not be on this processor.

    for (int i = bs.firstOnBit(); i != -1; i = bs.nextOnBit(i)) {
      center = stuntdoubles_[i];
      centerPos = center->getPos();
      for (unsigned int j = 0; j < stuntdoubles_.size(); ++j) {
	Vector3d r =centerPos - stuntdoubles_[j]->getPos();
	currSnapshot->wrapVector(r);
	if (r.length() <= distance) {
	  bsResult.setBitOn(j);
	}
      }
    }

    return bsResult;
  }

}
