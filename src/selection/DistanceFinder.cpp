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

#ifdef IS_MPI
#include <mpi.h>
#endif

#include "selection/DistanceFinder.hpp"
#include "primitives/Molecule.hpp"

namespace OpenMD {
  
  DistanceFinder::DistanceFinder(SimInfo* info) : info_(info) {
    nObjects_.push_back(info_->getNGlobalAtoms()+info_->getNGlobalRigidBodies());
    nObjects_.push_back(info_->getNGlobalBonds());
    nObjects_.push_back(info_->getNGlobalBends());
    nObjects_.push_back(info_->getNGlobalTorsions());
    nObjects_.push_back(info_->getNGlobalInversions());
    nObjects_.push_back(info_->getNGlobalMolecules());

    stuntdoubles_.resize(nObjects_[STUNTDOUBLE]);
    bonds_.resize(nObjects_[BOND]);
    bends_.resize(nObjects_[BEND]);
    torsions_.resize(nObjects_[TORSION]);
    inversions_.resize(nObjects_[INVERSION]);
    molecules_.resize(nObjects_[MOLECULE]);
    
    SimInfo::MoleculeIterator mi;
    Molecule::AtomIterator ai;
    Molecule::RigidBodyIterator rbIter;
    Molecule::BondIterator bondIter;
    Molecule::BendIterator bendIter;
    Molecule::TorsionIterator torsionIter;
    Molecule::InversionIterator inversionIter;

    Molecule* mol;
    Atom* atom;
    RigidBody* rb;
    Bond* bond;
    Bend* bend;
    Torsion* torsion;
    Inversion* inversion;    
    
    for (mol = info_->beginMolecule(mi); mol != NULL; 
         mol = info_->nextMolecule(mi)) {

      molecules_[mol->getGlobalIndex()] = mol;
      
      for(atom = mol->beginAtom(ai); atom != NULL; 
          atom = mol->nextAtom(ai)) {
	stuntdoubles_[atom->getGlobalIndex()] = atom;
      }
      for (rb = mol->beginRigidBody(rbIter); rb != NULL; 
           rb = mol->nextRigidBody(rbIter)) {
	stuntdoubles_[rb->getGlobalIndex()] = rb;
      }
      for (bond = mol->beginBond(bondIter); bond != NULL; 
           bond = mol->nextBond(bondIter)) {
        bonds_[bond->getGlobalIndex()] = bond;
      }   
      for (bend = mol->beginBend(bendIter); bend != NULL; 
           bend = mol->nextBend(bendIter)) {
        bends_[bend->getGlobalIndex()] = bend;
      }   
      for (torsion = mol->beginTorsion(torsionIter); torsion != NULL; 
           torsion = mol->nextTorsion(torsionIter)) {
        torsions_[torsion->getGlobalIndex()] = torsion;
      }   
      for (inversion = mol->beginInversion(inversionIter); inversion != NULL; 
           inversion = mol->nextInversion(inversionIter)) {
        inversions_[inversion->getGlobalIndex()] = inversion;
      }   

    }
  }

  SelectionSet DistanceFinder::find(const SelectionSet& bs, RealType distance) {
    StuntDouble* center;
    Vector3d centerPos;
    Snapshot* currSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
    SelectionSet bsResult(nObjects_);   
    assert(bsResult.size() == bs.size());
    
#ifdef IS_MPI
    int mol;
    int proc;
    RealType data[3];
    int worldRank;
    MPI_Comm_rank( MPI_COMM_WORLD, &worldRank);
#endif

    
    for (unsigned int j = 0; j < stuntdoubles_.size(); ++j) {
      if (stuntdoubles_[j] != NULL) {
        if (stuntdoubles_[j]->isRigidBody()) {
          RigidBody* rb = static_cast<RigidBody*>(stuntdoubles_[j]);
          rb->updateAtoms();
        }
      }
    }
            
    SelectionSet bsTemp = bs;
    bsTemp = bsTemp.parallelReduce();

    for(int i = 0; i < bsTemp.bitsets_[STUNTDOUBLE].size(); ++i) {
      if (bsTemp.bitsets_[STUNTDOUBLE][i]) {
             
#ifdef IS_MPI
        
        // Now, if we own stuntdouble i, we can use the position, but in
        // parallel, we'll need to let everyone else know what that
        // position is!
        
        mol = info_->getGlobalMolMembership(i);
        proc = info_->getMolToProc(mol);        
        
        if (proc == worldRank) {
          center = stuntdoubles_[i];
          centerPos = center->getPos();
          data[0] = centerPos.x();
          data[1] = centerPos.y();
          data[2] = centerPos.z();
          MPI_Bcast(data, 3, MPI_REALTYPE, proc, MPI_COMM_WORLD);
        } else {
          MPI_Bcast(data, 3, MPI_REALTYPE, proc, MPI_COMM_WORLD);
          centerPos = Vector3d(data);
        }
#else
        center = stuntdoubles_[i];
        centerPos = center->getPos();
#endif               
        for (unsigned int j = 0; j < molecules_.size(); ++j) {
          if (molecules_[j] != NULL) {
            Vector3d r =centerPos - molecules_[j]->getCom();
            currSnapshot->wrapVector(r);
            if (r.length() <= distance) {
              bsResult.bitsets_[MOLECULE].setBitOn(j);
            }
          }
        }
        for (unsigned int j = 0; j < stuntdoubles_.size(); ++j) {
          if (stuntdoubles_[j] != NULL) {
            Vector3d r =centerPos - stuntdoubles_[j]->getPos();
            currSnapshot->wrapVector(r);
            if (r.length() <= distance) {
              bsResult.bitsets_[STUNTDOUBLE].setBitOn(j);
            }
          }
        }
        for (unsigned int j = 0; j < bonds_.size(); ++j) {
          if (bonds_[j] != NULL) {
            Vector3d loc = bonds_[j]->getAtomA()->getPos();
            loc += bonds_[j]->getAtomB()->getPos();
            loc = loc / 2.0;
            Vector3d r = centerPos - loc;
            currSnapshot->wrapVector(r);
            if (r.length() <= distance) {
              bsResult.bitsets_[BOND].setBitOn(j);
            }
          }
        }
        for (unsigned int j = 0; j < bends_.size(); ++j) {
          if (bends_[j] != NULL) {          
            Vector3d loc = bends_[j]->getAtomA()->getPos();
            loc += bends_[j]->getAtomB()->getPos();
            loc += bends_[j]->getAtomC()->getPos();
            loc = loc / 3.0;
            Vector3d r = centerPos - loc;
            currSnapshot->wrapVector(r);
            if (r.length() <= distance) {
              bsResult.bitsets_[BEND].setBitOn(j);
            }
          }
        }
        for (unsigned int j = 0; j < torsions_.size(); ++j) {
          if (torsions_[j] != NULL) {
            Vector3d loc = torsions_[j]->getAtomA()->getPos();
            loc += torsions_[j]->getAtomB()->getPos();
            loc += torsions_[j]->getAtomC()->getPos();
            loc += torsions_[j]->getAtomD()->getPos();
            loc = loc / 4.0;
            Vector3d r = centerPos - loc;
            currSnapshot->wrapVector(r);
            if (r.length() <= distance) {
              bsResult.bitsets_[TORSION].setBitOn(j);
            }
          }
        }
        for (unsigned int j = 0; j < inversions_.size(); ++j) {
          if (inversions_[j] != NULL) {
            Vector3d loc = inversions_[j]->getAtomA()->getPos();
            loc += inversions_[j]->getAtomB()->getPos();
            loc += inversions_[j]->getAtomC()->getPos();
            loc += inversions_[j]->getAtomD()->getPos();
            loc = loc / 4.0;
            Vector3d r = centerPos - loc;
            currSnapshot->wrapVector(r);
            if (r.length() <= distance) {
              bsResult.bitsets_[INVERSION].setBitOn(j);
            }
          }
        }
      }
    }
    return bsResult;
  }
  
  
  SelectionSet DistanceFinder::find(const SelectionSet& bs, RealType distance, int frame ) {
    StuntDouble * center;
    Vector3d centerPos;
    Snapshot* currSnapshot = info_->getSnapshotManager()->getSnapshot(frame);
    SelectionSet bsResult(nObjects_);   
    assert(bsResult.size() == bs.size());

#ifdef IS_MPI
    int mol;
    int proc;
    RealType data[3];
    int worldRank;
    MPI_Comm_rank( MPI_COMM_WORLD, &worldRank);
#endif
 
    for (unsigned int j = 0; j < stuntdoubles_.size(); ++j) {
      if (stuntdoubles_[j] != NULL) {        
        if (stuntdoubles_[j]->isRigidBody()) {
          RigidBody* rb = static_cast<RigidBody*>(stuntdoubles_[j]);
          rb->updateAtoms(frame);
        }
      }
    }
    
    SelectionSet bsTemp = bs;
    bsTemp = bsTemp.parallelReduce();
    
    for(int i = 0; i < bsTemp.bitsets_[STUNTDOUBLE].size(); ++i) {
      if (bsTemp.bitsets_[STUNTDOUBLE][i]) {
        
        // Now, if we own stuntdouble i, we can use the position, but in
        // parallel, we'll need to let everyone else know what that
        // position is!
        
#ifdef IS_MPI
        mol = info_->getGlobalMolMembership(i);
        proc = info_->getMolToProc(mol);
        
        if (proc == worldRank) {
          center = stuntdoubles_[i];
          centerPos = center->getPos(frame);
          data[0] = centerPos.x();
          data[1] = centerPos.y();
          data[2] = centerPos.z();          
          MPI_Bcast(data, 3, MPI_REALTYPE, proc, MPI_COMM_WORLD);
        } else {
          MPI_Bcast(data, 3, MPI_REALTYPE, proc, MPI_COMM_WORLD);
          centerPos = Vector3d(data);
        }
#else
        center = stuntdoubles_[i];
        centerPos = center->getPos(frame);
#endif
        for (unsigned int j = 0; j < molecules_.size(); ++j) {
          if (molecules_[j] != NULL) {
            Vector3d r =centerPos - molecules_[j]->getCom(frame);
            currSnapshot->wrapVector(r);
            if (r.length() <= distance) {
              bsResult.bitsets_[MOLECULE].setBitOn(j);
            }
          }
        }
        
        for (unsigned int j = 0; j < stuntdoubles_.size(); ++j) {
          if (stuntdoubles_[j] != NULL) {          
            Vector3d r =centerPos - stuntdoubles_[j]->getPos(frame);
            currSnapshot->wrapVector(r);
            if (r.length() <= distance) {
              bsResult.bitsets_[STUNTDOUBLE].setBitOn(j);
            }
          }
        }
        for (unsigned int j = 0; j < bonds_.size(); ++j) {
          if (bonds_[j] != NULL) {
            Vector3d loc = bonds_[j]->getAtomA()->getPos(frame);
            loc += bonds_[j]->getAtomB()->getPos(frame);
            loc = loc / 2.0;
            Vector3d r = centerPos - loc;
            currSnapshot->wrapVector(r);
            if (r.length() <= distance) {
              bsResult.bitsets_[BOND].setBitOn(j);
            }
          }
        }
        for (unsigned int j = 0; j < bends_.size(); ++j) {
          if (bends_[j] != NULL) {
            Vector3d loc = bends_[j]->getAtomA()->getPos(frame);
            loc += bends_[j]->getAtomB()->getPos(frame);
            loc += bends_[j]->getAtomC()->getPos(frame);
            loc = loc / 3.0;
            Vector3d r = centerPos - loc;
            currSnapshot->wrapVector(r);
            if (r.length() <= distance) {
              bsResult.bitsets_[BEND].setBitOn(j);
            }
          }
        }
        for (unsigned int j = 0; j < torsions_.size(); ++j) {
          if (torsions_[j] != NULL) {
            Vector3d loc = torsions_[j]->getAtomA()->getPos(frame);
            loc += torsions_[j]->getAtomB()->getPos(frame);
            loc += torsions_[j]->getAtomC()->getPos(frame);
            loc += torsions_[j]->getAtomD()->getPos(frame);
            loc = loc / 4.0;
            Vector3d r = centerPos - loc;
            currSnapshot->wrapVector(r);
            if (r.length() <= distance) {
              bsResult.bitsets_[TORSION].setBitOn(j);
            }
          }
        }
        for (unsigned int j = 0; j < inversions_.size(); ++j) {
          if (inversions_[j] != NULL) {
            Vector3d loc = inversions_[j]->getAtomA()->getPos(frame);
            loc += inversions_[j]->getAtomB()->getPos(frame);
            loc += inversions_[j]->getAtomC()->getPos(frame);
            loc += inversions_[j]->getAtomD()->getPos(frame);
            loc = loc / 4.0;
            Vector3d r = centerPos - loc;
            currSnapshot->wrapVector(r);
            if (r.length() <= distance) {
              bsResult.bitsets_[INVERSION].setBitOn(j);
            }
          }
        }      
      }
    }
    return bsResult;    
  }
}
  
