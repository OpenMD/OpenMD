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
 
/**
 * @file SimInfo.cpp
 * @author    tlin
 * @date  11/02/2004
 * @version 1.0
 */

#include <algorithm>
#include <set>
#include <map>

#include "brains/SimInfo.hpp"
#include "math/Vector3.hpp"
#include "primitives/Molecule.hpp"
#include "primitives/StuntDouble.hpp"
#include "UseTheForce/DarkSide/neighborLists_interface.h"
#include "utils/MemoryUtils.hpp"
#include "utils/simError.h"
#include "selection/SelectionManager.hpp"
#include "io/ForceFieldOptions.hpp"
#include "UseTheForce/ForceField.hpp"
#include "nonbonded/SwitchingFunction.hpp"

#ifdef IS_MPI
#include "UseTheForce/mpiComponentPlan.h"
#include "UseTheForce/DarkSide/simParallel_interface.h"
#endif 

using namespace std;
namespace OpenMD {
  
  SimInfo::SimInfo(ForceField* ff, Globals* simParams) : 
    forceField_(ff), simParams_(simParams), 
    ndf_(0), fdf_local(0), ndfRaw_(0), ndfTrans_(0), nZconstraint_(0),
    nGlobalMols_(0), nGlobalAtoms_(0), nGlobalCutoffGroups_(0), 
    nGlobalIntegrableObjects_(0), nGlobalRigidBodies_(0),
    nAtoms_(0), nBonds_(0),  nBends_(0), nTorsions_(0), nInversions_(0), 
    nRigidBodies_(0), nIntegrableObjects_(0), nCutoffGroups_(0), 
    nConstraints_(0), sman_(NULL), fortranInitialized_(false), 
    calcBoxDipole_(false), useAtomicVirial_(true) {    
    
    MoleculeStamp* molStamp;
    int nMolWithSameStamp;
    int nCutoffAtoms = 0; // number of atoms belong to cutoff groups
    int nGroups = 0;       //total cutoff groups defined in meta-data file
    CutoffGroupStamp* cgStamp;    
    RigidBodyStamp* rbStamp;
    int nRigidAtoms = 0;
    
    vector<Component*> components = simParams->getComponents();
    
    for (vector<Component*>::iterator i = components.begin(); i !=components.end(); ++i) {
      molStamp = (*i)->getMoleculeStamp();
      nMolWithSameStamp = (*i)->getNMol();
      
      addMoleculeStamp(molStamp, nMolWithSameStamp);
      
      //calculate atoms in molecules
      nGlobalAtoms_ += molStamp->getNAtoms() *nMolWithSameStamp;   
      
      //calculate atoms in cutoff groups
      int nAtomsInGroups = 0;
      int nCutoffGroupsInStamp = molStamp->getNCutoffGroups();
      
      for (int j=0; j < nCutoffGroupsInStamp; j++) {
        cgStamp = molStamp->getCutoffGroupStamp(j);
        nAtomsInGroups += cgStamp->getNMembers();
      }
      
      nGroups += nCutoffGroupsInStamp * nMolWithSameStamp;
      
      nCutoffAtoms += nAtomsInGroups * nMolWithSameStamp;            
      
      //calculate atoms in rigid bodies
      int nAtomsInRigidBodies = 0;
      int nRigidBodiesInStamp = molStamp->getNRigidBodies();
      
      for (int j=0; j < nRigidBodiesInStamp; j++) {
        rbStamp = molStamp->getRigidBodyStamp(j);
        nAtomsInRigidBodies += rbStamp->getNMembers();
      }
      
      nGlobalRigidBodies_ += nRigidBodiesInStamp * nMolWithSameStamp;
      nRigidAtoms += nAtomsInRigidBodies * nMolWithSameStamp;            
      
    }
    
    //every free atom (atom does not belong to cutoff groups) is a cutoff 
    //group therefore the total number of cutoff groups in the system is 
    //equal to the total number of atoms minus number of atoms belong to 
    //cutoff group defined in meta-data file plus the number of cutoff 
    //groups defined in meta-data file
    nGlobalCutoffGroups_ = nGlobalAtoms_ - nCutoffAtoms + nGroups;
    
    //every free atom (atom does not belong to rigid bodies) is an 
    //integrable object therefore the total number of integrable objects 
    //in the system is equal to the total number of atoms minus number of 
    //atoms belong to rigid body defined in meta-data file plus the number 
    //of rigid bodies defined in meta-data file
    nGlobalIntegrableObjects_ = nGlobalAtoms_ - nRigidAtoms 
      + nGlobalRigidBodies_;
    
    nGlobalMols_ = molStampIds_.size();
    molToProcMap_.resize(nGlobalMols_);
  }
  
  SimInfo::~SimInfo() {
    map<int, Molecule*>::iterator i;
    for (i = molecules_.begin(); i != molecules_.end(); ++i) {
      delete i->second;
    }
    molecules_.clear();
       
    delete sman_;
    delete simParams_;
    delete forceField_;
  }


  bool SimInfo::addMolecule(Molecule* mol) {
    MoleculeIterator i;
    
    i = molecules_.find(mol->getGlobalIndex());
    if (i == molecules_.end() ) {
      
      molecules_.insert(make_pair(mol->getGlobalIndex(), mol));
      
      nAtoms_ += mol->getNAtoms();
      nBonds_ += mol->getNBonds();
      nBends_ += mol->getNBends();
      nTorsions_ += mol->getNTorsions();
      nInversions_ += mol->getNInversions();
      nRigidBodies_ += mol->getNRigidBodies();
      nIntegrableObjects_ += mol->getNIntegrableObjects();
      nCutoffGroups_ += mol->getNCutoffGroups();
      nConstraints_ += mol->getNConstraintPairs();
      
      addInteractionPairs(mol);
      
      return true;
    } else {
      return false;
    }
  }
  
  bool SimInfo::removeMolecule(Molecule* mol) {
    MoleculeIterator i;
    i = molecules_.find(mol->getGlobalIndex());

    if (i != molecules_.end() ) {

      assert(mol == i->second);
        
      nAtoms_ -= mol->getNAtoms();
      nBonds_ -= mol->getNBonds();
      nBends_ -= mol->getNBends();
      nTorsions_ -= mol->getNTorsions();
      nInversions_ -= mol->getNInversions();
      nRigidBodies_ -= mol->getNRigidBodies();
      nIntegrableObjects_ -= mol->getNIntegrableObjects();
      nCutoffGroups_ -= mol->getNCutoffGroups();
      nConstraints_ -= mol->getNConstraintPairs();

      removeInteractionPairs(mol);
      molecules_.erase(mol->getGlobalIndex());

      delete mol;
        
      return true;
    } else {
      return false;
    }
  }    

        
  Molecule* SimInfo::beginMolecule(MoleculeIterator& i) {
    i = molecules_.begin();
    return i == molecules_.end() ? NULL : i->second;
  }    

  Molecule* SimInfo::nextMolecule(MoleculeIterator& i) {
    ++i;
    return i == molecules_.end() ? NULL : i->second;    
  }


  void SimInfo::calcNdf() {
    int ndf_local;
    MoleculeIterator i;
    vector<StuntDouble*>::iterator j;
    Molecule* mol;
    StuntDouble* integrableObject;

    ndf_local = 0;
    
    for (mol = beginMolecule(i); mol != NULL; mol = nextMolecule(i)) {
      for (integrableObject = mol->beginIntegrableObject(j); integrableObject != NULL; 
	   integrableObject = mol->nextIntegrableObject(j)) {

	ndf_local += 3;

	if (integrableObject->isDirectional()) {
	  if (integrableObject->isLinear()) {
	    ndf_local += 2;
	  } else {
	    ndf_local += 3;
	  }
	}
            
      }
    }
    
    // n_constraints is local, so subtract them on each processor
    ndf_local -= nConstraints_;

#ifdef IS_MPI
    MPI_Allreduce(&ndf_local,&ndf_,1,MPI_INT,MPI_SUM, MPI_COMM_WORLD);
#else
    ndf_ = ndf_local;
#endif

    // nZconstraints_ is global, as are the 3 COM translations for the 
    // entire system:
    ndf_ = ndf_ - 3 - nZconstraint_;

  }

  int SimInfo::getFdf() {
#ifdef IS_MPI
    MPI_Allreduce(&fdf_local,&fdf_,1,MPI_INT,MPI_SUM, MPI_COMM_WORLD);
#else
    fdf_ = fdf_local;
#endif
    return fdf_;
  }
    
  void SimInfo::calcNdfRaw() {
    int ndfRaw_local;

    MoleculeIterator i;
    vector<StuntDouble*>::iterator j;
    Molecule* mol;
    StuntDouble* integrableObject;

    // Raw degrees of freedom that we have to set
    ndfRaw_local = 0;
    
    for (mol = beginMolecule(i); mol != NULL; mol = nextMolecule(i)) {
      for (integrableObject = mol->beginIntegrableObject(j); integrableObject != NULL;
	   integrableObject = mol->nextIntegrableObject(j)) {

	ndfRaw_local += 3;

	if (integrableObject->isDirectional()) {
	  if (integrableObject->isLinear()) {
	    ndfRaw_local += 2;
	  } else {
	    ndfRaw_local += 3;
	  }
	}
            
      }
    }
    
#ifdef IS_MPI
    MPI_Allreduce(&ndfRaw_local,&ndfRaw_,1,MPI_INT,MPI_SUM, MPI_COMM_WORLD);
#else
    ndfRaw_ = ndfRaw_local;
#endif
  }

  void SimInfo::calcNdfTrans() {
    int ndfTrans_local;

    ndfTrans_local = 3 * nIntegrableObjects_ - nConstraints_;


#ifdef IS_MPI
    MPI_Allreduce(&ndfTrans_local,&ndfTrans_,1,MPI_INT,MPI_SUM, MPI_COMM_WORLD);
#else
    ndfTrans_ = ndfTrans_local;
#endif

    ndfTrans_ = ndfTrans_ - 3 - nZconstraint_;
 
  }

  void SimInfo::addInteractionPairs(Molecule* mol) {
    ForceFieldOptions& options_ = forceField_->getForceFieldOptions();
    vector<Bond*>::iterator bondIter;
    vector<Bend*>::iterator bendIter;
    vector<Torsion*>::iterator torsionIter;
    vector<Inversion*>::iterator inversionIter;
    Bond* bond;
    Bend* bend;
    Torsion* torsion;
    Inversion* inversion;
    int a;
    int b;
    int c;
    int d;

    // atomGroups can be used to add special interaction maps between
    // groups of atoms that are in two separate rigid bodies.
    // However, most site-site interactions between two rigid bodies
    // are probably not special, just the ones between the physically
    // bonded atoms.  Interactions *within* a single rigid body should
    // always be excluded.  These are done at the bottom of this
    // function.

    map<int, set<int> > atomGroups;
    Molecule::RigidBodyIterator rbIter;
    RigidBody* rb;
    Molecule::IntegrableObjectIterator ii;
    StuntDouble* integrableObject;
    
    for (integrableObject = mol->beginIntegrableObject(ii); 
         integrableObject != NULL;
         integrableObject = mol->nextIntegrableObject(ii)) {
      
      if (integrableObject->isRigidBody()) {
        rb = static_cast<RigidBody*>(integrableObject);
        vector<Atom*> atoms = rb->getAtoms();
        set<int> rigidAtoms;
        for (int i = 0; i < static_cast<int>(atoms.size()); ++i) {
          rigidAtoms.insert(atoms[i]->getGlobalIndex());
        }
        for (int i = 0; i < static_cast<int>(atoms.size()); ++i) {
          atomGroups.insert(map<int, set<int> >::value_type(atoms[i]->getGlobalIndex(), rigidAtoms));
        }      
      } else {
        set<int> oneAtomSet;
        oneAtomSet.insert(integrableObject->getGlobalIndex());
        atomGroups.insert(map<int, set<int> >::value_type(integrableObject->getGlobalIndex(), oneAtomSet));        
      }
    }  
           
    for (bond= mol->beginBond(bondIter); bond != NULL; 
         bond = mol->nextBond(bondIter)) {

      a = bond->getAtomA()->getGlobalIndex();
      b = bond->getAtomB()->getGlobalIndex();   
    
      if (options_.havevdw12scale() || options_.haveelectrostatic12scale()) {
        oneTwoInteractions_.addPair(a, b);
      } else {
        excludedInteractions_.addPair(a, b);
      }
    }

    for (bend= mol->beginBend(bendIter); bend != NULL; 
         bend = mol->nextBend(bendIter)) {

      a = bend->getAtomA()->getGlobalIndex();
      b = bend->getAtomB()->getGlobalIndex();        
      c = bend->getAtomC()->getGlobalIndex();
      
      if (options_.havevdw12scale() || options_.haveelectrostatic12scale()) {
        oneTwoInteractions_.addPair(a, b);      
        oneTwoInteractions_.addPair(b, c);
      } else {
        excludedInteractions_.addPair(a, b);
        excludedInteractions_.addPair(b, c);
      }

      if (options_.havevdw13scale() || options_.haveelectrostatic13scale()) {
        oneThreeInteractions_.addPair(a, c);      
      } else {
        excludedInteractions_.addPair(a, c);
      }
    }

    for (torsion= mol->beginTorsion(torsionIter); torsion != NULL; 
         torsion = mol->nextTorsion(torsionIter)) {

      a = torsion->getAtomA()->getGlobalIndex();
      b = torsion->getAtomB()->getGlobalIndex();        
      c = torsion->getAtomC()->getGlobalIndex();        
      d = torsion->getAtomD()->getGlobalIndex();      

      if (options_.havevdw12scale() || options_.haveelectrostatic12scale()) {
        oneTwoInteractions_.addPair(a, b);      
        oneTwoInteractions_.addPair(b, c);
        oneTwoInteractions_.addPair(c, d);
      } else {
        excludedInteractions_.addPair(a, b);
        excludedInteractions_.addPair(b, c);
        excludedInteractions_.addPair(c, d);
      }

      if (options_.havevdw13scale() || options_.haveelectrostatic13scale()) {
        oneThreeInteractions_.addPair(a, c);      
        oneThreeInteractions_.addPair(b, d);      
      } else {
        excludedInteractions_.addPair(a, c);
        excludedInteractions_.addPair(b, d);
      }

      if (options_.havevdw14scale() || options_.haveelectrostatic14scale()) {
        oneFourInteractions_.addPair(a, d);      
      } else {
        excludedInteractions_.addPair(a, d);
      }
    }

    for (inversion= mol->beginInversion(inversionIter); inversion != NULL; 
         inversion = mol->nextInversion(inversionIter)) {

      a = inversion->getAtomA()->getGlobalIndex();
      b = inversion->getAtomB()->getGlobalIndex();        
      c = inversion->getAtomC()->getGlobalIndex();        
      d = inversion->getAtomD()->getGlobalIndex();        

      if (options_.havevdw12scale() || options_.haveelectrostatic12scale()) {
        oneTwoInteractions_.addPair(a, b);      
        oneTwoInteractions_.addPair(a, c);
        oneTwoInteractions_.addPair(a, d);
      } else {
        excludedInteractions_.addPair(a, b);
        excludedInteractions_.addPair(a, c);
        excludedInteractions_.addPair(a, d);
      }

      if (options_.havevdw13scale() || options_.haveelectrostatic13scale()) {
        oneThreeInteractions_.addPair(b, c);     
        oneThreeInteractions_.addPair(b, d);     
        oneThreeInteractions_.addPair(c, d);      
      } else {
        excludedInteractions_.addPair(b, c);
        excludedInteractions_.addPair(b, d);
        excludedInteractions_.addPair(c, d);
      }
    }

    for (rb = mol->beginRigidBody(rbIter); rb != NULL; 
         rb = mol->nextRigidBody(rbIter)) {
      vector<Atom*> atoms = rb->getAtoms();
      for (int i = 0; i < static_cast<int>(atoms.size()) -1 ; ++i) {
	for (int j = i + 1; j < static_cast<int>(atoms.size()); ++j) {
	  a = atoms[i]->getGlobalIndex();
	  b = atoms[j]->getGlobalIndex();
	  excludedInteractions_.addPair(a, b);
	}
      }
    }        

  }

  void SimInfo::removeInteractionPairs(Molecule* mol) {
    ForceFieldOptions& options_ = forceField_->getForceFieldOptions();
    vector<Bond*>::iterator bondIter;
    vector<Bend*>::iterator bendIter;
    vector<Torsion*>::iterator torsionIter;
    vector<Inversion*>::iterator inversionIter;
    Bond* bond;
    Bend* bend;
    Torsion* torsion;
    Inversion* inversion;
    int a;
    int b;
    int c;
    int d;

    map<int, set<int> > atomGroups;
    Molecule::RigidBodyIterator rbIter;
    RigidBody* rb;
    Molecule::IntegrableObjectIterator ii;
    StuntDouble* integrableObject;
    
    for (integrableObject = mol->beginIntegrableObject(ii); 
         integrableObject != NULL;
         integrableObject = mol->nextIntegrableObject(ii)) {
      
      if (integrableObject->isRigidBody()) {
        rb = static_cast<RigidBody*>(integrableObject);
        vector<Atom*> atoms = rb->getAtoms();
        set<int> rigidAtoms;
        for (int i = 0; i < static_cast<int>(atoms.size()); ++i) {
          rigidAtoms.insert(atoms[i]->getGlobalIndex());
        }
        for (int i = 0; i < static_cast<int>(atoms.size()); ++i) {
          atomGroups.insert(map<int, set<int> >::value_type(atoms[i]->getGlobalIndex(), rigidAtoms));
        }      
      } else {
        set<int> oneAtomSet;
        oneAtomSet.insert(integrableObject->getGlobalIndex());
        atomGroups.insert(map<int, set<int> >::value_type(integrableObject->getGlobalIndex(), oneAtomSet));        
      }
    }  

    for (bond= mol->beginBond(bondIter); bond != NULL; 
         bond = mol->nextBond(bondIter)) {
      
      a = bond->getAtomA()->getGlobalIndex();
      b = bond->getAtomB()->getGlobalIndex();   
    
      if (options_.havevdw12scale() || options_.haveelectrostatic12scale()) {
        oneTwoInteractions_.removePair(a, b);
      } else {
        excludedInteractions_.removePair(a, b);
      }
    }

    for (bend= mol->beginBend(bendIter); bend != NULL; 
         bend = mol->nextBend(bendIter)) {

      a = bend->getAtomA()->getGlobalIndex();
      b = bend->getAtomB()->getGlobalIndex();        
      c = bend->getAtomC()->getGlobalIndex();
      
      if (options_.havevdw12scale() || options_.haveelectrostatic12scale()) {
        oneTwoInteractions_.removePair(a, b);      
        oneTwoInteractions_.removePair(b, c);
      } else {
        excludedInteractions_.removePair(a, b);
        excludedInteractions_.removePair(b, c);
      }

      if (options_.havevdw13scale() || options_.haveelectrostatic13scale()) {
        oneThreeInteractions_.removePair(a, c);      
      } else {
        excludedInteractions_.removePair(a, c);
      }
    }

    for (torsion= mol->beginTorsion(torsionIter); torsion != NULL; 
         torsion = mol->nextTorsion(torsionIter)) {

      a = torsion->getAtomA()->getGlobalIndex();
      b = torsion->getAtomB()->getGlobalIndex();        
      c = torsion->getAtomC()->getGlobalIndex();        
      d = torsion->getAtomD()->getGlobalIndex();      
  
      if (options_.havevdw12scale() || options_.haveelectrostatic12scale()) {
        oneTwoInteractions_.removePair(a, b);      
        oneTwoInteractions_.removePair(b, c);
        oneTwoInteractions_.removePair(c, d);
      } else {
        excludedInteractions_.removePair(a, b);
        excludedInteractions_.removePair(b, c);
        excludedInteractions_.removePair(c, d);
      }

      if (options_.havevdw13scale() || options_.haveelectrostatic13scale()) {
        oneThreeInteractions_.removePair(a, c);      
        oneThreeInteractions_.removePair(b, d);      
      } else {
        excludedInteractions_.removePair(a, c);
        excludedInteractions_.removePair(b, d);
      }

      if (options_.havevdw14scale() || options_.haveelectrostatic14scale()) {
        oneFourInteractions_.removePair(a, d);      
      } else {
        excludedInteractions_.removePair(a, d);
      }
    }

    for (inversion= mol->beginInversion(inversionIter); inversion != NULL; 
         inversion = mol->nextInversion(inversionIter)) {

      a = inversion->getAtomA()->getGlobalIndex();
      b = inversion->getAtomB()->getGlobalIndex();        
      c = inversion->getAtomC()->getGlobalIndex();        
      d = inversion->getAtomD()->getGlobalIndex();        

      if (options_.havevdw12scale() || options_.haveelectrostatic12scale()) {
        oneTwoInteractions_.removePair(a, b);      
        oneTwoInteractions_.removePair(a, c);
        oneTwoInteractions_.removePair(a, d);
      } else {
        excludedInteractions_.removePair(a, b);
        excludedInteractions_.removePair(a, c);
        excludedInteractions_.removePair(a, d);
      }

      if (options_.havevdw13scale() || options_.haveelectrostatic13scale()) {
        oneThreeInteractions_.removePair(b, c);     
        oneThreeInteractions_.removePair(b, d);     
        oneThreeInteractions_.removePair(c, d);      
      } else {
        excludedInteractions_.removePair(b, c);
        excludedInteractions_.removePair(b, d);
        excludedInteractions_.removePair(c, d);
      }
    }

    for (rb = mol->beginRigidBody(rbIter); rb != NULL; 
         rb = mol->nextRigidBody(rbIter)) {
      vector<Atom*> atoms = rb->getAtoms();
      for (int i = 0; i < static_cast<int>(atoms.size()) -1 ; ++i) {
	for (int j = i + 1; j < static_cast<int>(atoms.size()); ++j) {
	  a = atoms[i]->getGlobalIndex();
	  b = atoms[j]->getGlobalIndex();
	  excludedInteractions_.removePair(a, b);
	}
      }
    }        
    
  }
  
  
  void SimInfo::addMoleculeStamp(MoleculeStamp* molStamp, int nmol) {
    int curStampId;
    
    //index from 0
    curStampId = moleculeStamps_.size();

    moleculeStamps_.push_back(molStamp);
    molStampIds_.insert(molStampIds_.end(), nmol, curStampId);
  }


  /**
   * update
   *
   *  Performs the global checks and variable settings after the
   *  objects have been created.
   * 
   */
  void SimInfo::update() {   
    setupSimVariables();
    calcNdf();
    calcNdfRaw();
    calcNdfTrans();
  }
  
  /**
   * getSimulatedAtomTypes
   *
   * Returns an STL set of AtomType* that are actually present in this
   * simulation.  Must query all processors to assemble this information.
   * 
   */
  set<AtomType*> SimInfo::getSimulatedAtomTypes() {
    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::AtomIterator ai;
    Atom* atom;
    set<AtomType*> atomTypes;
    
    for(mol = beginMolecule(mi); mol != NULL; mol = nextMolecule(mi)) {      
      for(atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
	atomTypes.insert(atom->getAtomType());
      }      
    }    

#ifdef IS_MPI

    // loop over the found atom types on this processor, and add their
    // numerical idents to a vector:

    vector<int> foundTypes;
    set<AtomType*>::iterator i;
    for (i = atomTypes.begin(); i != atomTypes.end(); ++i) 
      foundTypes.push_back( (*i)->getIdent() );

    // count_local holds the number of found types on this processor
    int count_local = foundTypes.size();

    // count holds the total number of found types on all processors
    // (some will be redundant with the ones found locally):
    int count;
    MPI::COMM_WORLD.Allreduce(&count_local, &count, 1, MPI::INT, MPI::SUM);

    // create a vector to hold the globally found types, and resize it:
    vector<int> ftGlobal;
    ftGlobal.resize(count);
    vector<int> counts;

    int nproc = MPI::COMM_WORLD.Get_size();
    counts.resize(nproc);
    vector<int> disps;
    disps.resize(nproc);

    // now spray out the foundTypes to all the other processors:
    
    MPI::COMM_WORLD.Allgatherv(&foundTypes[0], count_local, MPI::INT, 
                               &ftGlobal[0], &counts[0], &disps[0], MPI::INT);

    // foundIdents is a stl set, so inserting an already found ident
    // will have no effect.
    set<int> foundIdents;
    vector<int>::iterator j;
    for (j = ftGlobal.begin(); j != ftGlobal.end(); ++j)
      foundIdents.insert((*j));
    
    // now iterate over the foundIdents and get the actual atom types 
    // that correspond to these:
    set<int>::iterator it;
    for (it = foundIdents.begin(); it != foundIdents.end(); ++it)
      atomTypes.insert( forceField_->getAtomType((*it)) );
 
#endif
    
    return atomTypes;        
  }

  void SimInfo::setupSimVariables() {
    useAtomicVirial_ = simParams_->getUseAtomicVirial();
    // we only call setAccumulateBoxDipole if the accumulateBoxDipole parameter is true
    calcBoxDipole_ = false;
    if ( simParams_->haveAccumulateBoxDipole() ) 
      if ( simParams_->getAccumulateBoxDipole() ) {
	calcBoxDipole_ = true;       
      }

    set<AtomType*>::iterator i;
    set<AtomType*> atomTypes;
    atomTypes = getSimulatedAtomTypes();    
    int usesElectrostatic = 0;
    int usesMetallic = 0;
    int usesDirectional = 0;
    //loop over all of the atom types
    for (i = atomTypes.begin(); i != atomTypes.end(); ++i) {
      usesElectrostatic |= (*i)->isElectrostatic();
      usesMetallic |= (*i)->isMetal();
      usesDirectional |= (*i)->isDirectional();
    }

#ifdef IS_MPI    
    int temp;
    temp = usesDirectional;
    MPI_Allreduce(&temp, &usesDirectionalAtoms_, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);    

    temp = usesMetallic;
    MPI_Allreduce(&temp, &usesMetallicAtoms_, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);    

    temp = usesElectrostatic;
    MPI_Allreduce(&temp, &usesElectrostaticAtoms_, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD); 
#endif
    fInfo_.SIM_uses_PBC = usesPeriodicBoundaries_;    
    fInfo_.SIM_uses_DirectionalAtoms = usesDirectionalAtoms_;
    fInfo_.SIM_uses_MetallicAtoms = usesMetallicAtoms_;
    fInfo_.SIM_requires_SkipCorrection = usesElectrostaticAtoms_;
    fInfo_.SIM_requires_SelfCorrection = usesElectrostaticAtoms_;
    fInfo_.SIM_uses_AtomicVirial = usesAtomicVirial_;
  }

  void SimInfo::setupFortran() {
    int isError;
    int nExclude, nOneTwo, nOneThree, nOneFour;
    vector<int> fortranGlobalGroupMembership;
    
    isError = 0;

    //globalGroupMembership_ is filled by SimCreator    
    for (int i = 0; i < nGlobalAtoms_; i++) {
      fortranGlobalGroupMembership.push_back(globalGroupMembership_[i] + 1);
    }

    //calculate mass ratio of cutoff group
    vector<RealType> mfact;
    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::CutoffGroupIterator ci;
    CutoffGroup* cg;
    Molecule::AtomIterator ai;
    Atom* atom;
    RealType totalMass;

    //to avoid memory reallocation, reserve enough space for mfact
    mfact.reserve(getNCutoffGroups());
    
    for(mol = beginMolecule(mi); mol != NULL; mol = nextMolecule(mi)) {        
      for (cg = mol->beginCutoffGroup(ci); cg != NULL; cg = mol->nextCutoffGroup(ci)) {

	totalMass = cg->getMass();
	for(atom = cg->beginAtom(ai); atom != NULL; atom = cg->nextAtom(ai)) {
	  // Check for massless groups - set mfact to 1 if true
	  if (totalMass != 0)
	    mfact.push_back(atom->getMass()/totalMass);
	  else
	    mfact.push_back( 1.0 );
	}
      }       
    }

    //fill ident array of local atoms (it is actually ident of
    //AtomType, it is so confusing !!!)
    vector<int> identArray;

    //to avoid memory reallocation, reserve enough space identArray
    identArray.reserve(getNAtoms());
    
    for(mol = beginMolecule(mi); mol != NULL; mol = nextMolecule(mi)) {        
      for(atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
	identArray.push_back(atom->getIdent());
      }
    }    

    //fill molMembershipArray
    //molMembershipArray is filled by SimCreator    
    vector<int> molMembershipArray(nGlobalAtoms_);
    for (int i = 0; i < nGlobalAtoms_; i++) {
      molMembershipArray[i] = globalMolMembership_[i] + 1;
    }
    
    //setup fortran simulation

    nExclude = excludedInteractions_.getSize();
    nOneTwo = oneTwoInteractions_.getSize();
    nOneThree = oneThreeInteractions_.getSize();
    nOneFour = oneFourInteractions_.getSize();

    int* excludeList = excludedInteractions_.getPairList();
    int* oneTwoList = oneTwoInteractions_.getPairList();
    int* oneThreeList = oneThreeInteractions_.getPairList();
    int* oneFourList = oneFourInteractions_.getPairList();

    setFortranSim( &fInfo_, &nGlobalAtoms_, &nAtoms_, &identArray[0], 
                   &nExclude, excludeList, 
                   &nOneTwo, oneTwoList,
                   &nOneThree, oneThreeList,
                   &nOneFour, oneFourList,
                   &molMembershipArray[0], &mfact[0], &nCutoffGroups_, 
                   &fortranGlobalGroupMembership[0], &isError); 
    
    if( isError ){
      
      sprintf( painCave.errMsg,
	       "There was an error setting the simulation information in fortran.\n" );
      painCave.isFatal = 1;
      painCave.severity = OPENMD_ERROR;
      simError();
    }
    
    
    sprintf( checkPointMsg,
	     "succesfully sent the simulation information to fortran.\n");
    
    errorCheckPoint();
    
    // Setup number of neighbors in neighbor list if present
    if (simParams_->haveNeighborListNeighbors()) {
      int nlistNeighbors = simParams_->getNeighborListNeighbors();
      setNeighbors(&nlistNeighbors);
    }
   
#ifdef IS_MPI    
    //SimInfo is responsible for creating localToGlobalAtomIndex and
    //localToGlobalGroupIndex
    vector<int> localToGlobalAtomIndex(getNAtoms(), 0);
    vector<int> localToGlobalCutoffGroupIndex;
    mpiSimData parallelData;

    for (mol = beginMolecule(mi); mol != NULL; mol  = nextMolecule(mi)) {

      //local index(index in DataStorge) of atom is important
      for (atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
	localToGlobalAtomIndex[atom->getLocalIndex()] = atom->getGlobalIndex() + 1;
      }

      //local index of cutoff group is trivial, it only depends on the order of travesing
      for (cg = mol->beginCutoffGroup(ci); cg != NULL; cg = mol->nextCutoffGroup(ci)) {
	localToGlobalCutoffGroupIndex.push_back(cg->getGlobalIndex() + 1);
      }        
        
    }

    //fill up mpiSimData struct
    parallelData.nMolGlobal = getNGlobalMolecules();
    parallelData.nMolLocal = getNMolecules();
    parallelData.nAtomsGlobal = getNGlobalAtoms();
    parallelData.nAtomsLocal = getNAtoms();
    parallelData.nGroupsGlobal = getNGlobalCutoffGroups();
    parallelData.nGroupsLocal = getNCutoffGroups();
    parallelData.myNode = worldRank;
    MPI_Comm_size(MPI_COMM_WORLD, &(parallelData.nProcessors));

    //pass mpiSimData struct and index arrays to fortran
    setFsimParallel(&parallelData, &(parallelData.nAtomsLocal),
                    &localToGlobalAtomIndex[0],  &(parallelData.nGroupsLocal),
                    &localToGlobalCutoffGroupIndex[0], &isError);

    if (isError) {
      sprintf(painCave.errMsg,
	      "mpiRefresh errror: fortran didn't like something we gave it.\n");
      painCave.isFatal = 1;
      simError();
    }

    sprintf(checkPointMsg, " mpiRefresh successful.\n");
    errorCheckPoint();
#endif
    fortranInitialized_ = true;
  }

  void SimInfo::addProperty(GenericData* genData) {
    properties_.addProperty(genData);  
  }

  void SimInfo::removeProperty(const string& propName) {
    properties_.removeProperty(propName);  
  }

  void SimInfo::clearProperties() {
    properties_.clearProperties(); 
  }

  vector<string> SimInfo::getPropertyNames() {
    return properties_.getPropertyNames();  
  }
      
  vector<GenericData*> SimInfo::getProperties() { 
    return properties_.getProperties(); 
  }

  GenericData* SimInfo::getPropertyByName(const string& propName) {
    return properties_.getPropertyByName(propName); 
  }

  void SimInfo::setSnapshotManager(SnapshotManager* sman) {
    if (sman_ == sman) {
      return;
    }    
    delete sman_;
    sman_ = sman;

    Molecule* mol;
    RigidBody* rb;
    Atom* atom;
    SimInfo::MoleculeIterator mi;
    Molecule::RigidBodyIterator rbIter;
    Molecule::AtomIterator atomIter;;
 
    for (mol = beginMolecule(mi); mol != NULL; mol = nextMolecule(mi)) {
        
      for (atom = mol->beginAtom(atomIter); atom != NULL; atom = mol->nextAtom(atomIter)) {
	atom->setSnapshotManager(sman_);
      }
        
      for (rb = mol->beginRigidBody(rbIter); rb != NULL; rb = mol->nextRigidBody(rbIter)) {
	rb->setSnapshotManager(sman_);
      }
    }    
    
  }

  Vector3d SimInfo::getComVel(){ 
    SimInfo::MoleculeIterator i;
    Molecule* mol;

    Vector3d comVel(0.0);
    RealType totalMass = 0.0;
    
 
    for (mol = beginMolecule(i); mol != NULL; mol = nextMolecule(i)) {
      RealType mass = mol->getMass();
      totalMass += mass;
      comVel += mass * mol->getComVel();
    }  

#ifdef IS_MPI
    RealType tmpMass = totalMass;
    Vector3d tmpComVel(comVel);    
    MPI_Allreduce(&tmpMass,&totalMass,1,MPI_REALTYPE,MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(tmpComVel.getArrayPointer(), comVel.getArrayPointer(),3,MPI_REALTYPE,MPI_SUM, MPI_COMM_WORLD);
#endif

    comVel /= totalMass;

    return comVel;
  }

  Vector3d SimInfo::getCom(){ 
    SimInfo::MoleculeIterator i;
    Molecule* mol;

    Vector3d com(0.0);
    RealType totalMass = 0.0;
     
    for (mol = beginMolecule(i); mol != NULL; mol = nextMolecule(i)) {
      RealType mass = mol->getMass();
      totalMass += mass;
      com += mass * mol->getCom();
    }  

#ifdef IS_MPI
    RealType tmpMass = totalMass;
    Vector3d tmpCom(com);    
    MPI_Allreduce(&tmpMass,&totalMass,1,MPI_REALTYPE,MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(tmpCom.getArrayPointer(), com.getArrayPointer(),3,MPI_REALTYPE,MPI_SUM, MPI_COMM_WORLD);
#endif

    com /= totalMass;

    return com;

  }        

  ostream& operator <<(ostream& o, SimInfo& info) {

    return o;
  }
   
   
   /* 
   Returns center of mass and center of mass velocity in one function call.
   */
   
   void SimInfo::getComAll(Vector3d &com, Vector3d &comVel){ 
      SimInfo::MoleculeIterator i;
      Molecule* mol;
      
    
      RealType totalMass = 0.0;
    

      for (mol = beginMolecule(i); mol != NULL; mol = nextMolecule(i)) {
         RealType mass = mol->getMass();
         totalMass += mass;
         com += mass * mol->getCom();
         comVel += mass * mol->getComVel();           
      }  
      
#ifdef IS_MPI
      RealType tmpMass = totalMass;
      Vector3d tmpCom(com);  
      Vector3d tmpComVel(comVel);
      MPI_Allreduce(&tmpMass,&totalMass,1,MPI_REALTYPE,MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(tmpCom.getArrayPointer(), com.getArrayPointer(),3,MPI_REALTYPE,MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(tmpComVel.getArrayPointer(), comVel.getArrayPointer(),3,MPI_REALTYPE,MPI_SUM, MPI_COMM_WORLD);
#endif
      
      com /= totalMass;
      comVel /= totalMass;
   }        
   
   /* 
   Return intertia tensor for entire system and angular momentum Vector.


       [  Ixx -Ixy  -Ixz ]
    J =| -Iyx  Iyy  -Iyz |
       [ -Izx -Iyz   Izz ]
    */

   void SimInfo::getInertiaTensor(Mat3x3d &inertiaTensor, Vector3d &angularMomentum){
      
 
      RealType xx = 0.0;
      RealType yy = 0.0;
      RealType zz = 0.0;
      RealType xy = 0.0;
      RealType xz = 0.0;
      RealType yz = 0.0;
      Vector3d com(0.0);
      Vector3d comVel(0.0);
      
      getComAll(com, comVel);
      
      SimInfo::MoleculeIterator i;
      Molecule* mol;
      
      Vector3d thisq(0.0);
      Vector3d thisv(0.0);

      RealType thisMass = 0.0;
     
      
      
   
      for (mol = beginMolecule(i); mol != NULL; mol = nextMolecule(i)) {
        
         thisq = mol->getCom()-com;
         thisv = mol->getComVel()-comVel;
         thisMass = mol->getMass();
         // Compute moment of intertia coefficients.
         xx += thisq[0]*thisq[0]*thisMass;
         yy += thisq[1]*thisq[1]*thisMass;
         zz += thisq[2]*thisq[2]*thisMass;
         
         // compute products of intertia
         xy += thisq[0]*thisq[1]*thisMass;
         xz += thisq[0]*thisq[2]*thisMass;
         yz += thisq[1]*thisq[2]*thisMass;
            
         angularMomentum += cross( thisq, thisv ) * thisMass;
            
      }  
      
      
      inertiaTensor(0,0) = yy + zz;
      inertiaTensor(0,1) = -xy;
      inertiaTensor(0,2) = -xz;
      inertiaTensor(1,0) = -xy;
      inertiaTensor(1,1) = xx + zz;
      inertiaTensor(1,2) = -yz;
      inertiaTensor(2,0) = -xz;
      inertiaTensor(2,1) = -yz;
      inertiaTensor(2,2) = xx + yy;
      
#ifdef IS_MPI
      Mat3x3d tmpI(inertiaTensor);
      Vector3d tmpAngMom;
      MPI_Allreduce(tmpI.getArrayPointer(), inertiaTensor.getArrayPointer(),9,MPI_REALTYPE,MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(tmpAngMom.getArrayPointer(), angularMomentum.getArrayPointer(),3,MPI_REALTYPE,MPI_SUM, MPI_COMM_WORLD);
#endif
               
      return;
   }

   //Returns the angular momentum of the system
   Vector3d SimInfo::getAngularMomentum(){
      
      Vector3d com(0.0);
      Vector3d comVel(0.0);
      Vector3d angularMomentum(0.0);
      
      getComAll(com,comVel);
      
      SimInfo::MoleculeIterator i;
      Molecule* mol;
      
      Vector3d thisr(0.0);
      Vector3d thisp(0.0);
      
      RealType thisMass;
      
      for (mol = beginMolecule(i); mol != NULL; mol = nextMolecule(i)) {         
        thisMass = mol->getMass(); 
	thisr = mol->getCom()-com;
	thisp = (mol->getComVel()-comVel)*thisMass;
         
	angularMomentum += cross( thisr, thisp );
         
      }  
       
#ifdef IS_MPI
      Vector3d tmpAngMom;
      MPI_Allreduce(tmpAngMom.getArrayPointer(), angularMomentum.getArrayPointer(),3,MPI_REALTYPE,MPI_SUM, MPI_COMM_WORLD);
#endif
      
      return angularMomentum;
   }
   
  StuntDouble* SimInfo::getIOIndexToIntegrableObject(int index) {
    return IOIndexToIntegrableObject.at(index);
  }
  
  void SimInfo::setIOIndexToIntegrableObject(const vector<StuntDouble*>& v) {
    IOIndexToIntegrableObject= v;
  }

  /* Returns the Volume of the simulation based on a ellipsoid with semi-axes 
     based on the radius of gyration V=4/3*Pi*R_1*R_2*R_3
     where R_i are related to the principle inertia moments R_i = sqrt(C*I_i/N), this reduces to 
     V = 4/3*Pi*(C/N)^3/2*sqrt(det(I)). See S.E. Baltazar et. al. Comp. Mat. Sci. 37 (2006) 526-536.
  */
  void SimInfo::getGyrationalVolume(RealType &volume){
    Mat3x3d intTensor;
    RealType det;
    Vector3d dummyAngMom; 
    RealType sysconstants;
    RealType geomCnst;

    geomCnst = 3.0/2.0;
    /* Get the inertial tensor and angular momentum for free*/
    getInertiaTensor(intTensor,dummyAngMom);
    
    det = intTensor.determinant();
    sysconstants = geomCnst/(RealType)nGlobalIntegrableObjects_;
    volume = 4.0/3.0*NumericConstant::PI*pow(sysconstants,3.0/2.0)*sqrt(det);
    return;
  }

  void SimInfo::getGyrationalVolume(RealType &volume, RealType &detI){
    Mat3x3d intTensor;
    Vector3d dummyAngMom; 
    RealType sysconstants;
    RealType geomCnst;

    geomCnst = 3.0/2.0;
    /* Get the inertial tensor and angular momentum for free*/
    getInertiaTensor(intTensor,dummyAngMom);
    
    detI = intTensor.determinant();
    sysconstants = geomCnst/(RealType)nGlobalIntegrableObjects_;
    volume = 4.0/3.0*NumericConstant::PI*pow(sysconstants,3.0/2.0)*sqrt(detI);
    return;
  }
/*
   void SimInfo::setStuntDoubleFromGlobalIndex(vector<StuntDouble*> v) {
      assert( v.size() == nAtoms_ + nRigidBodies_);
      sdByGlobalIndex_ = v;
    }

    StuntDouble* SimInfo::getStuntDoubleFromGlobalIndex(int index) {
      //assert(index < nAtoms_ + nRigidBodies_);
      return sdByGlobalIndex_.at(index);
    }   
*/   
  int SimInfo::getNGlobalConstraints() {
    int nGlobalConstraints;
#ifdef IS_MPI
    MPI_Allreduce(&nConstraints_, &nGlobalConstraints, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);    
#else
    nGlobalConstraints =  nConstraints_;
#endif
    return nGlobalConstraints;
  }

}//end namespace OpenMD

