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
 * @file MoleculeCreator.cpp
 * @author tlin
 * @date 11/04/2004
 * @time 13:44am
 * @version 1.0
 */

#include <cassert>
#include <typeinfo>
#include <set>

#include "brains/MoleculeCreator.hpp"
#include "primitives/GhostBend.hpp"
#include "primitives/GhostTorsion.hpp"
#include "types/AtomType.hpp"
#include "types/FixedBondType.hpp"
#include "utils/simError.h"
#include "utils/StringUtils.hpp"

namespace OpenMD {
  
  Molecule* MoleculeCreator::createMolecule(ForceField* ff, 
                                            MoleculeStamp *molStamp,
					    int stampId, int globalIndex, 
                                            LocalIndexManager* localIndexMan) {
    Molecule* mol = new Molecule(stampId, globalIndex, molStamp->getName());
    
    //create atoms
    Atom* atom;
    AtomStamp* currentAtomStamp;
    int nAtom = molStamp->getNAtoms();
    for (int i = 0; i < nAtom; ++i) {
      currentAtomStamp = molStamp->getAtomStamp(i);
      atom = createAtom(ff, mol, currentAtomStamp, localIndexMan);
      mol->addAtom(atom);
    }

    //create rigidbodies
    RigidBody* rb;
    RigidBodyStamp * currentRigidBodyStamp;
    int nRigidbodies = molStamp->getNRigidBodies();

    for (int i = 0; i < nRigidbodies; ++i) {
      currentRigidBodyStamp = molStamp->getRigidBodyStamp(i);
      rb = createRigidBody(molStamp, mol, currentRigidBodyStamp, 
                           localIndexMan);
      mol->addRigidBody(rb);
    }
    
    //create bonds
    Bond* bond;
    BondStamp* currentBondStamp;
    int nBonds = molStamp->getNBonds();

    for (int i = 0; i < nBonds; ++i) {
      currentBondStamp = molStamp->getBondStamp(i);
      bond = createBond(ff, mol, currentBondStamp);
      mol->addBond(bond);
    }

    //create bends
    Bend* bend;
    BendStamp* currentBendStamp;
    int nBends = molStamp->getNBends();
    for (int i = 0; i < nBends; ++i) {
      currentBendStamp = molStamp->getBendStamp(i);
      bend = createBend(ff, mol, currentBendStamp);
      mol->addBend(bend);
    }

    //create torsions
    Torsion* torsion;
    TorsionStamp* currentTorsionStamp;
    int nTorsions = molStamp->getNTorsions();
    for (int i = 0; i < nTorsions; ++i) {
      currentTorsionStamp = molStamp->getTorsionStamp(i);
      torsion = createTorsion(ff, mol, currentTorsionStamp);
      mol->addTorsion(torsion);
    }

    //create inversions
    Inversion* inversion;
    InversionStamp* currentInversionStamp;
    int nInversions = molStamp->getNInversions();
    for (int i = 0; i < nInversions; ++i) {
      currentInversionStamp = molStamp->getInversionStamp(i);
      inversion = createInversion(ff, mol, currentInversionStamp);
      if (inversion != NULL ) {
        mol->addInversion(inversion);
      }
    }

    //create cutoffGroups
    CutoffGroup* cutoffGroup;
    CutoffGroupStamp* currentCutoffGroupStamp;
    int nCutoffGroups = molStamp->getNCutoffGroups();
    for (int i = 0; i < nCutoffGroups; ++i) {
      currentCutoffGroupStamp = molStamp->getCutoffGroupStamp(i);
      cutoffGroup = createCutoffGroup(mol, currentCutoffGroupStamp, localIndexMan);
      mol->addCutoffGroup(cutoffGroup);
    }

    //every free atom is a cutoff group    
    std::vector<Atom*> freeAtoms;
    std::vector<Atom*>::iterator ai;
    std::vector<Atom*>::iterator fai;

    //add all atoms into allAtoms set
    for(atom = mol->beginAtom(fai); atom != NULL; atom = mol->nextAtom(fai)) {
      freeAtoms.push_back(atom);
    }

    Molecule::CutoffGroupIterator ci;
    CutoffGroup* cg;
    
    for (cg = mol->beginCutoffGroup(ci); cg != NULL; 
         cg = mol->nextCutoffGroup(ci)) {
      
      for(atom = cg->beginAtom(ai); atom != NULL; atom = cg->nextAtom(ai)) {
        //erase the atoms belong to cutoff groups from freeAtoms vector
        freeAtoms.erase(std::remove(freeAtoms.begin(), freeAtoms.end(), atom),
                        freeAtoms.end());
      }      
    }       
    
    // loop over the free atoms and then create one cutoff group for
    // every single free atom
    
    for (fai = freeAtoms.begin(); fai != freeAtoms.end(); ++fai) {
      cutoffGroup = createCutoffGroup(mol, *fai, localIndexMan);
      mol->addCutoffGroup(cutoffGroup);
    }
    //create constraints
    createConstraintPair(mol);
    createConstraintElem(mol);
    
    //the construction of this molecule is finished
    mol->complete();
    
    return mol;
  }    


  Atom* MoleculeCreator::createAtom(ForceField* ff, Molecule* mol, 
                                    AtomStamp* stamp,
                                    LocalIndexManager* localIndexMan) {
    AtomType * atomType;
    Atom* atom;

    atomType =  ff->getAtomType(stamp->getType());
    
    if (atomType == NULL) {
      sprintf(painCave.errMsg, "Can not find Matching Atom Type for[%s]",
	      stamp->getType().c_str());

      painCave.isFatal = 1;
      simError();
    }
    
    //below code still have some kind of hard-coding smell
    if (atomType->isDirectional()){
      DirectionalAtom* dAtom;
      dAtom = new DirectionalAtom(atomType);
      atom = dAtom;    
    }
    else{
      atom = new Atom(atomType);
    }

    atom->setLocalIndex(localIndexMan->getNextAtomIndex());

    return atom;
  }
  
  RigidBody* MoleculeCreator::createRigidBody(MoleculeStamp *molStamp, 
                                              Molecule* mol,
					      RigidBodyStamp* rbStamp, 
					      LocalIndexManager* localIndexMan) {
    Atom* atom;
    int nAtoms;
    Vector3d refCoor;
    AtomStamp* atomStamp;
    
    RigidBody* rb = new RigidBody();
    nAtoms = rbStamp->getNMembers();    
    for (int i = 0; i < nAtoms; ++i) {
      //rbStamp->getMember(i) return the local index of current atom
      //inside the molecule.  It is not the same as local index of
      //atom which is the index of atom at DataStorage class
      atom = mol->getAtomAt(rbStamp->getMemberAt(i));
      atomStamp= molStamp->getAtomStamp(rbStamp->getMemberAt(i));    
      rb->addAtom(atom, atomStamp);
    }

    //after all of the atoms are added, we need to calculate the
    //reference coordinates
    rb->calcRefCoords();

    //set the local index of this rigid body, global index will be set later
    rb->setLocalIndex(localIndexMan->getNextRigidBodyIndex());

    //the rule for naming rigidbody MoleculeName_RB_Integer
    //The first part is the name of the molecule
    //The second part is alway fixed as "RB"
    //The third part is the index of the rigidbody defined in meta-data file
    //For example, Butane_RB_0 is a valid rigid body name of butane molecule
    /**@todo replace itoa by lexi_cast */
    std::string s = OpenMD_itoa(mol->getNRigidBodies(), 10);
    rb->setType(mol->getType() + "_RB_" + s.c_str());

    return rb;
  }    

  Bond* MoleculeCreator::createBond(ForceField* ff, Molecule* mol, 
                                    BondStamp* stamp) {
    BondType* bondType;
    Atom* atomA;
    Atom* atomB;
    
    atomA = mol->getAtomAt(stamp->getA());
    atomB = mol->getAtomAt(stamp->getB());
    
    assert( atomA && atomB);
    
    bondType = ff->getBondType(atomA->getType(), atomB->getType());

    if (bondType == NULL) {
      sprintf(painCave.errMsg, "Can not find Matching Bond Type for[%s, %s]",
	      atomA->getType().c_str(),
	      atomB->getType().c_str());
      
      painCave.isFatal = 1;
      simError();
    }
    return new Bond(atomA, atomB, bondType);    
  }    
  
  Bend* MoleculeCreator::createBend(ForceField* ff, Molecule* mol, 
                                    BendStamp* stamp) {
    Bend* bend = NULL; 
    std::vector<int> bendAtoms = stamp->getMembers(); 
    if (bendAtoms.size() == 3) {
      Atom* atomA = mol->getAtomAt(bendAtoms[0]);
      Atom* atomB = mol->getAtomAt(bendAtoms[1]);
      Atom* atomC = mol->getAtomAt(bendAtoms[2]);
      
      assert( atomA && atomB && atomC);
      
      BendType* bendType = ff->getBendType(atomA->getType().c_str(), 
                                           atomB->getType().c_str(), 
                                           atomC->getType().c_str());
      
      if (bendType == NULL) {
        sprintf(painCave.errMsg, "Can not find Matching Bend Type for[%s, %s, %s]",
                atomA->getType().c_str(),
                atomB->getType().c_str(),
                atomC->getType().c_str());
        
        painCave.isFatal = 1;
        simError();
      }
      
      bend = new Bend(atomA, atomB, atomC, bendType);
    } else if ( bendAtoms.size() == 2 && stamp->haveGhostVectorSource()) {
      int ghostIndex = stamp->getGhostVectorSource();
      int normalIndex = ghostIndex != bendAtoms[0] ? bendAtoms[0] : bendAtoms[1]; 
      Atom* normalAtom = mol->getAtomAt(normalIndex) ;        
      DirectionalAtom* ghostAtom = dynamic_cast<DirectionalAtom*>(mol->getAtomAt(ghostIndex));
      if (ghostAtom == NULL) {
	sprintf(painCave.errMsg, "Can not cast Atom to DirectionalAtom");
	painCave.isFatal = 1;
	simError();
      }
                
      BendType* bendType = ff->getBendType(normalAtom->getType(), ghostAtom->getType(), "GHOST");

      if (bendType == NULL) {
	sprintf(painCave.errMsg, "Can not find Matching Bend Type for[%s, %s, %s]",
		normalAtom->getType().c_str(),
		ghostAtom->getType().c_str(),
		"GHOST");

	painCave.isFatal = 1;
	simError();
      }
      
      bend = new GhostBend(normalAtom, ghostAtom, bendType);       
      
    } 
    
    return bend;
  }    

  Torsion* MoleculeCreator::createTorsion(ForceField* ff, Molecule* mol, 
                                          TorsionStamp* stamp) {

    Torsion* torsion = NULL;
    std::vector<int> torsionAtoms = stamp->getMembers();
    if (torsionAtoms.size() < 3) {
	return torsion;
    }

    Atom* atomA = mol->getAtomAt(torsionAtoms[0]);
    Atom* atomB = mol->getAtomAt(torsionAtoms[1]);
    Atom* atomC = mol->getAtomAt(torsionAtoms[2]);

    if (torsionAtoms.size() == 4) {
      Atom* atomD = mol->getAtomAt(torsionAtoms[3]);

      assert(atomA && atomB && atomC && atomD);
        
      TorsionType* torsionType = ff->getTorsionType(atomA->getType(), 
                                                    atomB->getType(), 
						    atomC->getType(), 
                                                    atomD->getType());
      if (torsionType == NULL) {
	sprintf(painCave.errMsg, "Can not find Matching Torsion Type for[%s, %s, %s, %s]",
		atomA->getType().c_str(),
		atomB->getType().c_str(),
		atomC->getType().c_str(),
		atomD->getType().c_str());
        
	painCave.isFatal = 1;
	simError();
      }
      
      torsion = new Torsion(atomA, atomB, atomC, atomD, torsionType);       
    }
    else {
      
      DirectionalAtom* dAtom = dynamic_cast<DirectionalAtom*>(mol->getAtomAt(stamp->getGhostVectorSource()));
      if (dAtom == NULL) {
	sprintf(painCave.errMsg, "Can not cast Atom to DirectionalAtom");
	painCave.isFatal = 1;
	simError();
      }        
      
      TorsionType* torsionType = ff->getTorsionType(atomA->getType(), atomB->getType(), 
						    atomC->getType(), "GHOST");
      
      if (torsionType == NULL) {
	sprintf(painCave.errMsg, "Can not find Matching Torsion Type for[%s, %s, %s, %s]",
		atomA->getType().c_str(),
		atomB->getType().c_str(),
		atomC->getType().c_str(),
		"GHOST");
        
	painCave.isFatal = 1;
	simError();
      }
      
      torsion = new GhostTorsion(atomA, atomB, dAtom, torsionType);               
    }
    
    return torsion;
  }    

  Inversion* MoleculeCreator::createInversion(ForceField* ff, Molecule* mol, 
                                              InversionStamp* stamp) {
    
    Inversion* inversion = NULL;
    int center = stamp->getCenter();
    std::vector<int> satellites = stamp->getSatellites();
    if (satellites.size() != 3) {
	return inversion;
    }

    Atom* atomA = mol->getAtomAt(center);
    Atom* atomB = mol->getAtomAt(satellites[0]);
    Atom* atomC = mol->getAtomAt(satellites[1]);
    Atom* atomD = mol->getAtomAt(satellites[2]);
      
    assert(atomA && atomB && atomC && atomD);
    
    InversionType* inversionType = ff->getInversionType(atomA->getType(), 
                                                        atomB->getType(), 
                                                        atomC->getType(), 
                                                        atomD->getType());

    if (inversionType == NULL) {
      sprintf(painCave.errMsg, "No Matching Inversion Type for[%s, %s, %s, %s]\n"
              "\t(May not be a problem: not all inversions are parametrized)\n",
              atomA->getType().c_str(),
              atomB->getType().c_str(),
              atomC->getType().c_str(),
              atomD->getType().c_str());
      
      painCave.isFatal = 0;
      painCave.severity = OPENMD_INFO;
      simError();
      return NULL;
    } else {
      
      inversion = new Inversion(atomA, atomB, atomC, atomD, inversionType);
      return inversion;
    }
  }
  

  CutoffGroup* MoleculeCreator::createCutoffGroup(Molecule* mol, 
                                                  CutoffGroupStamp* stamp,
                                                  LocalIndexManager* localIndexMan) {
    int nAtoms;
    CutoffGroup* cg;
    Atom* atom;
    cg = new CutoffGroup();
    
    nAtoms = stamp->getNMembers();
    for (int i =0; i < nAtoms; ++i) {
      atom = mol->getAtomAt(stamp->getMemberAt(i));
      assert(atom);
      cg->addAtom(atom);
    }
    
    //set the local index of this cutoffGroup, global index will be set later
    cg->setLocalIndex(localIndexMan->getNextCutoffGroupIndex());
    
    return cg;
  }    
  
  CutoffGroup* MoleculeCreator::createCutoffGroup(Molecule * mol, Atom* atom,
                                                  LocalIndexManager* localIndexMan) {
    CutoffGroup* cg;
    cg  = new CutoffGroup();
    cg->addAtom(atom);

    //set the local index of this cutoffGroup, global index will be set later
    cg->setLocalIndex(localIndexMan->getNextCutoffGroupIndex());

    return cg;
  }

  void MoleculeCreator::createConstraintPair(Molecule* mol) {

    //add bond constraints
    Molecule::BondIterator bi;
    Bond* bond;
    for (bond = mol->beginBond(bi); bond != NULL; bond = mol->nextBond(bi)) {
        
      BondType* bt = bond->getBondType();

      //class Parent1 {};
      //class Child1 : public Parent {};
      //class Child2 : public Parent {};
      //Child1* ch1 = new Child1();
      //Child2* ch2 = dynamic_cast<Child2*>(ch1); 
      //the dynamic_cast is succeed in above line. A compiler bug?        

      if (typeid(FixedBondType) == typeid(*bt)) {
	FixedBondType* fbt = dynamic_cast<FixedBondType*>(bt);

	ConstraintElem* consElemA = new ConstraintElem(bond->getAtomA());
	ConstraintElem* consElemB = new ConstraintElem(bond->getAtomB());            
	ConstraintPair* consPair = new ConstraintPair(consElemA, consElemB, fbt->getEquilibriumBondLength());
	mol->addConstraintPair(consPair);
      }
    }

    //rigidbody -- rigidbody constraint is not support yet
  }

  void MoleculeCreator::createConstraintElem(Molecule* mol) {

    ConstraintPair* consPair;
    Molecule::ConstraintPairIterator cpi;
    std::set<StuntDouble*> sdSet;
    for (consPair = mol->beginConstraintPair(cpi); consPair != NULL; consPair = mol->nextConstraintPair(cpi)) {

      StuntDouble* sdA = consPair->getConsElem1()->getStuntDouble();            
      if (sdSet.find(sdA) == sdSet.end()){
	sdSet.insert(sdA);
	mol->addConstraintElem(new ConstraintElem(sdA));
      }

      StuntDouble* sdB = consPair->getConsElem2()->getStuntDouble();            
      if (sdSet.find(sdB) == sdSet.end()){
	sdSet.insert(sdB);
	mol->addConstraintElem(new ConstraintElem(sdB));
      }
        
    }

  }
    
}
