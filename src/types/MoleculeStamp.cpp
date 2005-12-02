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
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>

#include "types/MoleculeStamp.hpp"

namespace oopse {
MoleculeStamp::MoleculeStamp() {
    DefineParameter(Name, "name");
    
    deprecatedKeywords_.insert("nAtoms");
    deprecatedKeywords_.insert("nBonds");
    deprecatedKeywords_.insert("nBends");
    deprecatedKeywords_.insert("nTorsions");
    deprecatedKeywords_.insert("nRigidBodies");
    deprecatedKeywords_.insert("nCutoffGroups");
    
}

MoleculeStamp::~MoleculeStamp() {

}

bool MoleculeStamp::addAtomStamp( AtomStamp* atom) {
    bool ret = addIndexSensitiveStamp(atomStamps_, atom);
    if (!ret) {
        std::cout << "multiple atoms have the same index: " << atom->getIndex() <<" in " << getName()  << " Molecule\n";
    }
    return ret;
    
}

bool MoleculeStamp::addBondStamp( BondStamp* bond) {
    bondStamps_.push_back(bond);
    return true;
}

bool MoleculeStamp::addBendStamp( BendStamp* bend) {
    bendStamps_.push_back(bend);
    return true;
}

bool MoleculeStamp::addTorsionStamp( TorsionStamp* torsion) {
    torsionStamps_.push_back(torsion);
    return true;
}

bool MoleculeStamp::addRigidBodyStamp( RigidBodyStamp* rigidbody) {
    bool ret = addIndexSensitiveStamp(rigidBodyStamps_, rigidbody);
    if (!ret) {
        std::cout << "multiple rigidbodies have the same index: " << rigidbody->getIndex() <<" in " << getName()  << " Molecule\n";
    }
    return ret;
}

bool MoleculeStamp::addCutoffGroupStamp( CutoffGroupStamp* cutoffgroup) {
    cutoffGroupStamps_.push_back(cutoffgroup);
    return true;
}

bool MoleculeStamp::addFragmentStamp( FragmentStamp* fragment) {
    return addIndexSensitiveStamp(fragmentStamps_, fragment);
}
    
void MoleculeStamp::validate() {
    DataHolder::validate();

    std::vector<AtomStamp*>::iterator ai = std::find(atomStamps_.begin(), atomStamps_.end(), static_cast<AtomStamp*>(NULL));
    if (ai != atomStamps_.end()) {
        std::cout << "Error in Molecule " << getName() << ": atom[" << ai - atomStamps_.begin()<< "] is missing\n";
    }

     std::vector<RigidBodyStamp*>::iterator ri = std::find(rigidBodyStamps_.begin(), rigidBodyStamps_.end(), static_cast<RigidBodyStamp*>(NULL));
     if (ri != rigidBodyStamps_.end()) {
         std::cout << "Error in Molecule " << getName() << ":rigidBody[" <<  ri - rigidBodyStamps_.begin()<< "] is missing\n";
     }
    
    std::vector<FragmentStamp*>::iterator fi = std::find(fragmentStamps_.begin(), fragmentStamps_.end(), static_cast<FragmentStamp*>(NULL));
    if (fi != fragmentStamps_.end()) {
        std::cout << "Error in Molecule " << getName() << ":fragment[" <<  fi - fragmentStamps_.begin()<< "] is missing\n";
    }

    //make sure index is not out of range
    int natoms = getNAtoms();
    for(int i = 0; i < getNBonds(); ++i) {
        BondStamp* bondStamp = getBondStamp(i);
        if (bondStamp->getA() >=  natoms && bondStamp->getB() >= natoms) {
            std::cout << "Error in Molecule " << getName() <<  ": bond between " << bondStamp->getA() << " and " << bondStamp->getB() << " is invalid\n";
        }
    }
    for(int i = 0; i < getNBends(); ++i) {
        BendStamp* bendStamp = getBendStamp(i);
        std::vector<int> bendAtoms =  bendStamp->getMembers();
        std::vector<int>::iterator j =std::find_if(bendAtoms.begin(), bendAtoms.end(), std::bind2nd(std::greater<int>(), natoms-1));
        if (j != bendAtoms.end()) {
            std::cout << "Error in Molecule " << getName();
        }

        if (bendAtoms.size() == 2 && !bendStamp->haveGhostVectorSource()) {
            std::cout << "Error in Molecule " << getName() << ": ghostVectorSouce is missing";
        }
    }    
    for(int i = 0; i < getNBends(); ++i) {
        TorsionStamp* torsionStamp = getTorsionStamp(i);
        std::vector<int> torsionAtoms =  torsionStamp ->getMembers();
        std::vector<int>::iterator j =std::find_if(torsionAtoms.begin(), torsionAtoms.end(), std::bind2nd(std::greater<int>(), natoms-1));
        if (j != torsionAtoms.end()) {
            std::cout << "Error in Molecule " << getName();
        }
    }
    for(int i = 0; i < getNCutoffGroups(); ++i) {
        CutoffGroupStamp* cutoffGroupStamp = getCutoffGroupStamp(i);
        std::vector<int> cutoffGroupAtoms =  cutoffGroupStamp ->getMembers();
        std::vector<int>::iterator j =std::find_if(cutoffGroupAtoms.begin(), cutoffGroupAtoms.end(), std::bind2nd(std::greater<int>(), natoms-1));
        if (j != cutoffGroupAtoms.end()) {
            std::cout << "Error in Molecule " << getName();
        }
    }
        
    atom2Rigidbody.resize(natoms); 
    // negative number means atom is a free atom, does not belong to rigidbody
    //every element in atom2Rigidbody has unique negative number at the very beginning
    for(int i = 0; i < atom2Rigidbody.size(); ++i) {
        atom2Rigidbody[i] = -1 - i;
    }

    for (int i = 0; i < getNRigidBodies(); ++i) {
        RigidBodyStamp* rbStamp = getRigidBodyStamp(i);
        std::vector<int> members = rbStamp->getMembers();
        for(std::vector<int>::iterator j = members.begin(); j != members.end(); ++j) {
            atom2Rigidbody[*j] = i;                 
        }
    }
    //make sure atoms belong to same rigidbody do not bond to each other
    for(int i = 0; i < getNBonds(); ++i) {
        BondStamp* bondStamp = getBondStamp(i);
        if (atom2Rigidbody[bondStamp->getA()] == atom2Rigidbody[bondStamp->getB()])
            std::cout << "Error in Molecule " << getName() << ": "<<"bond between " << bondStamp->getA() << " and " << bondStamp->getB() << "belong to same rigidbody " << atom2Rigidbody[bondStamp->getA()] << "\n";
        }
        
    for(int i = 0; i < getNBends(); ++i) {
        BendStamp* bendStamp = getBendStamp(i);
        std::vector<int> bendAtoms =  bendStamp->getMembers();
        std::vector<int> rigidSet(getNRigidBodies(), 0);
        std::vector<int>::iterator j;
        for( j = bendAtoms.begin(); j != bendAtoms.end(); ++j) {
            int rigidbodyIndex = atom2Rigidbody[*j];
            if (rigidbodyIndex >= 0) {
                ++rigidSet[rigidbodyIndex];
                if (rigidSet[rigidbodyIndex] > 1) {
                    std::cout << "Error in Molecule " << getName() << ": ";
                    //std::cout << "atoms of bend " <<  << "belong to same rigidbody " << rigidbodyIndex << "\n";                    
                }
            }
        }
    }      
    for(int i = 0; i < getNTorsions(); ++i) {
        TorsionStamp* torsionStamp = getTorsionStamp(i);
        std::vector<int> torsionAtoms =  torsionStamp->getMembers();
        std::vector<int> rigidSet(getNRigidBodies(), 0);
        std::vector<int>::iterator j;
        for( j = torsionAtoms.begin(); j != torsionAtoms.end(); ++j) {
            int rigidbodyIndex = atom2Rigidbody[*j];
            if (rigidbodyIndex >= 0) {
                ++rigidSet[rigidbodyIndex];
                if (rigidSet[rigidbodyIndex] > 1) {
                    std::cout << "Error in Molecule " << getName() << ": ";
                    //std::cout << "atoms of torsion " <<  << "belong to same rigidbody " << rigidbodyIndex << "\n";                    
                }
            }
        }
    } 


    //fill in bond information into atom
    fillBondInfo();
    findBends();
    findTorsions();

    int nrigidAtoms = 0;
    for (int i = 0; i < getNRigidBodies(); ++i) {
        RigidBodyStamp* rbStamp = getRigidBodyStamp(i);
        nrigidAtoms += rbStamp->getNMembers();
    }
    nintegrable_ = getNAtoms()+ getNRigidBodies() - nrigidAtoms;

 }

void MoleculeStamp::fillBondInfo() {
}

void MoleculeStamp::findBends() {

}

void MoleculeStamp::findTorsions() {

}

//Function Name: isBondInSameRigidBody
//Return true is both atoms of the bond belong to the same rigid body, otherwise return false
bool MoleculeStamp::isBondInSameRigidBody(BondStamp* bond){
  int rbA;
  int rbB;
  int consAtomA;
  int consAtomB;

  if (!isAtomInRigidBody(bond->getA(),rbA, consAtomA))
    return false;

  if(!isAtomInRigidBody(bond->getB(),rbB, consAtomB) )
    return false;

  if(rbB == rbA)
    return true;
  else
    return false;
}

// Function Name: isAtomInRigidBody 
//return false if atom does not belong to a rigid body, otherwise return true 
bool MoleculeStamp::isAtomInRigidBody(int atomIndex){
  int whichRigidBody;
  int consAtomIndex;

  return isAtomInRigidBody(atomIndex, whichRigidBody, consAtomIndex);
   
}

// Function Name: isAtomInRigidBody 
//return false if atom does not belong to a rigid body otherwise return true and set whichRigidBody 
//and consAtomIndex
//atomIndex : the index of atom in component
//whichRigidBody: the index of rigidbody in component
//consAtomIndex:  the position of joint atom apears in  rigidbody's definition
bool MoleculeStamp::isAtomInRigidBody(int atomIndex, int& whichRigidBody, int& consAtomIndex){
  RigidBodyStamp* rbStamp;
  int numRb;
  int numAtom;

  whichRigidBody = -1;
  consAtomIndex = -1;

  numRb = this->getNRigidBodies();
  
  for(int i = 0 ; i < numRb; i++){
    rbStamp = this->getRigidBodyStamp(i);
    numAtom = rbStamp->getNMembers();
    for(int j = 0; j < numAtom; j++)
      if (rbStamp->getMemberAt(j) == atomIndex){
        whichRigidBody = i;
        consAtomIndex = j;
        return true;
      }
  }

  return false;
   
}

//return the position of joint atom apears in  rigidbody's definition
//for the time being, we will use the most inefficient algorithm, the complexity is O(N2)
//actually we could improve the complexity to O(NlgN) by sorting the atom index in rigid body first
std::vector<std::pair<int, int> > MoleculeStamp::getJointAtoms(int rb1, int rb2){
  RigidBodyStamp* rbStamp1;
  RigidBodyStamp* rbStamp2;
  int natomInRb1;
  int natomInRb2;
  int atomIndex1;
  int atomIndex2;
  std::vector<std::pair<int, int> > jointAtomIndexPair;
  
  rbStamp1 = this->getRigidBodyStamp(rb1);
  natomInRb1 =rbStamp1->getNMembers();

  rbStamp2 = this->getRigidBodyStamp(rb2);
  natomInRb2 =rbStamp2->getNMembers();

  for(int i = 0; i < natomInRb1; i++){
    atomIndex1 = rbStamp1->getMemberAt(i);
      
    for(int j= 0; j < natomInRb1; j++){
      atomIndex2 = rbStamp2->getMemberAt(j);

      if(atomIndex1 == atomIndex2){
        jointAtomIndexPair.push_back(std::make_pair(i, j));
        break;
      }
      
    }

  }

  return jointAtomIndexPair;
}

}
