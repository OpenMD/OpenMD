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

#include <functional> 
#include <iostream>
#include <sstream>
#include "types/MoleculeStamp.hpp"
#include "utils/Tuple.hpp"

namespace oopse {

template<class ContainerType>
std::string containerToString(ContainerType& cont) {
    std::ostringstream oss;
    oss << "(";
    typename ContainerType::iterator i = cont.begin();
    if (i != cont.end()) {
        oss << *i;
        ++i;
    }
    for (; i != cont.end();++i) {
        oss << ", ";
        oss << *i;
    }
    oss << ")";
    return oss.str();
}

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
         std::cout<< "Error in Molecule " << getName()  << ": multiple atoms have the same indices"<< atom->getIndex() <<"\n";
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
        std::cout<< "Error in Molecule " << getName()  << ": multiple rigidbodies have the same indices: " << rigidbody->getIndex() <<"\n";
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

    atom2Rigidbody.resize(getNAtoms()); 
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

    checkAtoms();
    checkBonds();
    fillBondInfo();
    checkBends();
    checkTorsions();
    checkRigidBodies();
    checkCutoffGroups();
    checkFragments();

    int nrigidAtoms = 0;
    for (int i = 0; i < getNRigidBodies(); ++i) {
        RigidBodyStamp* rbStamp = getRigidBodyStamp(i);
        nrigidAtoms += rbStamp->getNMembers();
    }
    nintegrable_ = getNAtoms()+ getNRigidBodies() - nrigidAtoms;

 }

void MoleculeStamp::checkAtoms() {
    std::vector<AtomStamp*>::iterator ai = std::find(atomStamps_.begin(), atomStamps_.end(), static_cast<AtomStamp*>(NULL));
    if (ai != atomStamps_.end()) {
        std::cout << "Error in Molecule " << getName() << ": atom[" << ai - atomStamps_.begin()<< "] is missing\n";
    }

}

void MoleculeStamp::checkBonds() {
    //make sure index is not out of range
    int natoms = getNAtoms();
    for(int i = 0; i < getNBonds(); ++i) {
        BondStamp* bondStamp = getBondStamp(i);
        if (bondStamp->getA() >=  natoms && bondStamp->getB() >= natoms) {
            std::cout << "Error in Molecule " << getName() <<  ": bond(" << bondStamp->getA() << ", " << bondStamp->getB() << ") is invalid\n";
        }
    }
    
    //make sure bonds are unique
    std::set<std::pair<int, int> > allBonds;
    for(int i = 0; i < getNBonds(); ++i) {
        BondStamp* bondStamp= getBondStamp(i);        
        std::pair<int, int> bondPair(bondStamp->getA(), bondStamp->getB());
        //make sure bondTuple.first is always less than or equal to bondTuple.third
        if (bondPair.first > bondPair.second) {
            std::swap(bondPair.first, bondPair.second);
        }
        
        std::set<std::pair<int, int> >::iterator iter = allBonds.find(bondPair);
        if ( iter != allBonds.end()) {
            std::cout << "Error in Molecule " << getName() << ": " << "bond(" <<iter->first << ", "<< iter->second << ")appears multiple times\n";
        } else {
            allBonds.insert(bondPair);
        }
    }
    
    //make sure atoms belong to same rigidbody do not bond to each other
    for(int i = 0; i < getNBonds(); ++i) {
        BondStamp* bondStamp = getBondStamp(i);
        if (atom2Rigidbody[bondStamp->getA()] == atom2Rigidbody[bondStamp->getB()]) {
            std::cout << "Error in Molecule " << getName() << ": "<<"bond(" << bondStamp->getA() << ", " << bondStamp->getB() << ") belong to same rigidbody " << atom2Rigidbody[bondStamp->getA()] << "\n";
        }
    }
    
}

struct BendLessThan : public std::binary_function<IntTuple4, IntTuple4, bool> {
    bool operator()(IntTuple3 b1, IntTuple3 b2) {
        return b1.first < b2.first
             || (!(b2.first < b1.first) && b1.second < b2.second)
             || (!(b2.first < b1.first) && !(b2.second < b2.second) && b1.third < b2.third);
    }
};

void MoleculeStamp::checkBends() {
    for(int i = 0; i < getNBends(); ++i) {
        BendStamp* bendStamp = getBendStamp(i);
        std::vector<int> bendAtoms =  bendStamp->getMembers();
        std::vector<int>::iterator j =std::find_if(bendAtoms.begin(), bendAtoms.end(), std::bind2nd(std::greater<int>(), getNAtoms()-1));
        if (j != bendAtoms.end()) {
            std::cout << "Error in Molecule " << getName() << " : atoms of bend" << containerToString(bendAtoms) << "have invalid indices\n";
        }

        if (bendAtoms.size() == 2 ) {
            if (!bendStamp->haveGhostVectorSource()) {
                std::cout << "Error in Molecule " << getName() << ": ghostVectorSouce is missing\n";
            }else{
                int ghostIndex = bendStamp->getGhostVectorSource();
                if (ghostIndex < getNAtoms()) {
                    if (std::find(bendAtoms.begin(), bendAtoms.end(), ghostIndex) == bendAtoms.end()) {
                      std::cout <<  "Error in Molecule " << getName() << ": ghostVectorSouce "<< ghostIndex<<"is invalid\n"; 
                    }
                    if (!getAtomStamp(ghostIndex)->haveOrientation()) {
                        std::cout <<  "Error in Molecule " << getName() << ": ghost atom must be a directioanl atom\n"; 
                    }
                }else {
                    std::cout << "Error in Molecule " << getName() <<  ": ghostVectorsource " << ghostIndex<< "  is invalid\n";
                }
            }
        } else if (bendAtoms.size() == 3 && bendStamp->haveGhostVectorSource()) {
            std::cout <<  "Error in Molecule " << getName() << ": normal bend should not have ghostVectorSouce\n"; 
        }
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
                    std::cout << "Error in Molecule " << getName() << ": bend" << containerToString(bendAtoms) << " belong to same rigidbody " << rigidbodyIndex << "\n";                    
                }
            }
        }
    } 
    
    
    std::set<IntTuple3, BendLessThan> allBends;
    std::set<IntTuple3, BendLessThan>::iterator iter;
    for(int i = 0; i < getNBends(); ++i) {
        BendStamp* bendStamp= getBendStamp(i);
        std::vector<int> bend = bendStamp->getMembers();
        if (bend.size() == 2) {
        // in case we have two ghost bend. For example, 
        // bend {
        // members (0, 1);
        //   ghostVectorSource = 0;
        // }
        // and
        // bend {
        //   members (0, 1);
        // ghostVectorSource = 0;
        // }
        // In order to distinguish them. we expand them to Tuple3.
        // the first one is expanded to (0, 0, 1) while the second one is expaned to (0, 1, 1)
             int ghostIndex = bendStamp->getGhostVectorSource();
             std::vector<int>::iterator j = std::find(bend.begin(), bend.end(), ghostIndex);
             if (j != bend.end()) {
                bend.insert(j, ghostIndex);
             }
        }
        
        IntTuple3 bendTuple(bend[0], bend[1], bend[2]);
        //make sure bendTuple.first is always less than or equal to bendTuple.third
        if (bendTuple.first > bendTuple.third) {
            std::swap(bendTuple.first, bendTuple.third);
        }
        
        iter = allBends.find(bendTuple);
        if ( iter != allBends.end()) {
            std::cout << "Error in Molecule " << getName() << ": " << "Bend appears multiple times\n";
        } else {
            allBends.insert(bendTuple);
        }
    }

    for (int i = 0; i < getNBonds(); ++i) {
        BondStamp* bondStamp = getBondStamp(i);
        int a = bondStamp->getA();
        int b = bondStamp->getB();

        AtomStamp* atomA = getAtomStamp(a);
        AtomStamp* atomB = getAtomStamp(b);

        //find bend c--a--b
        AtomStamp::AtomIter ai;
        for(int c= atomA->getFirstBonedAtom(ai);c != -1;c = atomA->getNextBonedAtom(ai))
        {
            if(b == c)
                continue;          
            
            IntTuple3 newBend(c, a, b);
            if (newBend.first > newBend.third) {
                std::swap(newBend.first, newBend.third);
            }

            if (allBends.find(newBend) == allBends.end() ) {                
                allBends.insert(newBend);
                BendStamp * newBendStamp = new BendStamp();
                newBendStamp->setMembers(newBend);
                addBendStamp(newBendStamp);
            }
        }        

        //find bend a--b--c
        for(int c= atomB->getFirstBonedAtom(ai);c != -1;c = atomB->getNextBonedAtom(ai))
        {
            if(a == c)
                continue;          

            IntTuple3 newBend( a, b, c);
            if (newBend.first > newBend.third) {
                std::swap(newBend.first, newBend.third);
            }            
            if (allBends.find(newBend) == allBends.end() ) {                
                allBends.insert(newBend);
                BendStamp * newBendStamp = new BendStamp();
                newBendStamp->setMembers(newBend);
                addBendStamp(newBendStamp);
            }
        }        
    }

}

struct TorsionLessThan : public std::binary_function<IntTuple4, IntTuple4, bool> {
    bool operator()(IntTuple4 t1, IntTuple4 t2) {

        return t1.first < t2.first
             || (!(t2.first < t1.first) && t1.second < t2.second)
             || (!(t2.first < t1.first) && !(t2.second < t2.second) && t1.third < t2.third)
             ||(!(t2.first < t1.first) && !(t2.second < t2.second) && !(t2.third < t1.third) && t1.fourth < t2.fourth);
    }



};


void MoleculeStamp::checkTorsions() {
    for(int i = 0; i < getNBends(); ++i) {
        TorsionStamp* torsionStamp = getTorsionStamp(i);
        std::vector<int> torsionAtoms =  torsionStamp ->getMembers();
        std::vector<int>::iterator j =std::find_if(torsionAtoms.begin(), torsionAtoms.end(), std::bind2nd(std::greater<int>(), getNAtoms()-1));
        if (j != torsionAtoms.end()) {
            std::cout << "Error in Molecule " << getName() << ": atoms of torsion" << containerToString(torsionAtoms) << " have invalid indices\n"; 
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
                    std::cout << "Error in Molecule " << getName() << ": torsion" << containerToString(torsionAtoms) << "is invalid\n";                  
                }
            }
        }
    }     

    std::set<IntTuple4, TorsionLessThan> allTorsions;
    std::set<IntTuple4, TorsionLessThan>::iterator iter;
     for(int i = 0; i < getNTorsions(); ++i) {
         TorsionStamp* torsionStamp= getTorsionStamp(i);
         std::vector<int> torsion = torsionStamp->getMembers();
         if (torsion.size() == 3) {
             int ghostIndex = torsionStamp->getGhostVectorSource();
             std::vector<int>::iterator j = std::find(torsion.begin(), torsion.end(), ghostIndex);
             if (j != torsion.end()) {
                torsion.insert(j, ghostIndex);
             }
         }

        IntTuple4 torsionTuple(torsion[0], torsion[1], torsion[2], torsion[3]);
        if (torsionTuple.first > torsionTuple.fourth) {
            std::swap(torsionTuple.first, torsionTuple.fourth);
            std::swap(torsionTuple.second, torsionTuple.third);                    
        }                

         iter = allTorsions.find(torsionTuple);
         if ( iter == allTorsions.end()) {
            allTorsions.insert(torsionTuple);
         } else {
            std::cout << "Error in Molecule " << getName() << ": " << "Torsion appears multiple times\n";
         }
     }

    for (int i = 0; i < getNBonds(); ++i) {
        BondStamp* bondStamp = getBondStamp(i);
        int b = bondStamp->getA();
        int c = bondStamp->getB();

        AtomStamp* atomB = getAtomStamp(b);
        AtomStamp* atomC = getAtomStamp(c);

        AtomStamp::AtomIter ai2;
        AtomStamp::AtomIter ai3;

        for(int a = atomB->getFirstBonedAtom(ai2);a != -1;a = atomB->getNextBonedAtom(ai2))
        {
            if(a == c)
                continue;

            for(int d = atomC->getFirstBonedAtom(ai3);d != -1;d = atomC->getNextBonedAtom(ai3))
            {
                if(d == b)
                    continue;
                
                IntTuple4 newTorsion(a, b, c, d);
                //make sure the first element is always less than or equal to the fourth element in IntTuple4
                if (newTorsion.first > newTorsion.fourth) {
                    std::swap(newTorsion.first, newTorsion.fourth);
                    std::swap(newTorsion.second, newTorsion.third);                    
                }                
                if (allTorsions.find(newTorsion) == allTorsions.end() ) {                
                    allTorsions.insert(newTorsion);
                    TorsionStamp * newTorsionStamp = new TorsionStamp();
                    newTorsionStamp->setMembers(newTorsion);
                    addTorsionStamp(newTorsionStamp);                    
                }            
            }
        }    
    }

}

void MoleculeStamp::checkRigidBodies() {
     std::vector<RigidBodyStamp*>::iterator ri = std::find(rigidBodyStamps_.begin(), rigidBodyStamps_.end(), static_cast<RigidBodyStamp*>(NULL));
     if (ri != rigidBodyStamps_.end()) {
         std::cout << "Error in Molecule " << getName() << ":rigidBody[" <<  ri - rigidBodyStamps_.begin()<< "] is missing\n";
     }

    for (int i = 0; i < getNRigidBodies(); ++i) {
        RigidBodyStamp* rbStamp = getRigidBodyStamp(i);
        std::vector<int> rigidAtoms =  rbStamp ->getMembers();
        std::vector<int>::iterator j =std::find_if(rigidAtoms.begin(), rigidAtoms.end(), std::bind2nd(std::greater<int>(), getNAtoms()-1));
        if (j != rigidAtoms.end()) {
            std::cout << "Error in Molecule " << getName();
        }
        
    }    
}

void MoleculeStamp::checkCutoffGroups() {

    for(int i = 0; i < getNCutoffGroups(); ++i) {
        CutoffGroupStamp* cutoffGroupStamp = getCutoffGroupStamp(i);
        std::vector<int> cutoffGroupAtoms =  cutoffGroupStamp ->getMembers();
        std::vector<int>::iterator j =std::find_if(cutoffGroupAtoms.begin(), cutoffGroupAtoms.end(), std::bind2nd(std::greater<int>(), getNAtoms()-1));
        if (j != cutoffGroupAtoms.end()) {
            std::cout << "Error in Molecule " << getName() << ": cutoffGroup" << " is out of range\n"; 
        }
    }    
}

void MoleculeStamp::checkFragments() {

    std::vector<FragmentStamp*>::iterator fi = std::find(fragmentStamps_.begin(), fragmentStamps_.end(), static_cast<FragmentStamp*>(NULL));
    if (fi != fragmentStamps_.end()) {
        std::cout << "Error in Molecule " << getName() << ":fragment[" <<  fi - fragmentStamps_.begin()<< "] is missing\n";
    }
    
}

void MoleculeStamp::fillBondInfo() {

    for (int i = 0; i < getNBonds(); ++i) {
        BondStamp* bondStamp = getBondStamp(i);
        int a = bondStamp->getA();
        int b = bondStamp->getB();
        AtomStamp* atomA = getAtomStamp(a);
        AtomStamp* atomB = getAtomStamp(b);
        atomA->addBond(i);
        atomA->addBondedAtom(b);
        atomB->addBond(i);        
        atomB->addBondedAtom(a);

    }
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
  return atom2Rigidbody[atomIndex] >=0 ;
   
}

// Function Name: isAtomInRigidBody 
//return false if atom does not belong to a rigid body otherwise return true and set whichRigidBody 
//and consAtomIndex
//atomIndex : the index of atom in component
//whichRigidBody: the index of rigidbody in component
//consAtomIndex:  the position of joint atom apears in  rigidbody's definition
bool MoleculeStamp::isAtomInRigidBody(int atomIndex, int& whichRigidBody, int& consAtomIndex){

  

  whichRigidBody = -1;
  consAtomIndex = -1;

  if (atom2Rigidbody[atomIndex] >=0) {
    whichRigidBody = atom2Rigidbody[atomIndex];
    RigidBodyStamp* rbStamp = getRigidBodyStamp(whichRigidBody);
    int numAtom = rbStamp->getNMembers();
    for(int j = 0; j < numAtom; j++) {
      if (rbStamp->getMemberAt(j) == atomIndex){
        consAtomIndex = j;
        return true;
      }
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
