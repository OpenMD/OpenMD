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
#include <algorithm>
#include <functional> 
#include <iostream>
#include <sstream>
#include "types/MoleculeStamp.hpp"
#include "utils/Tuple.hpp"
#include "utils/MemoryUtils.hpp"
namespace OpenMD {
  
  template<class ContainerType>
  bool hasDuplicateElement(const ContainerType& cont) {
    ContainerType tmp = cont;
    std::sort(tmp.begin(), tmp.end());
    tmp.erase(std::unique(tmp.begin(), tmp.end()), tmp.end());
    return tmp.size() != cont.size();
  }
  
  MoleculeStamp::MoleculeStamp() {
    DefineParameter(Name, "name");
    
    deprecatedKeywords_.insert("nAtoms");
    deprecatedKeywords_.insert("nBonds");
    deprecatedKeywords_.insert("nBends");
    deprecatedKeywords_.insert("nTorsions");
    deprecatedKeywords_.insert("nInversions");
    deprecatedKeywords_.insert("nRigidBodies");
    deprecatedKeywords_.insert("nCutoffGroups");
    
  }
  
  MoleculeStamp::~MoleculeStamp() {
    MemoryUtils::deletePointers(atomStamps_);
    MemoryUtils::deletePointers(bondStamps_);
    MemoryUtils::deletePointers(bendStamps_);
    MemoryUtils::deletePointers(torsionStamps_);
    MemoryUtils::deletePointers(inversionStamps_);
    MemoryUtils::deletePointers(rigidBodyStamps_);
    MemoryUtils::deletePointers(cutoffGroupStamps_);
    MemoryUtils::deletePointers(fragmentStamps_);    
  }
  
  bool MoleculeStamp::addAtomStamp( AtomStamp* atom) {
    bool ret = addIndexSensitiveStamp(atomStamps_, atom);
    if (!ret) {
      std::ostringstream oss;
      oss<< "Error in Molecule " << getName()  << 
        ": multiple atoms have the same indices"<< atom->getIndex() <<"\n";
      throw OpenMDException(oss.str());
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

  bool MoleculeStamp::addInversionStamp( InversionStamp* inversion) {
    inversionStamps_.push_back(inversion);
    return true;
  }
  
  bool MoleculeStamp::addRigidBodyStamp( RigidBodyStamp* rigidbody) {
    bool ret = addIndexSensitiveStamp(rigidBodyStamps_, rigidbody);
    if (!ret) {
      std::ostringstream oss;
      oss << "Error in Molecule " << getName()  << 
        ": multiple rigidbodies have the same indices: " << 
        rigidbody->getIndex() <<"\n";
      throw OpenMDException(oss.str());
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

    // A negative number means the atom is a free atom, and does not
    // belong to rigidbody. Every element in atom2Rigidbody has unique
    // negative number at the very beginning

    for(int i = 0; i < atom2Rigidbody.size(); ++i) {
      atom2Rigidbody[i] = -1 - i;
    }
    for (int i = 0; i < getNRigidBodies(); ++i) {
      RigidBodyStamp* rbStamp = getRigidBodyStamp(i);
      std::vector<int> members = rbStamp->getMembers();
      for(std::vector<int>::iterator j = members.begin(); 
          j != members.end(); ++j) {
        atom2Rigidbody[*j] = i;                 
      }
    }
    
    checkAtoms();
    checkBonds();
    fillBondInfo();
    checkBends();
    checkTorsions();
    checkInversions();
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
    std::vector<AtomStamp*>::iterator ai = std::find(atomStamps_.begin(), 
                                                     atomStamps_.end(), 
                                                     static_cast<AtomStamp*>(NULL));
    if (ai != atomStamps_.end()) {
      std::ostringstream oss;
      oss << "Error in Molecule " << getName() << ": atom[" << 
        ai - atomStamps_.begin()<< "] is missing\n";
      throw OpenMDException(oss.str());
    }
    
  }

  void MoleculeStamp::checkBonds() {
    std::ostringstream oss;
    //make sure index is not out of range
    int natoms = getNAtoms();
    for(int i = 0; i < getNBonds(); ++i) {
      BondStamp* bondStamp = getBondStamp(i);
      if (bondStamp->getA() > natoms-1 ||  bondStamp->getA() < 0 || 
          bondStamp->getB() > natoms-1 || bondStamp->getB() < 0 || 
          bondStamp->getA() == bondStamp->getB()) {
        
        oss << "Error in Molecule " << getName() <<  ": bond(" << 
          bondStamp->getA() << ", " << bondStamp->getB() << ") is invalid\n";
        throw OpenMDException(oss.str());
      }
    }
    
    //make sure bonds are unique
    std::set<std::pair<int, int> > allBonds;
    for(int i = 0; i < getNBonds(); ++i) {
      BondStamp* bondStamp= getBondStamp(i);        
      std::pair<int, int> bondPair(bondStamp->getA(), bondStamp->getB());
      //make sure bondTuple.first is always less than or equal to
      //bondTuple.third
      if (bondPair.first > bondPair.second) {
        std::swap(bondPair.first, bondPair.second);
      }
      
      std::set<std::pair<int, int> >::iterator iter = allBonds.find(bondPair);
      if ( iter != allBonds.end()) {
        
        oss << "Error in Molecule " << getName() << ": " << "bond(" <<
          iter->first << ", "<< iter->second << ") appears multiple times\n";
        throw OpenMDException(oss.str());
      } else {
        allBonds.insert(bondPair);
      }
    }
    
    //make sure atoms belong to same rigidbody do not bond to each other
    for(int i = 0; i < getNBonds(); ++i) {
      BondStamp* bondStamp = getBondStamp(i);
      if (atom2Rigidbody[bondStamp->getA()] == atom2Rigidbody[bondStamp->getB()]) {
        
        oss << "Error in Molecule " << getName() << ": "<<"bond(" << 
          bondStamp->getA() << ", " << bondStamp->getB() << 
          ") belong to same rigidbody " << 
          atom2Rigidbody[bondStamp->getA()] << "\n";
        throw OpenMDException(oss.str());
      }
    }    
  }
  
  void MoleculeStamp::checkBends() {
    std::ostringstream oss;
    for(int i = 0; i < getNBends(); ++i) {
      BendStamp* bendStamp = getBendStamp(i);
      std::vector<int> bendAtoms =  bendStamp->getMembers();
      std::vector<int>::iterator j =std::find_if(bendAtoms.begin(), 
                                                 bendAtoms.end(), 
                                                 std::bind2nd(std::greater<int>(), 
                                                              getNAtoms()-1));
      std::vector<int>::iterator k =std::find_if(bendAtoms.begin(), 
                                                 bendAtoms.end(), 
                                                 std::bind2nd(std::less<int>(),
                                                              0));
      
      if (j != bendAtoms.end() || k != bendAtoms.end()) {
        
        oss << "Error in Molecule " << getName() << " : atoms of bend" << 
          containerToString(bendAtoms) << " have invalid indices\n";
        throw OpenMDException(oss.str());
      }
      
      if (hasDuplicateElement(bendAtoms)) {
        oss << "Error in Molecule " << getName() << " : atoms of bend" << 
          containerToString(bendAtoms) << " have duplicated indices\n";    
        throw OpenMDException(oss.str());            
      }
      
      if (bendAtoms.size() == 2 ) {
        if (!bendStamp->haveGhostVectorSource()) {
          
          oss << "Error in Molecule " << getName() << 
            ": ghostVectorSouce is missing\n";
          throw OpenMDException(oss.str());
        }else{
          int ghostIndex = bendStamp->getGhostVectorSource();
          if (ghostIndex < getNAtoms()) {
            if (std::find(bendAtoms.begin(), bendAtoms.end(), 
                          ghostIndex) == bendAtoms.end()) {
              
              oss <<  "Error in Molecule " << getName() << 
                ": ghostVectorSouce "<< ghostIndex<<"is invalid\n"; 
              throw OpenMDException(oss.str());
            }
            if (!getAtomStamp(ghostIndex)->haveOrientation()) {
              
              oss <<  "Error in Molecule " << getName() << 
                ": ghost atom must be a directioanl atom\n"; 
              throw OpenMDException(oss.str());
            }
          } else {
            oss << "Error in Molecule " << getName() <<  
              ": ghostVectorSource " << ghostIndex<< "  is invalid\n";
            throw OpenMDException(oss.str());
          }
        }
      } else if (bendAtoms.size() == 3 && bendStamp->haveGhostVectorSource()) {
        oss <<  "Error in Molecule " << getName() << 
          ": normal bend should not have ghostVectorSouce\n"; 
        throw OpenMDException(oss.str());
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
            oss << "Error in Molecule " << getName() << ": bend" << 
              containerToString(bendAtoms) << " belong to same rigidbody " << 
              rigidbodyIndex << "\n";  
            throw OpenMDException(oss.str());
          }
        }
      }
    } 
    
    
    std::set<IntTuple3> allBends;
    std::set<IntTuple3>::iterator iter;
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
        // the first one is expanded to (0, 0, 1) while the second one
        // is expaned to (0, 1, 1)
        int ghostIndex = bendStamp->getGhostVectorSource();
        std::vector<int>::iterator j = std::find(bend.begin(), bend.end(), 
                                                 ghostIndex);
        if (j != bend.end()) {
          bend.insert(j, ghostIndex);
        }
      }
      
      IntTuple3 bendTuple(bend[0], bend[1], bend[2]);

      // make sure bendTuple.first is always less than or equal to
      // bendTuple.third
      if (bendTuple.first > bendTuple.third) {
        std::swap(bendTuple.first, bendTuple.third);
      }
      
      iter = allBends.find(bendTuple);
      if ( iter != allBends.end()) {
        oss << "Error in Molecule " << getName() << ": " << "Bend" << 
          containerToString(bend)<< " appears multiple times\n";
        throw OpenMDException(oss.str());
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
      for(int c= atomA->getFirstBondedAtom(ai);c != -1;
          c = atomA->getNextBondedAtom(ai)) {
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
      for(int c= atomB->getFirstBondedAtom(ai);c != -1;
          c = atomB->getNextBondedAtom(ai)) {
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
  
  void MoleculeStamp::checkTorsions() {
    std::ostringstream oss;
    for(int i = 0; i < getNTorsions(); ++i) {
      TorsionStamp* torsionStamp = getTorsionStamp(i);
      std::vector<int> torsionAtoms =  torsionStamp ->getMembers();
      std::vector<int>::iterator j =std::find_if(torsionAtoms.begin(), 
                                                 torsionAtoms.end(), 
                                                 std::bind2nd(std::greater<int>(), 
                                                              getNAtoms()-1));
      std::vector<int>::iterator k =std::find_if(torsionAtoms.begin(), 
                                                 torsionAtoms.end(), 
                                                 std::bind2nd(std::less<int>(), 0));
      
      if (j != torsionAtoms.end() || k != torsionAtoms.end()) {
        oss << "Error in Molecule " << getName() << ": atoms of torsion" << 
          containerToString(torsionAtoms) << " have invalid indices\n"; 
        throw OpenMDException(oss.str());
      }
      if (hasDuplicateElement(torsionAtoms)) {
        oss << "Error in Molecule " << getName() << " : atoms of torsion" << 
          containerToString(torsionAtoms) << " have duplicated indices\n";    
        throw OpenMDException(oss.str());            
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
            oss << "Error in Molecule " << getName() << ": torsion" << 
              containerToString(torsionAtoms) << "is invalid\n";           
            throw OpenMDException(oss.str());
          }
        }
      }
    }     
    
    std::set<IntTuple4> allTorsions;
    std::set<IntTuple4>::iterator iter;
    for(int i = 0; i < getNTorsions(); ++i) {
      TorsionStamp* torsionStamp= getTorsionStamp(i);
      std::vector<int> torsion = torsionStamp->getMembers();
      if (torsion.size() == 3) {
        int ghostIndex = torsionStamp->getGhostVectorSource();
        std::vector<int>::iterator j = std::find(torsion.begin(), 
                                                 torsion.end(), ghostIndex);
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
        oss << "Error in Molecule " << getName() << ": " << "Torsion" << 
          containerToString(torsion)<< " appears multiple times\n";
        throw OpenMDException(oss.str());
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
      
      for(int a = atomB->getFirstBondedAtom(ai2);a != -1;
          a = atomB->getNextBondedAtom(ai2)) {
        if(a == c)
          continue;
        
        for(int d = atomC->getFirstBondedAtom(ai3);d != -1;
            d = atomC->getNextBondedAtom(ai3)) {          
          if(d == b)
            continue;
          
          IntTuple4 newTorsion(a, b, c, d);
          //make sure the first element is always less than or equal
          //to the fourth element in IntTuple4
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

  void MoleculeStamp::checkInversions() {
    std::ostringstream oss;

    // first we automatically find the other three atoms that 
    // are satellites of an inversion center:

    for(int i = 0; i < getNInversions(); ++i) {
      InversionStamp* inversionStamp = getInversionStamp(i);
      int center = inversionStamp->getCenter();
      std::vector<int> satellites;

      for (int j = 0; j < getNBonds(); ++j) {
	BondStamp* bondStamp = getBondStamp(j);
	int a = bondStamp->getA();
	int b = bondStamp->getB();
      
	if(a == center) {
	  satellites.push_back(b);
	}
	if (b == center) {
	  satellites.push_back(a);
	}
      }

      if (satellites.size() == 3) {
	std::sort(satellites.begin(), satellites.end());
	inversionStamp->setSatellites(satellites);	
      } else {
	oss << "Error in Molecule " << getName() << ": found wrong number" << 
	  " of bonds for inversion center " << center;
        throw OpenMDException(oss.str());
      }
    }
    

    // then we do some sanity checking on the inversions:
    
    for(int i = 0; i < getNInversions(); ++i) {
      InversionStamp* inversionStamp = getInversionStamp(i);

      std::vector<int> inversionAtoms =  inversionStamp->getSatellites();
      // add the central atom to the beginning of the list:
      inversionAtoms.insert(inversionAtoms.begin(),inversionStamp->getCenter());

      std::vector<int>::iterator j =std::find_if(inversionAtoms.begin(), 
                                                 inversionAtoms.end(), 
                                                 std::bind2nd(std::greater<int>(), 
                                                              getNAtoms()-1));
      std::vector<int>::iterator k =std::find_if(inversionAtoms.begin(), 
                                                 inversionAtoms.end(), 
                                                 std::bind2nd(std::less<int>(), 0));
      
      if (j != inversionAtoms.end() || k != inversionAtoms.end()) {
        oss << "Error in Molecule " << getName() << ": atoms of inversion" << 
          containerToString(inversionAtoms) << " have invalid indices\n"; 
        throw OpenMDException(oss.str());
      }


      if (hasDuplicateElement(inversionAtoms)) {
        oss << "Error in Molecule " << getName() << " : atoms of inversion" << 
          containerToString(inversionAtoms) << " have duplicated indices\n";
        throw OpenMDException(oss.str());            
      }        
    }
    


    for(int i = 0; i < getNInversions(); ++i) {
      InversionStamp* inversionStamp = getInversionStamp(i);
      std::vector<int> inversionAtoms =  inversionStamp->getSatellites();
      inversionAtoms.push_back(inversionStamp->getCenter());
      std::vector<int> rigidSet(getNRigidBodies(), 0);
      std::vector<int>::iterator j;
      for( j = inversionAtoms.begin(); j != inversionAtoms.end(); ++j) {
        int rigidbodyIndex = atom2Rigidbody[*j];
        if (rigidbodyIndex >= 0) {
          ++rigidSet[rigidbodyIndex];
          if (rigidSet[rigidbodyIndex] > 1) {
            oss << "Error in Molecule " << getName() << ": inversion centered on atom " << 
              inversionStamp->getCenter() << " has atoms that belong to same rigidbody " << 
              rigidbodyIndex << "\n";  
            throw OpenMDException(oss.str());
          }
        }
      }
    } 

    std::set<IntTuple4> allInversions;
    std::set<IntTuple4>::iterator iter;
    for(int i = 0; i < getNInversions(); ++i) {
      InversionStamp* inversionStamp= getInversionStamp(i);
      int cent = inversionStamp->getCenter();
      std::vector<int> inversion = inversionStamp->getSatellites();
      
      IntTuple4 inversionTuple(cent, inversion[0], inversion[1], inversion[2]);

      // In OpenMD, the Central atom in an inversion comes first, and
      // has a special position.  The other three atoms can come in
      // random order, and should be sorted in increasing numerical
      // order to check for duplicates.  This requires three pairwise
      // swaps:

      if (inversionTuple.third > inversionTuple.fourth) 
	std::swap(inversionTuple.third, inversionTuple.fourth);

      if (inversionTuple.second > inversionTuple.third) 
	std::swap(inversionTuple.second, inversionTuple.third);

      if (inversionTuple.third > inversionTuple.fourth) 
	std::swap(inversionTuple.third, inversionTuple.fourth);

      
      iter = allInversions.find(inversionTuple);
      if ( iter == allInversions.end()) {
        allInversions.insert(inversionTuple);
      } else {
        oss << "Error in Molecule " << getName() << ": " << "Inversion" << 
          containerToString(inversion)<< " appears multiple times\n";
        throw OpenMDException(oss.str());
      }
    }


    // Next we automatically find the inversion centers that weren't
    // explicitly created.  An inversion center is any atom that has
    // exactly three bonds to it.  Not all inversion centers have
    // potentials associated with them.    
    
    for (int i = 0; i < getNAtoms(); ++i) {
      AtomStamp* ai = getAtomStamp(i);
      if (ai->getBondCount() == 3) {
        AtomStamp::AtomIter ai2;
        std::vector<int> satellites;
        for(int a = ai->getFirstBondedAtom(ai2);a != -1;
            a = ai->getNextBondedAtom(ai2)) {
          satellites.push_back(a);
        }
        if (satellites.size() == 3) {
          int cent = ai->getIndex();
          std::sort(satellites.begin(), satellites.end());
          IntTuple4 newInversion(cent, satellites[0], satellites[1], satellites[2]);

          if (newInversion.third > newInversion.fourth) 
            std::swap(newInversion.third, newInversion.fourth);
          
          if (newInversion.second > newInversion.third) 
            std::swap(newInversion.second, newInversion.third);
          
          if (newInversion.third > newInversion.fourth) 
            std::swap(newInversion.third, newInversion.fourth);          

          if (allInversions.find(newInversion) == allInversions.end() ) {
            allInversions.insert(newInversion);
            InversionStamp * newInversionStamp = new InversionStamp();
            newInversionStamp->setCenter(cent);
            newInversionStamp->setSatellites(satellites);
            addInversionStamp(newInversionStamp);        
          }            

        } else {
          oss << "Error in Molecule " << getName() << ": found bond mismatch" << 
            " when detecting inversion centers.";
          throw OpenMDException(oss.str());
        }
        
      }
    }
  }
  
  void MoleculeStamp::checkRigidBodies() {
    std::ostringstream oss;
    std::vector<RigidBodyStamp*>::iterator ri = std::find(rigidBodyStamps_.begin(), 
                                                          rigidBodyStamps_.end(), 
                                                          static_cast<RigidBodyStamp*>(NULL));
    if (ri != rigidBodyStamps_.end()) {
      oss << "Error in Molecule " << getName() << ":rigidBody[" <<  
        ri - rigidBodyStamps_.begin()<< "] is missing\n";
      throw OpenMDException(oss.str());
    }
    
    for (int i = 0; i < getNRigidBodies(); ++i) {
      RigidBodyStamp* rbStamp = getRigidBodyStamp(i);
      std::vector<int> rigidAtoms =  rbStamp ->getMembers();
      std::vector<int>::iterator j =std::find_if(rigidAtoms.begin(), 
                                                 rigidAtoms.end(), 
                                                 std::bind2nd(std::greater<int>(), 
                                                              getNAtoms()-1));
      if (j != rigidAtoms.end()) {
        oss << "Error in Molecule " << getName();
        throw OpenMDException(oss.str());
      }      
    }    
  }
  
  void MoleculeStamp::checkCutoffGroups() {
    std::vector<AtomStamp*>::iterator ai;
    std::vector<int>::iterator fai;

    //add all atoms into freeAtoms_ set
    for(ai = atomStamps_.begin(); ai != atomStamps_.end(); ++ai) {
      freeAtoms_.push_back( (*ai)->getIndex() );
    }

    for(int i = 0; i < getNCutoffGroups(); ++i) {
      CutoffGroupStamp* cutoffGroupStamp = getCutoffGroupStamp(i);
      std::vector<int> cutoffGroupAtoms =  cutoffGroupStamp ->getMembers();
      std::vector<int>::iterator j =std::find_if(cutoffGroupAtoms.begin(), 
                                                 cutoffGroupAtoms.end(), 
                                                 std::bind2nd(std::greater<int>(), 
                                                              getNAtoms()-1));
      if (j != cutoffGroupAtoms.end()) {
        std::ostringstream oss;
        oss << "Error in Molecule " << getName() << ": cutoffGroup" << 
          " is out of range\n"; 
        throw OpenMDException(oss.str());
      }

      for(fai = cutoffGroupAtoms.begin(); fai != cutoffGroupAtoms.end(); 
          ++fai) {

        // erase the atoms belonging to cutoff groups from freeAtoms_ vector        
        freeAtoms_.erase(std::remove(freeAtoms_.begin(),
                                     freeAtoms_.end(), (*fai)), 
                         freeAtoms_.end());
      }
    }  
  }

  void MoleculeStamp::checkFragments() {

    std::vector<FragmentStamp*>::iterator fi = std::find(fragmentStamps_.begin(), 
                                                         fragmentStamps_.end(),
                                                         static_cast<FragmentStamp*>(NULL));
    if (fi != fragmentStamps_.end()) {
      std::ostringstream oss;
      oss << "Error in Molecule " << getName() << ":fragment[" <<  
        fi - fragmentStamps_.begin()<< "] is missing\n";
      throw OpenMDException(oss.str());
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



  // Function Name: isBondInSameRigidBody
  // Returns true is both atoms of the bond belong to the same rigid
  // body, otherwise return false
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
  // Returns false if atom does not belong to a rigid body, otherwise
  // returns true
  bool MoleculeStamp::isAtomInRigidBody(int atomIndex){
    return atom2Rigidbody[atomIndex] >=0 ;    
  }
  
  // Function Name: isAtomInRigidBody 
  // Returns false if atom does not belong to a rigid body otherwise
  // returns true and sets whichRigidBody and consAtomIndex
  // atomIndex : the index of atom in component
  // whichRigidBody: the index of the rigidbody in the component
  // consAtomIndex:  the position the joint atom appears in the rigidbody's 
  //                 definition 
  bool MoleculeStamp::isAtomInRigidBody(int atomIndex, int& whichRigidBody, 
                                        int& consAtomIndex){
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
  
  // Returns the position of joint atoms apearing in a rigidbody's definition
  // For the time being, we will use the most inefficient algorithm,
  // the complexity is O(N^2).  We could improve the
  // complexity to O(NlogN) by sorting the atom index in rigid body
  // first
  std::vector<std::pair<int, int> > MoleculeStamp::getJointAtoms(int rb1, 
                                                                 int rb2){
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
