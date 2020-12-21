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
 
/**
 * @file Molecule.cpp
 * @author    tlin
 * @date  10/28/2004
 * @version 1.0
 */ 

#include <algorithm>
#include <memory>
#include <set>

#include "primitives/Molecule.hpp"
#include "utils/MemoryUtils.hpp"
#include "utils/simError.h"
#include "utils/StringUtils.hpp"

namespace OpenMD {
  Molecule::Molecule(int stampId, int globalIndex, const std::string& molName, 
                     int region) : 
                                   globalIndex_(globalIndex), 
                                   stampId_(stampId),
				   region_(region),
				   moleculeName_(molName),
                                   constrainTotalCharge_(false) {
  }
  
  Molecule::~Molecule() {
    
    MemoryUtils::deletePointers(atoms_);
    MemoryUtils::deletePointers(bonds_);
    MemoryUtils::deletePointers(bends_);
    MemoryUtils::deletePointers(torsions_);
    MemoryUtils::deletePointers(inversions_);
    MemoryUtils::deletePointers(rigidBodies_);
    MemoryUtils::deletePointers(cutoffGroups_);
    MemoryUtils::deletePointers(constraintPairs_);
    MemoryUtils::deletePointers(constraintElems_);
    MemoryUtils::deletePointers(hBondDonors_);

    // integrableObjects_ don't own the objects
    integrableObjects_.clear();
    fluctuatingCharges_.clear();
  }
  
  void Molecule::addAtom(Atom* atom) {
    if (std::find(atoms_.begin(), atoms_.end(), atom) == atoms_.end()) {
      atoms_.push_back(atom);
    }
  }
  
  void Molecule::addBond(Bond* bond) {
    if (std::find(bonds_.begin(), bonds_.end(), bond) == bonds_.end()) {
      bonds_.push_back(bond);
    }
  }
  
  void Molecule::addBend(Bend* bend) {
    if (std::find(bends_.begin(), bends_.end(), bend) == bends_.end()) {
      bends_.push_back(bend);
    }
  }
  
  void Molecule::addTorsion(Torsion* torsion) {
    if (std::find(torsions_.begin(), torsions_.end(), torsion) == 
        torsions_.end()) {
      torsions_.push_back(torsion);
    }
  }

  void Molecule::addInversion(Inversion* inversion) {
    if (std::find(inversions_.begin(), inversions_.end(), inversion) == 
        inversions_.end()) {
      inversions_.push_back(inversion);
    }
  }
  
  void Molecule::addRigidBody(RigidBody *rb) {
    if (std::find(rigidBodies_.begin(), rigidBodies_.end(), rb) == 
        rigidBodies_.end()) {
      rigidBodies_.push_back(rb);
    }
  }
  
  void Molecule::addCutoffGroup(CutoffGroup* cp) {
    if (std::find(cutoffGroups_.begin(), cutoffGroups_.end(), cp) == 
        cutoffGroups_.end()) {
      cutoffGroups_.push_back(cp);
    }    
  }
  
  void Molecule::addConstraintPair(ConstraintPair* cp) {
    if (std::find(constraintPairs_.begin(), constraintPairs_.end(), cp) == 
        constraintPairs_.end()) {
      constraintPairs_.push_back(cp);
    }    
  }
  
  void Molecule::addConstraintElem(ConstraintElem* cp) {
    if (std::find(constraintElems_.begin(), constraintElems_.end(), cp) == 
        constraintElems_.end()) {
      constraintElems_.push_back(cp);
    }
  }
  
  void Molecule::complete() {
    
    std::set<Atom*> rigidAtoms;
    Atom* atom;
    Atom* atom1;
    Atom* atom2;
    AtomIterator ai, aj;
    RigidBody* rb;
    RigidBodyIterator rbIter;
    Bond* bond;
    BondIterator bi;

    // Get list of all the atoms that are part of rigid bodies

    for (rb = beginRigidBody(rbIter); rb != NULL; rb = nextRigidBody(rbIter)) {
      rigidAtoms.insert(rb->getBeginAtomIter(), rb->getEndAtomIter());
    }
    
    // add any atom that wasn't part of a rigid body to the list of integrableObjects

    for (atom = beginAtom(ai); atom != NULL; atom = nextAtom(ai)) {
      
      if (rigidAtoms.find(atom) == rigidAtoms.end()) {

	// If an atom does not belong to a rigid body, it is an
	// integrable object

	integrableObjects_.push_back(atom);
      }
    }
    
    // then add the rigid bodies themselves to the integrableObjects

    for (rb = beginRigidBody(rbIter); rb != NULL; rb = nextRigidBody(rbIter)) {
      integrableObjects_.push_back(rb);
    } 

    // find the atoms that are fluctuating charges and add them to the 
    // fluctuatingCharges_ vector

    for (atom = beginAtom(ai); atom != NULL; atom = nextAtom(ai)) {
      if ( atom->isFluctuatingCharge() )
        fluctuatingCharges_.push_back( atom );      
    }


    // find the electronegative atoms and add them to the
    // hBondAcceptors_ vector:
    
    for (atom = beginAtom(ai); atom != NULL; atom = nextAtom(ai)) {
      AtomType* at = atom->getAtomType();
      // get the chain of base types for this atom type:
      std::vector<AtomType*> ayb = at->allYourBase();
      // use the last type in the chain of base types for the name:
      std::string bn = UpperCase(ayb[ayb.size()-1]->getName());

        if (bn.compare("O")==0 || bn.compare("N")==0
            || bn.compare("F")==0) 
          hBondAcceptors_.push_back( atom );
      
    }
    
    // find electronegative atoms that are either bonded to
    // hydrogens or are present in the same rigid bodies:
    
    for (bond = beginBond(bi); bond != NULL; bond = nextBond(bi)) {
      Atom* atom1 = bond->getAtomA();
      Atom* atom2 = bond->getAtomB();
      AtomType* at1 = atom1->getAtomType();
      AtomType* at2 = atom1->getAtomType();
      // get the chain of base types for this atom type:
      std::vector<AtomType*> ayb1 = at1->allYourBase();
      std::vector<AtomType*> ayb2 = at2->allYourBase();
      // use the last type in the chain of base types for the name:
      std::string bn1 = UpperCase(ayb1[ayb1.size()-1]->getName());
      std::string bn2 = UpperCase(ayb2[ayb2.size()-1]->getName());
      
      if (bn1.compare("H")==0) {
        if (bn2.compare("O")==0 || bn2.compare("N")==0
            || bn2.compare("F")==0) {
          HBondDonor* donor = new HBondDonor();
          donor->donorAtom = atom2;
          donor->donatedHydrogen = atom1;
          hBondDonors_.push_back( donor );
        }
      }
      if (bn2.compare("H")==0) {
        if (bn1.compare("O")==0 || bn1.compare("N")==0
            || bn1.compare("F")==0) {
          HBondDonor* donor = new HBondDonor();
          donor->donorAtom = atom1;
          donor->donatedHydrogen = atom2;
            hBondDonors_.push_back( donor );
        }
      } 
    }
    
    for (rb = beginRigidBody(rbIter); rb != NULL;
         rb = nextRigidBody(rbIter)) {
      for(atom1 = rb->beginAtom(ai); atom1 != NULL;
          atom1 = rb->nextAtom(ai)) {
        AtomType* at1 = atom1->getAtomType();
        // get the chain of base types for this atom type:
        std::vector<AtomType*> ayb1 = at1->allYourBase();
        // use the last type in the chain of base types for the name:
        std::string bn1 = UpperCase(ayb1[ayb1.size()-1]->getName());
        
        if (bn1.compare("O")==0 || bn1.compare("N")==0
            || bn1.compare("F")==0) {
          for(atom2 = rb->beginAtom(aj); atom2 != NULL;
              atom2 = rb->nextAtom(aj)) {
            AtomType* at2 = atom2->getAtomType();
            // get the chain of base types for this atom type:              
            std::vector<AtomType*> ayb2 = at2->allYourBase();
            // use the last type in the chain of base types for the name:
            std::string bn2 = UpperCase(ayb2[ayb2.size()-1]->getName());
            if (bn2.compare("H")==0) {              
              HBondDonor* donor = new HBondDonor();
              donor->donorAtom = atom1;
              donor->donatedHydrogen = atom2;
              hBondDonors_.push_back( donor );
            }
          }
        }
      }
    }
  }
  
  RealType Molecule::getMass() {
    StuntDouble* sd;
    std::vector<StuntDouble*>::iterator i;
    RealType mass = 0.0;
    
    for (sd = beginIntegrableObject(i); sd != NULL; sd = 
           nextIntegrableObject(i)){
      mass += sd->getMass();
    }
    
    return mass;    
  }

  Vector3d Molecule::getCom() {
    StuntDouble* sd;
    std::vector<StuntDouble*>::iterator i;
    Vector3d com;
    RealType totalMass = 0;
    RealType mass;
    
    for (sd = beginIntegrableObject(i); sd != NULL; sd = 
           nextIntegrableObject(i)){
      mass = sd->getMass();
      totalMass += mass;
      com += sd->getPos() * mass;    
    }
    
    com /= totalMass;

    return com;
  }
  
  Vector3d Molecule::getCom(int snapshotNo) {
    StuntDouble* sd;
    std::vector<StuntDouble*>::iterator i;
    Vector3d com;
    RealType totalMass = 0;
    RealType mass;
    
    for (sd = beginIntegrableObject(i); sd != NULL; sd = 
           nextIntegrableObject(i)){
      mass = sd->getMass();
      totalMass += mass;
      com += sd->getPos(snapshotNo) * mass;    
    }
    
    com /= totalMass;

    return com;
  }

  void Molecule::moveCom(const Vector3d& delta) {
    StuntDouble* sd;
    std::vector<StuntDouble*>::iterator i;
    
    for (sd = beginIntegrableObject(i); sd != NULL; sd = 
           nextIntegrableObject(i)){
      sd->setPos(sd->getPos() + delta);
    }    
  }
  
  Vector3d Molecule::getComVel() {
    StuntDouble* sd;
    std::vector<StuntDouble*>::iterator i;
    Vector3d velCom;
    RealType totalMass = 0;
    RealType mass;
    
    for (sd = beginIntegrableObject(i); sd != NULL; sd = 
           nextIntegrableObject(i)){
      mass = sd->getMass();
      totalMass += mass;
      velCom += sd->getVel() * mass;    
    }
    
    velCom /= totalMass;
    
    return velCom;
  }

  RealType Molecule::getPotential() {

    Bond* bond;
    Bend* bend;
    Torsion* torsion;
    Inversion* inversion;
    Molecule::BondIterator bondIter;;
    Molecule::BendIterator  bendIter;
    Molecule::TorsionIterator  torsionIter;
    Molecule::InversionIterator  inversionIter;

    RealType potential = 0.0;

    for (bond = beginBond(bondIter); bond != NULL; bond = nextBond(bondIter)) {
      potential += bond->getPotential();
    }

    for (bend = beginBend(bendIter); bend != NULL; bend = nextBend(bendIter)) {
      potential += bend->getPotential();
    }

    for (torsion = beginTorsion(torsionIter); torsion != NULL; torsion = 
           nextTorsion(torsionIter)) {
      potential += torsion->getPotential();
    }
    
    for (inversion = beginInversion(inversionIter); inversion != NULL; 
         inversion =  nextInversion(inversionIter)) {
      potential += inversion->getPotential();
    }
    
    return potential;
    
  }
  
  void Molecule::addProperty(std::shared_ptr<GenericData> genData) {
    properties_.addProperty(genData);  
  }

  void Molecule::removeProperty(const std::string& propName) {
    properties_.removeProperty(propName);  
  }

  std::vector<std::string> Molecule::getPropertyNames() {
    return properties_.getPropertyNames();  
  }
      
  std::vector<std::shared_ptr<GenericData> > Molecule::getProperties() { 
    return properties_.getProperties(); 
  }

  std::shared_ptr<GenericData> Molecule::getPropertyByName(const std::string& propName) {
    return properties_.getPropertyByName(propName); 
  }

  std::ostream& operator <<(std::ostream& o, Molecule& mol) {
    o << std::endl;
    o << "Molecule " << mol.getGlobalIndex() << "has: " << std::endl;
    o << mol.getNAtoms() << " atoms" << std::endl;
    o << mol.getNBonds() << " bonds" << std::endl;
    o << mol.getNBends() << " bends" << std::endl;
    o << mol.getNTorsions() << " torsions" << std::endl;
    o << mol.getNInversions() << " inversions" << std::endl;
    o << mol.getNRigidBodies() << " rigid bodies" << std::endl;
    o << mol.getNIntegrableObjects() << " integrable objects" << std::endl;
    o << mol.getNCutoffGroups() << " cutoff groups" << std::endl;
    o << mol.getNConstraintPairs() << " constraint pairs" << std::endl;
    o << mol.getNFluctuatingCharges() << " fluctuating charges" << std::endl;
    return o;
  }
  
}//end namespace OpenMD
