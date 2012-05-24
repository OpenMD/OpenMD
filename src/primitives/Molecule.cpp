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
 * @file Molecule.cpp
 * @author    tlin
 * @date  10/28/2004
 * @version 1.0
 */ 

#include <algorithm>
#include <set>

#include "primitives/Molecule.hpp"
#include "utils/MemoryUtils.hpp"
#include "utils/simError.h"

namespace OpenMD {
  Molecule::Molecule(int stampId, int globalIndex, const std::string& molName) 
    : stampId_(stampId), globalIndex_(globalIndex), moleculeName_(molName), 
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
    AtomIterator ai;
    RigidBody* rb;
    RigidBodyIterator rbIter;

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
    
    for (inversion = beginInversion(inversionIter); torsion != NULL; 
         inversion =  nextInversion(inversionIter)) {
      potential += inversion->getPotential();
    }
    
    return potential;
    
  }
  
  void Molecule::addProperty(GenericData* genData) {
    properties_.addProperty(genData);  
  }

  void Molecule::removeProperty(const std::string& propName) {
    properties_.removeProperty(propName);  
  }

  void Molecule::clearProperties() {
    properties_.clearProperties(); 
  }

  std::vector<std::string> Molecule::getPropertyNames() {
    return properties_.getPropertyNames();  
  }
      
  std::vector<GenericData*> Molecule::getProperties() { 
    return properties_.getProperties(); 
  }

  GenericData* Molecule::getPropertyByName(const std::string& propName) {
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
