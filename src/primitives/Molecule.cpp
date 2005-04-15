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

namespace oopse {
  Molecule::Molecule(int stampId, int globalIndex, const std::string& molName) 
    : stampId_(stampId), globalIndex_(globalIndex), moleculeName_(molName) {

    }

  Molecule::~Molecule() {

    MemoryUtils::deletePointers(atoms_);
    MemoryUtils::deletePointers(bonds_);
    MemoryUtils::deletePointers(bends_);
    MemoryUtils::deletePointers(torsions_);
    MemoryUtils::deletePointers(rigidBodies_);
    MemoryUtils::deletePointers(cutoffGroups_);
    MemoryUtils::deletePointers(constraintPairs_);
    MemoryUtils::deletePointers(constraintElems_);
    //integrableObjects_ don't own the objects
    integrableObjects_.clear();
    
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
    if (std::find(torsions_.begin(), torsions_.end(), torsion) == torsions_.end()) {
      torsions_.push_back(torsion);
    }
  }

  void Molecule::addRigidBody(RigidBody *rb) {
    if (std::find(rigidBodies_.begin(), rigidBodies_.end(), rb) == rigidBodies_.end()) {
      rigidBodies_.push_back(rb);
    }
  }

  void Molecule::addCutoffGroup(CutoffGroup* cp) {
    if (std::find(cutoffGroups_.begin(), cutoffGroups_.end(), cp) == cutoffGroups_.end()) {
      cutoffGroups_.push_back(cp);
    }

  }

  void Molecule::addConstraintPair(ConstraintPair* cp) {
    if (std::find(constraintPairs_.begin(), constraintPairs_.end(), cp) == constraintPairs_.end()) {
      constraintPairs_.push_back(cp);
    }

  }

  void Molecule::addConstraintElem(ConstraintElem* cp) {
    if (std::find(constraintElems_.begin(), constraintElems_.end(), cp) == constraintElems_.end()) {
      constraintElems_.push_back(cp);
    }

  }

  void Molecule::complete() {
    
    std::set<Atom*> rigidAtoms;
    RigidBody* rb;
    std::vector<RigidBody*>::iterator rbIter;

    
    for (rb = beginRigidBody(rbIter); rb != NULL; rb = nextRigidBody(rbIter)) {
      rigidAtoms.insert(rb->getBeginAtomIter(), rb->getEndAtomIter());
    }

    Atom* atom;
    AtomIterator ai;
    for (atom = beginAtom(ai); atom != NULL; atom = nextAtom(ai)) {
   
      if (rigidAtoms.find(*ai) == rigidAtoms.end()) {
	//if an atom does not belong to a rigid body, it is an integrable object
	integrableObjects_.push_back(*ai);
      }
    }

    //find all free atoms (which do not belong to rigid bodies)  
    //performs the "difference" operation from set theory,  the output range contains a copy of every
    //element that is contained in [allAtoms.begin(), allAtoms.end()) and not contained in 
    //[rigidAtoms.begin(), rigidAtoms.end()).
    //std::set_difference(allAtoms.begin(), allAtoms.end(), rigidAtoms.begin(), rigidAtoms.end(),
    //                        std::back_inserter(integrableObjects_));

    //if (integrableObjects_.size() != allAtoms.size() - rigidAtoms.size()) {
    //    //Some atoms in rigidAtoms are not in allAtoms, something must be wrong
    //    sprintf(painCave.errMsg, "Atoms in rigidbody are not in the atom list of the same molecule");
    //
    //    painCave.isFatal = 1;
    //    simError();        
    //}
    for (rb = beginRigidBody(rbIter); rb != NULL; rb = nextRigidBody(rbIter)) {
      integrableObjects_.push_back(rb);
    } 
    //integrableObjects_.insert(integrableObjects_.end(), rigidBodies_.begin(), rigidBodies_.end());
  }

  double Molecule::getMass() {
    StuntDouble* sd;
    std::vector<StuntDouble*>::iterator i;
    double mass = 0.0;

    for (sd = beginIntegrableObject(i); sd != NULL; sd = nextIntegrableObject(i)){
      mass += sd->getMass();
    }

    return mass;

  }

  Vector3d Molecule::getCom() {
    StuntDouble* sd;
    std::vector<StuntDouble*>::iterator i;
    Vector3d com;
    double totalMass = 0;
    double mass;
    
    for (sd = beginIntegrableObject(i); sd != NULL; sd = nextIntegrableObject(i)){
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
    
    for (sd = beginIntegrableObject(i); sd != NULL; sd = nextIntegrableObject(i)){
      sd->setPos(sd->getPos() + delta);
    }

  }

  Vector3d Molecule::getComVel() {
    StuntDouble* sd;
    std::vector<StuntDouble*>::iterator i;
    Vector3d velCom;
    double totalMass = 0;
    double mass;
    
    for (sd = beginIntegrableObject(i); sd != NULL; sd = nextIntegrableObject(i)){
      mass = sd->getMass();
      totalMass += mass;
      velCom += sd->getVel() * mass;    
    }

    velCom /= totalMass;

    return velCom;
  }

  double Molecule::getPotential() {

    Bond* bond;
    Bend* bend;
    Torsion* torsion;
    Molecule::BondIterator bondIter;;
    Molecule::BendIterator  bendIter;
    Molecule::TorsionIterator  torsionIter;

    double potential = 0.0;

    for (bond = beginBond(bondIter); bond != NULL; bond = nextBond(bondIter)) {
      potential += bond->getPotential();
    }

    for (bend = beginBend(bendIter); bend != NULL; bend = nextBend(bendIter)) {
      potential += bend->getPotential();
    }

    for (torsion = beginTorsion(torsionIter); torsion != NULL; torsion = nextTorsion(torsionIter)) {
      potential += torsion->getPotential();
    }

    return potential;

  }

  std::ostream& operator <<(std::ostream& o, Molecule& mol) {
    o << std::endl;
    o << "Molecule " << mol.getGlobalIndex() << "has: " << std::endl;
    o << mol.getNAtoms() << " atoms" << std::endl;
    o << mol.getNBonds() << " bonds" << std::endl;
    o << mol.getNBends() << " bends" << std::endl;
    o << mol.getNTorsions() << " torsions" << std::endl;
    o << mol.getNRigidBodies() << " rigid bodies" << std::endl;
    o << mol.getNIntegrableObjects() << "integrable objects" << std::endl;
    o << mol.getNCutoffGroups() << "cutoff groups" << std::endl;
    o << mol.getNConstraintPairs() << "constraint pairs" << std::endl;
    return o;
  }

}//end namespace oopse
