/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the following paper when you publish your work:
 *
 * [1] Drisko et al., J. Open Source Softw. 9, 7004 (2024).
 *
 * Good starting points for code and simulation methodology are:
 *
 * [2] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [3] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [4] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [7] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 * [9] Drisko & Gezelter, J. Chem. Theory Comput. 20, 4986-4997 (2024).
 */

/**
 * @file Molecule.cpp
 * @author    tlin
 * @date  10/28/2004
 * @version 1.0
 */

#include "primitives/Molecule.hpp"

#include <algorithm>
#include <memory>
#include <set>

#include "types/FixedChargeAdapter.hpp"
#include "utils/MemoryUtils.hpp"
#include "utils/StringUtils.hpp"
#include "utils/simError.h"

namespace OpenMD {
  Molecule::Molecule(int globalIndex, MoleculeStamp* molStamp) :
      globalIndex_(globalIndex), molStamp_ {molStamp},
      constrainTotalCharge_(false) {}

  Molecule::~Molecule() {
    Utils::deletePointers(atoms_);
    Utils::deletePointers(bonds_);
    Utils::deletePointers(bends_);
    Utils::deletePointers(torsions_);
    Utils::deletePointers(inversions_);
    Utils::deletePointers(rigidBodies_);
    Utils::deletePointers(cutoffGroups_);
    Utils::deletePointers(constraintPairs_);
    Utils::deletePointers(constraintElems_);
    Utils::deletePointers(hBondDonors_);

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

  void Molecule::addRigidBody(RigidBody* rb) {
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

    // add any atom that wasn't part of a rigid body to the list of
    // integrableObjects

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
      if (atom->isFluctuatingCharge()) fluctuatingCharges_.push_back(atom);
    }

    // find the electronegative atoms and add them to the
    // hBondAcceptors_ vector:

    for (atom = beginAtom(ai); atom != NULL; atom = nextAtom(ai)) {
      AtomType* at = atom->getAtomType();
      // get the chain of base types for this atom type:
      std::vector<AtomType*> ayb = at->allYourBase();
      // use the last type in the chain of base types for the name:
      std::string bn = UpperCase(ayb[ayb.size() - 1]->getName());

      if (bn.compare("O") == 0 || bn.compare("N") == 0 || bn.compare("F") == 0)
        hBondAcceptors_.push_back(atom);
    }

    // find electronegative atoms that are either bonded to
    // hydrogens or are present in the same rigid bodies:

    for (bond = beginBond(bi); bond != NULL; bond = nextBond(bi)) {
      Atom* atom1   = bond->getAtomA();
      Atom* atom2   = bond->getAtomB();
      AtomType* at1 = atom1->getAtomType();
      AtomType* at2 = atom2->getAtomType();
      // get the chain of base types for this atom type:
      std::vector<AtomType*> ayb1 = at1->allYourBase();
      std::vector<AtomType*> ayb2 = at2->allYourBase();
      // use the last type in the chain of base types for the name:
      std::string bn1 = UpperCase(ayb1[ayb1.size() - 1]->getName());
      std::string bn2 = UpperCase(ayb2[ayb2.size() - 1]->getName());

      if (bn1.compare("H") == 0) {
        if (bn2.compare("O") == 0 || bn2.compare("N") == 0 ||
            bn2.compare("F") == 0) {
          HBondDonor* donor      = new HBondDonor();
          donor->donorAtom       = atom2;
          donor->donatedHydrogen = atom1;
          hBondDonors_.push_back(donor);
        }
      }
      if (bn2.compare("H") == 0) {
        if (bn1.compare("O") == 0 || bn1.compare("N") == 0 ||
            bn1.compare("F") == 0) {
          HBondDonor* donor      = new HBondDonor();
          donor->donorAtom       = atom1;
          donor->donatedHydrogen = atom2;
          hBondDonors_.push_back(donor);
        }
      }
    }

    for (rb = beginRigidBody(rbIter); rb != NULL; rb = nextRigidBody(rbIter)) {
      for (atom1 = rb->beginAtom(ai); atom1 != NULL; atom1 = rb->nextAtom(ai)) {
        AtomType* at1 = atom1->getAtomType();
        // get the chain of base types for this atom type:
        std::vector<AtomType*> ayb1 = at1->allYourBase();
        // use the last type in the chain of base types for the name:
        std::string bn1 = UpperCase(ayb1[ayb1.size() - 1]->getName());

        if (bn1.compare("O") == 0 || bn1.compare("N") == 0 ||
            bn1.compare("F") == 0) {
          for (atom2 = rb->beginAtom(aj); atom2 != NULL;
               atom2 = rb->nextAtom(aj)) {
            AtomType* at2 = atom2->getAtomType();
            // get the chain of base types for this atom type:
            std::vector<AtomType*> ayb2 = at2->allYourBase();
            // use the last type in the chain of base types for the name:
            std::string bn2 = UpperCase(ayb2[ayb2.size() - 1]->getName());
            if (bn2.compare("H") == 0) {
              HBondDonor* donor      = new HBondDonor();
              donor->donorAtom       = atom1;
              donor->donatedHydrogen = atom2;
              hBondDonors_.push_back(donor);
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

    for (sd = beginIntegrableObject(i); sd != NULL;
         sd = nextIntegrableObject(i)) {
      mass += sd->getMass();
    }

    return mass;
  }

  RealType Molecule::getFixedCharge() {
    Atom* atom;
    std::vector<Atom*>::iterator i;
    RealType q = 0.0;

    for (atom = beginAtom(i); atom != NULL; atom = nextAtom(i)) {
      AtomType* atomType = atom->getAtomType();

      FixedChargeAdapter fca = FixedChargeAdapter(atomType);
      if (fca.isFixedCharge()) { q += fca.getCharge(); }
    }

    return q;
  }

  Vector3d Molecule::getPrevCom() {
    StuntDouble* sd;
    std::vector<StuntDouble*>::iterator i;
    Vector3d com {};
    RealType totalMass {};
    RealType mass {};

    for (sd = beginIntegrableObject(i); sd != NULL;
         sd = nextIntegrableObject(i)) {
      mass = sd->getMass();
      totalMass += mass;
      com += sd->getPrevPos() * mass;
    }

    com /= totalMass;

    return com;
  }

  Vector3d Molecule::getCom() {
    StuntDouble* sd;
    std::vector<StuntDouble*>::iterator i;
    Vector3d com {};
    RealType totalMass {};
    RealType mass {};

    for (sd = beginIntegrableObject(i); sd != NULL;
         sd = nextIntegrableObject(i)) {
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
    Vector3d com {};
    RealType totalMass {};
    RealType mass {};

    for (sd = beginIntegrableObject(i); sd != NULL;
         sd = nextIntegrableObject(i)) {
      mass = sd->getMass();
      totalMass += mass;
      com += sd->getPos(snapshotNo) * mass;
    }

    com /= totalMass;

    return com;
  }

  void Molecule::setCom(const Vector3d& newCom) {
    Vector3d delta = newCom - getCom();
    moveCom(delta);
  }

  void Molecule::moveCom(const Vector3d& delta) {
    StuntDouble* sd;
    std::vector<StuntDouble*>::iterator i;

    for (sd = beginIntegrableObject(i); sd != NULL;
         sd = nextIntegrableObject(i)) {
      sd->setPos(sd->getPos() + delta);
    }
  }

  Vector3d Molecule::getComVel() {
    StuntDouble* sd;
    std::vector<StuntDouble*>::iterator i;
    Vector3d velCom;
    RealType totalMass = 0;
    RealType mass;

    for (sd = beginIntegrableObject(i); sd != NULL;
         sd = nextIntegrableObject(i)) {
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
    Molecule::BondIterator bondIter;
    ;
    Molecule::BendIterator bendIter;
    Molecule::TorsionIterator torsionIter;
    Molecule::InversionIterator inversionIter;

    RealType potential = 0.0;

    for (bond = beginBond(bondIter); bond != NULL; bond = nextBond(bondIter)) {
      potential += bond->getPotential();
    }

    for (bend = beginBend(bendIter); bend != NULL; bend = nextBend(bendIter)) {
      potential += bend->getPotential();
    }

    for (torsion = beginTorsion(torsionIter); torsion != NULL;
         torsion = nextTorsion(torsionIter)) {
      potential += torsion->getPotential();
    }

    for (inversion = beginInversion(inversionIter); inversion != NULL;
         inversion = nextInversion(inversionIter)) {
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

  std::vector<std::shared_ptr<GenericData>> Molecule::getProperties() {
    return properties_.getProperties();
  }

  std::shared_ptr<GenericData> Molecule::getPropertyByName(
      const std::string& propName) {
    return properties_.getPropertyByName(propName);
  }

  std::ostream& operator<<(std::ostream& o, Molecule& mol) {
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

}  // namespace OpenMD
