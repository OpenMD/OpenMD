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

#include "primitives/Molecule.hpp"
#include <algorithm>

namespace OpenMD {

Molecule::~Molecule() {

    deleteVectorOfPointer(atoms_);
    deleteVectorOfPointer(bonds_);
    deleteVectorOfPointer(bends_);
    deleteVectorOfPointer(torsions_);
    deleteVectorOfPointer(rigidbodies_);
    deleteVectorOfPointer(cutoffGroups_);

    integrableObjects_.clear();
    
}

void Molecule::addAtom(Atom* atom) {
    if (atoms_.find(atom) == atoms_.end()) {
        atoms_.push_back(atom);
    }
}

void Molecule::addBond(Bond* bond) {
    if (bonds_.find(bond) == bonds_.end()) {
        bonds_.push_back(bond);
    }
}

void Molecule::addBend(Bend* bend) {
    if (bends_.find(bend) == bends_.end()) {
        bends_.push_back(bend);
    }
}

void Molecule::addTorsion(Torsion* torsion) {
    if (torsions_.find(torsion) == torsions_.end()) {
        torsions_.push_back(torsion);
    }
}

void Molecule::addRigidBody(RigidBody *rb) {
    if (rigidBodies_.find(bond) == bonds_.end()) {
        rigidBodies_.push_back(rb);
    }
}

void Molecule::addCutoffGroup(CutoffGroup* cp) {
    if (cutoffGroups_.find(bond) == bonds_.end()) {
        cutoffGroups_.push_back(cp);
    }

}

void Molecule::complete() {
    
    std::set<Atom*> allAtoms;
    allAtoms.insert(atoms_.begin(), atoms_.end());

    std::set<Atom*> rigidAtoms;
    RigidBody* rb;
    std::vector<RigidBody*> rbIter;

    
    for (rb = beginRigidBody(rbIter); rb != NULL; rb = nextRigidBody(rbIter)) {
        rigidAtoms.insert(rb->beginAtomIter(), rb->endAtomIter());
    }

    //find all free atoms (which do not belong to rigid bodies)  
    //performs the "difference" operation from set theory,  the output range contains a copy of every
    //element that is contained in [allAtoms.begin(), allAtoms.end()) and not contained in 
    //[rigidAtoms.begin(), rigidAtoms.end()).
    std::set_difference(allAtoms.begin(), allAtoms.end(), rigidAtoms.begin(), rigidAtoms.end(),
                            std::back_inserter(integrableObjects_));

    if (integrableObjects_.size() != allAtoms.size() - rigidAtoms.size()) {
        //Some atoms in rigidAtoms are not in allAtoms, something must be wrong
    }

    integrableObjects_.insert(integrableObjects_.end(), rigidBodies_.begin(), rigidBodies_.end());
}

Atom* Molecule::beginAtom(std::vector<Atom*>::iterator& i) {
    i = atoms_.begin();
    return (i == atoms_.end()) ? NULL : *i;
}

Atom* Molecule::nextAtom(std::vector<Atom*>::iterator& i) {
    ++i;
    return (i == atoms_.end()) ? NULL : *i;    
}

Bond* Molecule::beginBond(std::vector<Bond*>::iterator& i) {
    i = bonds_.begin();
    return (i == bonds_.end()) ? NULL : *i;
}

Bond* Molecule::nextBond(std::vector<Bond*>::iterator& i) {
    ++i;
    return (i == bonds_.end()) ? NULL : *i;    

}


Bend* Molecule::beginBend(std::vector<Bend*>::iterator& i) {
    i = bends_.begin();
    return (i == bends_.end()) ? NULL : *i;
}

Bend* Molecule::nextBend(std::vector<Bend*>::iterator& i) {
    ++i;
    return (i == bends_.end()) ? NULL : *i;    
}

Torsion* Molecule::beginTorsion(std::vector<Torsion*>::iterator& i) {
    i = torsions_.begin();
    return (i == torsions_.end()) ? NULL : *i;
}

Torsion* Molecule::nextTorsion(std::vector<Torsion*>::iterator& i) {
    ++i;
    return (i == torsions_.end()) ? NULL : *i;    
}    

RigidBody* Molecule::beginRigidBody(std::vector<RigidBody*>::iterator& i) {
    i = rigidBodies_.begin();
    return (i == rigidBodies_.end()) ? NULL : *i;
}

RigidBody* Molecule::nextRigidBody(std::vector<RigidBody*>::iterator& i) {
    ++i;
    return (i == rigidBodies_.end()) ? NULL : *i;    
}

StuntDouble* Molecule::beginIntegrableObject(std::vector<StuntDouble*>::iterator& i) {
    i = integrableObjects_.begin();
    return (i == integrableObjects_.end()) ? NULL : *i;
}

StuntDouble* Molecule::nextIntegrableObject(std::vector<StuntDouble*>::iterator& i) {
    ++i;
    return (i == integrableObjects_.end()) ? NULL : *i;    
}    

CutoffGroup* Molecule::beginCutoffGroup(std::vector<CutoffGroup*>::iterator& i) {
    i = cutoffGroups_.begin();
    return (i == cutoffGroups_.end()) ? NULL : *i;
}

CutoffGroup* Molecule::nextCutoffGroup(std::vector<CutoffGroup*>::iterator& i) {            
    ++i;
    return (i == cutoffGroups_.end()) ? NULL : *i;    
} 

void Molecule::calcForces() {
    RigidBody* rb;
    Bond* bond;
    Bend* bend;
    Torsion* torsion;
    std::vector<RigidBody*> rbIter;
    std::vector<Bond*> bondIter;;
    std::vector<Bend*> bendIter;
    std::vector<Torsion*> torsionIter;

    for (rb = beginRigidBody(rbIter); rb != NULL; rb = nextRigidBody(rbIter)) {
        rb->updateAtoms();
    }

    for (bond = beginBond(bondIter); bond != NULL; bond = nextBond(bondIter)) {
        bond->calcForce();
    }

    for (bend = beginBend(bendIter); bend != NULL; bend = nextBend(bendIter)) {
        bend->calcForce();
    }

    for (torsion = beginTorsion(torsionIter); torsion != NULL; torsion = nextTorsion(torsionIter)) {
        torsion->calcForce();
    }
    
}

double Molecule::getPotential() {
    //RigidBody* rb;
    Bond* bond;
    Bend* bend;
    Torsion* torsion;
    //std::vector<RigidBody*> rbIter;
    std::vector<Bond*> bondIter;;
    std::vector<Bend*> bendIter;
    std::vector<Torsion*> torsionIter;

    double potential = 0;
    
    //for (rb = beginRigidBody(rbIter); rb != NULL; rb = nextRigidBody(rbIter)) {
    //    rb->updateAtoms();
    //}

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

void Molecule::moveCom(const Vetor3d& delta) {
    StuntDouble* sd;
    std::vector<StuntDouble*>::iterator i;
    
    for (sd = beginIntegrableObject(i); sd != NULL; sd = nextIntegrableObject(i)){
        s->setPos(sd->getPos() + delta);
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

std::ostream& operator <<(std::ostream& o, const Molecule& mol) {
    o << std::endl;
    o << "Molecule " << mol.getGlobalIndex() << "has: " << std::endl;
    o << mol.getNAtoms() << " atoms" << std::endl;
    o << mol.getNBonds() << " bonds" << std::endl;
    o << mol.getNBends() << " bends" << std::endl;
    o << mol.getNTorsions() << " torsions" << std::endl;
    o << mol.getNRigidBodies() << " rigid bodies" << std::endl;
    o << mol.getNIntegrableObjects() << "integrable objects" << std::endl;
    o << mol.getNCutoffGroups() << "cutoff groups" << std::endl;

    return o;
}

}//end namespace OpenMD
