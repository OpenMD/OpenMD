/*
 * Copyright (C) 2000-2004  Object Oriented Parallel Simulation Engine (OOPSE) project
 * 
 * Contact: oopse@oopse.org
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 */

/**
 * @file Molecule.cpp
 * @author    tlin
 * @date  10/28/2004
 * @version 1.0
 */ 

#include "primitives/Molecule.hpp"
#include <algorithm>

namespace oopse {

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

}//end namespace oopse
