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
 * @file Molecule.hpp
 * @author    tlin
 * @date  10/25/2004
 * @version 1.0
 */ 

#include "primitives/Molecule.hpp"

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

}//end namespace oopse
