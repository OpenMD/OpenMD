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

#ifndef PRIMITIVES_MOLECULE_HPP
#define PRIMITIVES_MOLECULE_HPP
#include <vector>
#include <iostream>
#include "math/Vector3.hpp"

namespace oopse{

/**
 * @class Molecule Molecule.hpp "primitives/Molecule.hpp"
 * @brief
 */
class Molecule {
    public:

        Molecule();
        virtual ~Molecule();

        /**
         * Returns the global index of this molecule.
         * @return  the global index of this molecule 
         */
        int getGlobalIndex() {
            return globalIndex_;
        }

        /**
         * Sets the global index of this molecule.
         * @param new global index to be set
         */
        void setGlobalIndex(int index) {
            return globalIndex_;
        }
        
        /** 
         * Returns the local index of this molecule 
         * @return the local index of this molecule
         */
        int getLocalIndex() {
            return localIndex_;
        }

        /** add an atom into this molecule */
        void addAtom(Atom* atom);

        /** add a bond into this molecule */
        void addBond(Bond* bond);

        /** add a bend into this molecule */
        void addBend(Bend* bend);

        /** add a torsion into this molecule*/
        void addTorsion(Torsion* torsion);

        /** add a rigidbody into this molecule */
        void addRigidBody(RigidBody *rb);

        /** add a cutoff group into this molecule */
        void addCutoffGroup(CutoffGroup* cp)         

        /** */
        void complete();

        /** Returns the total number of atoms in this molecule */
        unsigned int getNAtoms() {
            return atoms_.size();
        }

        /** Returns the total number of bonds in this molecule */        
        unsigned int getNBonds(){
            return bonds_.size();
        }

        /** Returns the total number of bends in this molecule */        
        unsigned int getNBends() {
            return bends_.size();
        }

        /** Returns the total number of torsions in this molecule */        
        unsigned int getNTorsions() {
            return torsions_.size();
        }

        /** Returns the total number of rigid bodies in this molecule */        
        unsigned int getNRigidBodies() {
            return rigidBodies_.size();
        }

        /** Returns the total number of integrable objects in this molecule */
        unsigned int getNIntegrableObjects() {
            return integrableObjects_.size();
        }

        /** Returns the total number of cutoff groups in this molecule */
        unsigned int getNCutoffGroups() {
            return cutoffGroups_.size();
        }

        /**
         * Returns the first atom in this molecule and initialize the iterator.
         * @return the first atom, return NULL if there is not cut off group in this molecule
         * @param i iteraotr
         */        
        Atom* beginAtom(std::vector<Atom*>::iterator& i);

        Atom* nextAtom(std::vector<Atom*>::iterator& i);

        /**
         * Returns the first bond in this molecule and initialize the iterator.
         * @return the first bond, return NULL if there is not cut off group in this molecule
         * @param i iteraotr
         */
        Bond* beginBond(std::vector<Bond*>::iterator& i);

        Bond* nextBond(std::vector<Bond*>::iterator& i);

        /**
         * Returns the first bend in this molecule and initialize the iterator.
         * @return the first bend, return NULL if there is not cut off group in this molecule
         * @param i iteraotr
         */
        Bend* beginBend(std::vector<Bend*>::iterator& i);

        Bend* nextBend(std::vector<Bend*>::iterator& i);

        /**
         * Returns the first torsion in this molecule and initialize the iterator.
         * @return the first torsion, return NULL if there is not cut off group in this molecule
         * @param i iteraotr
         */
        Torsion* beginTorsion(std::vector<Torsion*>::iterator& i);
        Torsion* nextTorsion(std::vector<Torsion*>::iterator& i);

        /**
         * Returns the first rigid body in this molecule and initialize the iterator.
         * @return the first rigid body, return NULL if there is not cut off group in this molecule
         * @param i iteraotr
         */
        RigidBody* beginRigidBody(std::vector<RigidBody*>::iterator& i);

        RigidBody* nextRigidBody(std::vector<RigidBody*>::iterator& i);

        /**
         * Returns the first integrable object in this molecule and initialize the iterator.
         * @return the first integrable object, return NULL if there is not cut off group in this molecule
         * @param i iteraotr
         */
        StuntDouble* beginIntegrableObject(std::vector<StuntDouble*>::iterator& i);
        
        StuntDouble* nextIntegrableObject(std::vector<StuntDouble*>::iterator& i);

        /**
         * Returns the first cutoff group in this molecule and initialize the iterator.
         * @return the first cutoff group, return NULL if there is not cut off group in this molecule
         * @param i iteraotr
         */
        CutoffGroup* beginCutoffGroup(std::vector<CutoffGroup*>::iterator& i);

        /**
         * Returns next cutoff group based on the iterator 
         * @return next cutoff group
         * @param i 
         */        
        CutoffGroup* nextCutoffGroup(std::vector<CutoffGroup*>::iterator& i);

        //void setStampID( int info ) {stampID = info;} 

        void calcForces( void );

        void atoms2rigidBodies( void );

        /** return the total potential energy of short range interaction of this molecule */
        double getPotential();


        /** return the center of mass of this molecule */
        Vector3d getCom();

        /** Moves the center of this molecule */
        void moveCom(const Vetor3d& delta);

        /** Returns the velocity of center of mass of this molecule */
        Vector3d getComVel();

        /** Returns the total mass of this molecule */
        double getTotalMass();

        friend std::ostream& operator <<(std::ostream& o, const Molecule& mol);
        
    private:
        int localIndex_;
        int globalIndex_;

        std::vector<Atom*> atoms_;
        std::vector<Bond*> bonds_;
        std::vector<Bend*> bends_;
        std::vector<Torsion*> torsions_;
        std::vector<RigidBody*> rigidBodies_;
        std::vector<StuntDouble*> integrableObjects_;
        std::vector<CutoffGroup*> cutoffGroups_;
};

} //namespace oopse
#endif //
