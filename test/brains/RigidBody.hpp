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
 * @file RigidBody.hpp
 * @author    tlin
 * @date  10/23/2004
 * @version 1.0
 */ 

#ifndef PRIMITIVES_RIGIDBODY_HPP
#define PRIMITIVES_RIGIDBODY_HPP

#include <vector>

#include "primitives/StuntDouble.hpp"
#include "primitives/DirectionalAtom.hpp"
namespace oopse{
    class RigidBody : public StuntDouble {
        public:
            RigidBody();

           /**
             * Sets  the previous rotation matrix of this stuntdouble
             * @param a  new rotation matrix 
             */         
           virtual void setPrevA(const RotMat3x3d& a);
           
           /**
             * Sets  the current rotation matrix of this stuntdouble
             * @param a  new rotation matrix 
             */         
            virtual void setA(const RotMat3x3d& a);
           /**
             * Sets  the rotation matrix of this stuntdouble in specified snapshot
             * @param a rotation matrix to be set 
             * @param snapshotNo 
             * @see #getA
             */         
            virtual void setA(const RotMat3x3d& a, int snapshotNo);

            /**
             * Returns the inertia tensor of this stuntdouble
             * @return the inertia tensor of this stuntdouble
             */ 
            virtual Mat3x3d getI();


            /** Sets the internal unit frame of this stuntdouble by three euler angles */
            void setUnitFrameFromEuler(double phi, double theta, double psi);
            
            /**
             * Returns the gradient of this stuntdouble
             * @return the inertia tensor of this stuntdouble
             * @see #setI
             */ 
            virtual std::vector<double> getGrad();

            virtual void accept(BaseVisitor* v);

            void addAtom(Atom* atom);

            /** calculate the reference coordinates */
            void calcRefCoords();

            /** Convert Atomic forces and torques to total forces and torques */
            void calcForcesAndTorques();

            /** update the positions of atoms belong to this rigidbody */
            void updateAtoms();

            /** 
             * Returns the atoms of this rigid body
             * @return the atoms of this rigid body in a vector
             */           
            std::vector<Atom*> getAtoms() {
                return atoms_;
            }

            /** 
             * Returns the number of atoms in this rigid body
             * @return the number of atoms in this rigid body
             */
            int getNumAtoms() {
                return atoms_.size();
            }

            /**
             * Return the position of atom which belongs to this rigid body.
             * @return true if index is valid otherwise return false
             * @param pos the position of atom which will be set on return if index is valid
             * @param index the index of the atom in rigid body's private data member atoms_
             */
            bool getAtomPos(Vector3d& pos, unsigned int index);

            /**
             * Return the position of atom which belongs to this rigid body.
             * @return true if atom belongs to this rigid body,otherwise return false
             * @param pos position of atom which will be set on return if atom belongs to this rigid body
             * @param atom the pointer to an atom
             */            
            bool getAtomPos(Vector3d& pos, Atom* atom);

            /**
             * Return the velocity of atom which belongs to this rigid body.
             * @return true if index is valid otherwise return false
             * @param vel the velocity of atom which will be set on return if index is valid
             * @param index the index of the atom in rigid body's private data member atoms_
             */
            bool getAtomVel(Vector3d& vel, unsigned int index);

            /**
             * Return the velocity of atom which belongs to this rigid body.
             * @return true if atom belongs to this rigid body,otherwise return false
             * @param vel velocity of atom which will be set on return if atom belongs to this rigid body
             * @param atom the pointer to an atom
             */ 
            bool getAtomVel(Vector3d& vel, Atom*);

            /**
             * Return the reference coordinate of atom which belongs to this rigid body.
             * @return true if index is valid otherwise return false
             * @param coor the reference coordinate of atom which will be set on return if index is valid
             * @param index the index of the atom in rigid body's private data member atoms_
             */
            bool getAtomRefCoor(Vector3d& coor, unsigned int index);

            /**
             * Return the velocity of atom which belongs to this rigid body.
             * @return true if atom belongs to this rigid body,otherwise return false
             * @param coor velocity of atom which will be set on return if atom belongs to this rigid body
             * @param atom the pointer to an atom
             */ 
            bool getAtomRefCoor(Vector3d& coor, Atom* atom);

        private:
            
            Mat3x3d inertiaTensor_;     
            RotMat3x3d sU_;               /**< body fixed standard unit vector */
            
            std::vector<Atom*> atoms_;
            std::vector<Vector3d> refCoords_;
    };

}//namepace oopse

#endif //PRIMITIVES_RIGIDBODY_HPP

