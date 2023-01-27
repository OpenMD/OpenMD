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
 * @file RigidBody.hpp
 * @author    tlin
 * @date  10/23/2004
 * @version 1.0
 */

#ifndef PRIMITIVES_RIGIDBODY_HPP
#define PRIMITIVES_RIGIDBODY_HPP

#include <vector>

#include "primitives/DirectionalAtom.hpp"
#include "primitives/StuntDouble.hpp"
namespace OpenMD {
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

    /** Sets the internal unit frame of this stuntdouble by three euler angles
     */
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
    std::vector<Atom*> getAtoms() { return atoms_; }

    /**
     * Returns the number of atoms in this rigid body
     * @return the number of atoms in this rigid body
     */
    int getNumAtoms() { return atoms_.size(); }

    /**
     * Return the position of atom which belongs to this rigid body.
     * @return true if index is valid otherwise return false
     * @param pos the position of atom which will be set on return if index is
     * valid
     * @param index the index of the atom in rigid body's private data member
     * atoms_
     */
    bool getAtomPos(Vector3d& pos, unsigned int index);

    /**
     * Return the position of atom which belongs to this rigid body.
     * @return true if atom belongs to this rigid body,otherwise return false
     * @param pos position of atom which will be set on return if atom belongs
     * to this rigid body
     * @param atom the pointer to an atom
     */
    bool getAtomPos(Vector3d& pos, Atom* atom);

    /**
     * Return the velocity of atom which belongs to this rigid body.
     * @return true if index is valid otherwise return false
     * @param vel the velocity of atom which will be set on return if index is
     * valid
     * @param index the index of the atom in rigid body's private data member
     * atoms_
     */
    bool getAtomVel(Vector3d& vel, unsigned int index);

    /**
     * Return the velocity of atom which belongs to this rigid body.
     * @return true if atom belongs to this rigid body,otherwise return false
     * @param vel velocity of atom which will be set on return if atom belongs
     * to this rigid body
     * @param atom the pointer to an atom
     */
    bool getAtomVel(Vector3d& vel, Atom*);

    /**
     * Return the reference coordinate of atom which belongs to this rigid body.
     * @return true if index is valid otherwise return false
     * @param coor the reference coordinate of atom which will be set on return
     * if index is valid
     * @param index the index of the atom in rigid body's private data member
     * atoms_
     */
    bool getAtomRefCoor(Vector3d& coor, unsigned int index);

    /**
     * Return the velocity of atom which belongs to this rigid body.
     * @return true if atom belongs to this rigid body,otherwise return false
     * @param coor velocity of atom which will be set on return if atom belongs
     * to this rigid body
     * @param atom the pointer to an atom
     */
    bool getAtomRefCoor(Vector3d& coor, Atom* atom);

  private:
    Mat3x3d inertiaTensor_;
    RotMat3x3d sU_; /**< body fixed standard unit vector */

    std::vector<Atom*> atoms_;
    std::vector<Vector3d> refCoords_;
  };

}  // namespace OpenMD

#endif  // PRIMITIVES_RIGIDBODY_HPP
