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
#include "types/AtomStamp.hpp"

namespace OpenMD {
  class RigidBody : public StuntDouble {
  public:
    using AtomIterator = std::vector<Atom*>::iterator;

    RigidBody();

    virtual std::string getType() { return name_; }

    /** Sets the name of this stuntRealType*/
    virtual void setType(const std::string& name) { name_ = name; }

    /**
     * Sets  the previous rotation matrix of this stuntdouble
     * @param a  new rotation matrix
     */
    virtual void setPrevA(const RotMat3x3d& a);

    /**
     * Sets  the current rotation matrix of this stuntdouble
     * @param a  new rotation matrix
     * @note setA will not change the position and rotation matrix of
     * Directional atoms belong to this rigidbody. If you want to do that, use
     * #updateAtoms
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

    /**
     * Returns the gradient of this stuntdouble
     * @return the gradient of this stuntdouble
     */
    virtual std::vector<RealType> getGrad();

    virtual void accept(BaseVisitor* v);

    void addAtom(Atom* at, AtomStamp* ats);

    /** calculates the reference coordinates */
    void calcRefCoords();

    /** Converts Atomic forces and torques to total forces and torques */
    void calcForcesAndTorques();

    /**
        Converts Atomic forces and torques to total forces and torques
        and computes the rigid body contribution to the virial.
        Returns the rigid body contribution to the virial as a 3x3
        matrix.
        */
    Mat3x3d calcForcesAndTorquesAndVirial();

    /** update the positions of atoms belong to this rigidbody */
    void updateAtoms();

    void updateAtoms(int frame);

    void updateAtomVel();

    void updateAtomVel(int frame);

    Atom* beginAtom(std::vector<Atom*>::iterator& i) {
      i = atoms_.begin();
      return i != atoms_.end() ? *i : NULL;
    }

    Atom* nextAtom(std::vector<Atom*>::iterator& i) {
      ++i;
      return i != atoms_.end() ? *i : NULL;
    }

    std::vector<Atom*>::iterator getBeginAtomIter() { return atoms_.begin(); }

    std::vector<Atom*>::iterator getEndAtomIter() { return atoms_.end(); }

    /**
     * Returns the atoms of this rigid body
     * @return the atoms of this rigid body in a vector
     * @deprecated
     */
    std::vector<Atom*> getAtoms() { return atoms_; }

    /**
     * Returns the number of atoms in this rigid body
     * @return the number of atoms in this rigid body
     */
    size_t getNumAtoms() { return atoms_.size(); }

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

    Vector3d getRefCOM() { return refCOM_; }

  private:
    std::string name_;
    Mat3x3d inertiaTensor_;
    RotMat3x3d sU_; /**< body fixed standard unit vector */

    std::vector<Atom*> atoms_;
    Vector3d refCOM_; /**< center of mass relative to original coordinates */
    std::vector<Vector3d> refCoords_;
    std::vector<RotMat3x3d> refOrients_;
  };
}  // namespace OpenMD

#endif  // PRIMITIVES_RIGIDBODY_HPP
