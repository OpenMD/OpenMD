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
 * @file Molecule.hpp
 * @author    tlin
 * @date  10/25/2004
 * @version 1.0
 */

#ifndef PRIMITIVES_MOLECULE_HPP
#define PRIMITIVES_MOLECULE_HPP
#include <iostream>
#include <vector>

#include "math/Vector3.hpp"

namespace OpenMD {

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
    int getGlobalIndex() { return globalIndex_; }

    /**
     * Sets the global index of this molecule.
     * @param new global index to be set
     */
    void setGlobalIndex(int index) { return globalIndex_; }

    /**
     * Returns the local index of this molecule
     * @return the local index of this molecule
     */
    int getLocalIndex() { return localIndex_; }

    /** add an atom into this molecule */
    void addAtom(Atom* atom);

    /** add a bond into this molecule */
    void addBond(Bond* bond);

    /** add a bend into this molecule */
    void addBend(Bend* bend);

    /** add a torsion into this molecule*/
    void addTorsion(Torsion* torsion);

    /** add a rigidbody into this molecule */
    void addRigidBody(RigidBody* rb);

    /** add a cutoff group into this molecule */
    void addCutoffGroup(CutoffGroup* cp)

        /** */
        void complete();

    /** Returns the total number of atoms in this molecule */
    unsigned int getNAtoms() { return atoms_.size(); }

    /** Returns the total number of bonds in this molecule */
    unsigned int getNBonds() { return bonds_.size(); }

    /** Returns the total number of bends in this molecule */
    unsigned int getNBends() { return bends_.size(); }

    /** Returns the total number of torsions in this molecule */
    unsigned int getNTorsions() { return torsions_.size(); }

    /** Returns the total number of rigid bodies in this molecule */
    unsigned int getNRigidBodies() { return rigidBodies_.size(); }

    /** Returns the total number of integrable objects in this molecule */
    unsigned int getNIntegrableObjects() { return integrableObjects_.size(); }

    /** Returns the total number of cutoff groups in this molecule */
    unsigned int getNCutoffGroups() { return cutoffGroups_.size(); }

    /**
     * Returns the first atom in this molecule and initialize the iterator.
     * @return the first atom, return NULL if there is not cut off group in this
     * molecule
     * @param i iteraotr
     */
    Atom* beginAtom(std::vector<Atom*>::iterator& i);

    Atom* nextAtom(std::vector<Atom*>::iterator& i);

    /**
     * Returns the first bond in this molecule and initialize the iterator.
     * @return the first bond, return NULL if there is not cut off group in this
     * molecule
     * @param i iteraotr
     */
    Bond* beginBond(std::vector<Bond*>::iterator& i);

    Bond* nextBond(std::vector<Bond*>::iterator& i);

    /**
     * Returns the first bend in this molecule and initialize the iterator.
     * @return the first bend, return NULL if there is not cut off group in this
     * molecule
     * @param i iteraotr
     */
    Bend* beginBend(std::vector<Bend*>::iterator& i);

    Bend* nextBend(std::vector<Bend*>::iterator& i);

    /**
     * Returns the first torsion in this molecule and initialize the iterator.
     * @return the first torsion, return NULL if there is not cut off group in
     * this molecule
     * @param i iteraotr
     */
    Torsion* beginTorsion(std::vector<Torsion*>::iterator& i);
    Torsion* nextTorsion(std::vector<Torsion*>::iterator& i);

    /**
     * Returns the first rigid body in this molecule and initialize the
     * iterator.
     * @return the first rigid body, return NULL if there is not cut off group
     * in this molecule
     * @param i iteraotr
     */
    RigidBody* beginRigidBody(std::vector<RigidBody*>::iterator& i);

    RigidBody* nextRigidBody(std::vector<RigidBody*>::iterator& i);

    /**
     * Returns the first integrable object in this molecule and initialize the
     * iterator.
     * @return the first integrable object, return NULL if there is not cut off
     * group in this molecule
     * @param i iteraotr
     */
    StuntDouble* beginIntegrableObject(std::vector<StuntDouble*>::iterator& i);

    StuntDouble* nextIntegrableObject(std::vector<StuntDouble*>::iterator& i);

    /**
     * Returns the first cutoff group in this molecule and initialize the
     * iterator.
     * @return the first cutoff group, return NULL if there is not cut off group
     * in this molecule
     * @param i iteraotr
     */
    CutoffGroup* beginCutoffGroup(std::vector<CutoffGroup*>::iterator& i);

    /**
     * Returns next cutoff group based on the iterator
     * @return next cutoff group
     * @param i
     */
    CutoffGroup* nextCutoffGroup(std::vector<CutoffGroup*>::iterator& i);

    // void setStampID( int info ) {stampID = info;}

    void calcForces(void);

    void atoms2rigidBodies(void);

    /** return the total potential energy of short range interaction of this
     * molecule */
    double getPotential();

    /** return the center of mass of this molecule */
    Vector3d getCom();

    /** Moves the center of this molecule */
    void moveCom(const Vetor3d& delta);

    /** Returns the velocity of center of mass of this molecule */
    Vector3d getComVel();

    /** Returns the total mass of this molecule */
    double getTotalMass();

    friend std::ostream& operator<<(std::ostream& o, const Molecule& mol);

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

}  // namespace OpenMD
#endif  //
