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
 * @file Molecule.hpp
 * @author    tlin
 * @date  10/25/2004
 * @version 1.0
 */

#ifndef PRIMITIVES_MOLECULE_HPP
#define PRIMITIVES_MOLECULE_HPP

#include <iostream>
#include <memory>
#include <vector>

#include "constraints/ConstraintPair.hpp"
#include "math/Vector3.hpp"
#include "primitives/Atom.hpp"
#include "primitives/Bend.hpp"
#include "primitives/Bond.hpp"
#include "primitives/CutoffGroup.hpp"
#include "primitives/Inversion.hpp"
#include "primitives/RigidBody.hpp"
#include "primitives/Torsion.hpp"
#include "utils/PropertyMap.hpp"

namespace OpenMD {

  class Constraint;

  /**
   * @class Molecule Molecule.hpp "primitives/Molecule.hpp"
   * @brief
   */
  class Molecule {
  public:
    struct HBondDonor {
      Atom* donorAtom;
      Atom* donatedHydrogen;
    };

    using AtomIterator              = std::vector<Atom*>::iterator;
    using BondIterator              = std::vector<Bond*>::iterator;
    using BendIterator              = std::vector<Bend*>::iterator;
    using TorsionIterator           = std::vector<Torsion*>::iterator;
    using InversionIterator         = std::vector<Inversion*>::iterator;
    using RigidBodyIterator         = std::vector<RigidBody*>::iterator;
    using CutoffGroupIterator       = std::vector<CutoffGroup*>::iterator;
    using IntegrableObjectIterator  = std::vector<StuntDouble*>::iterator;
    using ConstraintPairIterator    = std::vector<ConstraintPair*>::iterator;
    using ConstraintElemIterator    = std::vector<ConstraintElem*>::iterator;
    using FluctuatingChargeIterator = std::vector<Atom*>::iterator;
    using HBondDonorIterator        = std::vector<HBondDonor*>::iterator;
    using HBondAcceptorIterator     = std::vector<Atom*>::iterator;

    Molecule(int globalIndex, MoleculeStamp* molStamp);
    virtual ~Molecule();

    /**
     * Returns the global index of this molecule.
     * @return  the global index of this molecule
     */
    int getGlobalIndex() { return globalIndex_; }

    /**
     * Returns the stamp id of this molecule
     * @note Ideally, every molecule should keep a pointer of its
     * molecule stamp instead of its stamp id. However, the pointer
     * will become invalid, if the molecule migrate to other
     * processor.
     */
    int getStampId() { return molStamp_->getIdent(); }
    int getRegion() { return molStamp_->getRegion(); }

    MoleculeStamp* getMolStamp() const { return molStamp_; }

    /** Returns the name of the molecule */
    std::string getType() { return molStamp_->getName(); }

    /**
     * Sets the global index of this molecule.
     * @param index new global index to be set
     */
    void setGlobalIndex(int index) { globalIndex_ = index; }

    void setConstrainTotalCharge(bool ctc) { constrainTotalCharge_ = ctc; }

    bool constrainTotalCharge() { return constrainTotalCharge_; }

    /** add an atom into this molecule */
    void addAtom(Atom* atom);

    /** add a bond into this molecule */
    void addBond(Bond* bond);

    /** add a bend into this molecule */
    void addBend(Bend* bend);

    /** add a torsion into this molecule*/
    void addTorsion(Torsion* torsion);

    /** add an improper torsion into this molecule*/
    void addInversion(Inversion* inversion);

    /** add a rigidbody into this molecule */
    void addRigidBody(RigidBody* rb);

    /** add a cutoff group into this molecule */
    void addCutoffGroup(CutoffGroup* cp);

    void addConstraintPair(ConstraintPair* consPair);

    void addConstraintElem(ConstraintElem* consElem);

    /** */
    void complete();

    /** Returns the total number of atoms in this molecule */
    size_t getNAtoms() { return atoms_.size(); }

    /** Returns the total number of bonds in this molecule */
    size_t getNBonds() { return bonds_.size(); }

    /** Returns the total number of bends in this molecule */
    size_t getNBends() { return bends_.size(); }

    /** Returns the total number of torsions in this molecule */
    size_t getNTorsions() { return torsions_.size(); }

    /** Returns the total number of improper torsions in this molecule */
    size_t getNInversions() { return inversions_.size(); }

    /** Returns the total number of rigid bodies in this molecule */
    size_t getNRigidBodies() { return rigidBodies_.size(); }

    /** Returns the total number of integrable objects in this molecule */
    size_t getNIntegrableObjects() { return integrableObjects_.size(); }

    /** Returns the total number of cutoff groups in this molecule */
    size_t getNCutoffGroups() { return cutoffGroups_.size(); }

    /** Returns the total number of constraints in this molecule */
    size_t getNConstraintPairs() { return constraintPairs_.size(); }

    /** Returns the total number of fluctuating charges in this molecule */
    size_t getNFluctuatingCharges() { return fluctuatingCharges_.size(); }
    /** Returns the total number of Hydrogen Bond donors in this molecule */
    size_t getNHBondDonors() { return hBondDonors_.size(); }

    /** Returns the total number of Hydrogen Bond acceptors in this molecule */
    size_t getNHBondAcceptors() { return hBondAcceptors_.size(); }

    Atom* getAtomAt(unsigned int i) {
      assert(i < atoms_.size());
      return atoms_[i];
    }

    RigidBody* getRigidBodyAt(unsigned int i) {
      assert(i < rigidBodies_.size());
      return rigidBodies_[i];
    }

    Atom* beginAtom(std::vector<Atom*>::iterator& i) {
      i = atoms_.begin();
      return (i == atoms_.end()) ? NULL : *i;
    }

    Atom* nextAtom(std::vector<Atom*>::iterator& i) {
      ++i;
      return (i == atoms_.end()) ? NULL : *i;
    }

    Bond* beginBond(std::vector<Bond*>::iterator& i) {
      i = bonds_.begin();
      return (i == bonds_.end()) ? NULL : *i;
    }

    Bond* nextBond(std::vector<Bond*>::iterator& i) {
      ++i;
      return (i == bonds_.end()) ? NULL : *i;
    }

    Bend* beginBend(std::vector<Bend*>::iterator& i) {
      i = bends_.begin();
      return (i == bends_.end()) ? NULL : *i;
    }

    Bend* nextBend(std::vector<Bend*>::iterator& i) {
      ++i;
      return (i == bends_.end()) ? NULL : *i;
    }

    Torsion* beginTorsion(std::vector<Torsion*>::iterator& i) {
      i = torsions_.begin();
      return (i == torsions_.end()) ? NULL : *i;
    }

    Torsion* nextTorsion(std::vector<Torsion*>::iterator& i) {
      ++i;
      return (i == torsions_.end()) ? NULL : *i;
    }

    Inversion* beginInversion(std::vector<Inversion*>::iterator& i) {
      i = inversions_.begin();
      return (i == inversions_.end()) ? NULL : *i;
    }

    Inversion* nextInversion(std::vector<Inversion*>::iterator& i) {
      ++i;
      return (i == inversions_.end()) ? NULL : *i;
    }

    RigidBody* beginRigidBody(std::vector<RigidBody*>::iterator& i) {
      i = rigidBodies_.begin();
      return (i == rigidBodies_.end()) ? NULL : *i;
    }

    RigidBody* nextRigidBody(std::vector<RigidBody*>::iterator& i) {
      ++i;
      return (i == rigidBodies_.end()) ? NULL : *i;
    }

    StuntDouble* beginIntegrableObject(std::vector<StuntDouble*>::iterator& i) {
      i = integrableObjects_.begin();
      return (i == integrableObjects_.end()) ? NULL : *i;
    }

    StuntDouble* nextIntegrableObject(std::vector<StuntDouble*>::iterator& i) {
      ++i;
      return (i == integrableObjects_.end()) ? NULL : *i;
    }

    CutoffGroup* beginCutoffGroup(std::vector<CutoffGroup*>::iterator& i) {
      i = cutoffGroups_.begin();
      return (i == cutoffGroups_.end()) ? NULL : *i;
    }

    CutoffGroup* nextCutoffGroup(std::vector<CutoffGroup*>::iterator& i) {
      ++i;
      return (i == cutoffGroups_.end()) ? NULL : *i;
    }

    ConstraintPair* beginConstraintPair(
        std::vector<ConstraintPair*>::iterator& i) {
      i = constraintPairs_.begin();
      return (i == constraintPairs_.end()) ? NULL : *i;
    }

    ConstraintPair* nextConstraintPair(
        std::vector<ConstraintPair*>::iterator& i) {
      ++i;
      return (i == constraintPairs_.end()) ? NULL : *i;
    }

    ConstraintElem* beginConstraintElem(
        std::vector<ConstraintElem*>::iterator& i) {
      i = constraintElems_.begin();
      return (i == constraintElems_.end()) ? NULL : *i;
    }

    ConstraintElem* nextConstraintElem(
        std::vector<ConstraintElem*>::iterator& i) {
      ++i;
      return (i == constraintElems_.end()) ? NULL : *i;
    }

    Atom* beginFluctuatingCharge(std::vector<Atom*>::iterator& i) {
      i = fluctuatingCharges_.begin();
      return (i == fluctuatingCharges_.end()) ? NULL : *i;
    }

    Atom* nextFluctuatingCharge(std::vector<Atom*>::iterator& i) {
      ++i;
      return (i == fluctuatingCharges_.end()) ? NULL : *i;
    }

    HBondDonor* beginHBondDonor(std::vector<HBondDonor*>::iterator& i) {
      i = hBondDonors_.begin();
      return (i == hBondDonors_.end()) ? NULL : *i;
    }

    HBondDonor* nextHBondDonor(std::vector<HBondDonor*>::iterator& i) {
      ++i;
      return (i == hBondDonors_.end()) ? NULL : *i;
    }

    Atom* beginHBondAcceptor(std::vector<Atom*>::iterator& i) {
      i = hBondAcceptors_.begin();
      return (i == hBondAcceptors_.end()) ? NULL : *i;
    }

    Atom* nextHBondAcceptor(std::vector<Atom*>::iterator& i) {
      ++i;
      return (i == hBondAcceptors_.end()) ? NULL : *i;
    }

    /**
     * Returns the total potential energy of short range interaction
     * of this molecule
     */
    RealType getPotential();

    /** get total mass of this molecule */
    RealType getMass();

    /** get total fixed charge of this molecule */
    RealType getFixedCharge();

    /**
     * Returns the center of mass position of this molecule in
     * the previous snapshot
     *
     * @return the center of mass position of this molecule.
     */
    Vector3d getPrevCom();

    /**
     * Returns the current center of mass position of this molecule.
     *
     * @return the center of mass position of this molecule.
     */
    Vector3d getCom();

    /**
     * Returns the center of mass position of this molecule in
     * specified snapshot
     *
     * @return the center of mass position of this molecule
     * @param snapshotNo
     */
    Vector3d getCom(int snapshotNo);

    /** Sets the center of this molecule */
    void setCom(const Vector3d& newCom);

    /** Moves the center of this molecule */
    void moveCom(const Vector3d& delta);

    /** Returns the velocity of center of mass of this molecule */
    Vector3d getComVel();

    friend std::ostream& operator<<(std::ostream& o, Molecule& mol);

    // below functions are just forward functions
    /**
     * Adds property into property map
     * @param genData GenericData to be added into PropertyMap
     */
    void addProperty(std::shared_ptr<GenericData> genData);

    /**
     * Removes property from PropertyMap by name
     * @param propName the name of property to be removed
     */
    void removeProperty(const std::string& propName);

    /**
     * Returns all names of properties
     * @return all names of properties
     */
    std::vector<std::string> getPropertyNames();

    /**
     * Returns all of the properties in PropertyMap
     * @return all of the properties in PropertyMap
     */
    std::vector<std::shared_ptr<GenericData>> getProperties();

    /**
     * Returns property
     * @param propName name of property
     * @return a pointer point to property with propName. If no property named
     * propName exists, return NULL
     */
    std::shared_ptr<GenericData> getPropertyByName(const std::string& propName);

  private:
    int globalIndex_;

    std::vector<Atom*> atoms_;
    std::vector<Bond*> bonds_;
    std::vector<Bend*> bends_;
    std::vector<Torsion*> torsions_;
    std::vector<Inversion*> inversions_;
    std::vector<RigidBody*> rigidBodies_;
    std::vector<StuntDouble*> integrableObjects_;
    std::vector<CutoffGroup*> cutoffGroups_;
    std::vector<ConstraintPair*> constraintPairs_;
    std::vector<ConstraintElem*> constraintElems_;
    std::vector<Atom*> fluctuatingCharges_;
    std::vector<HBondDonor*> hBondDonors_;
    std::vector<Atom*> hBondAcceptors_;

    PropertyMap properties_;
    bool constrainTotalCharge_;
    MoleculeStamp* molStamp_;
  };
}  // namespace OpenMD

#endif  //
