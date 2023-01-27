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

#include "selection/IndexFinder.hpp"

#include "primitives/Molecule.hpp"

namespace OpenMD {

  IndexFinder::IndexFinder(SimInfo* info) : info_(info) {
    nObjects_.push_back(info_->getNGlobalAtoms() +
                        info_->getNGlobalRigidBodies());
    nObjects_.push_back(info_->getNGlobalBonds());
    nObjects_.push_back(info_->getNGlobalBends());
    nObjects_.push_back(info_->getNGlobalTorsions());
    nObjects_.push_back(info_->getNGlobalInversions());
    nObjects_.push_back(info_->getNGlobalMolecules());

    selectionSets_.resize(info_->getNGlobalMolecules());
    init();
  }

  void IndexFinder::init() {
    SimInfo::MoleculeIterator mi;
    Molecule::AtomIterator ai;
    Molecule::RigidBodyIterator rbIter;
    Molecule::BondIterator bondIter;
    Molecule::BendIterator bendIter;
    Molecule::TorsionIterator torsionIter;
    Molecule::InversionIterator inversionIter;

    Molecule* mol;
    Atom* atom;
    RigidBody* rb;
    Bond* bond;
    Bend* bend;
    Torsion* torsion;
    Inversion* inversion;

    for (mol = info_->beginMolecule(mi); mol != NULL;
         mol = info_->nextMolecule(mi)) {
      SelectionSet ss(nObjects_);

      ss.bitsets_[MOLECULE].setBitOn(mol->getGlobalIndex());

      for (atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
        ss.bitsets_[STUNTDOUBLE].setBitOn(atom->getGlobalIndex());
      }
      for (rb = mol->beginRigidBody(rbIter); rb != NULL;
           rb = mol->nextRigidBody(rbIter)) {
        ss.bitsets_[STUNTDOUBLE].setBitOn(rb->getGlobalIndex());
      }
      for (bond = mol->beginBond(bondIter); bond != NULL;
           bond = mol->nextBond(bondIter)) {
        ss.bitsets_[BOND].setBitOn(bond->getGlobalIndex());
      }
      for (bend = mol->beginBend(bendIter); bend != NULL;
           bend = mol->nextBend(bendIter)) {
        ss.bitsets_[BEND].setBitOn(bend->getGlobalIndex());
      }
      for (torsion = mol->beginTorsion(torsionIter); torsion != NULL;
           torsion = mol->nextTorsion(torsionIter)) {
        ss.bitsets_[TORSION].setBitOn(torsion->getGlobalIndex());
      }
      for (inversion = mol->beginInversion(inversionIter); inversion != NULL;
           inversion = mol->nextInversion(inversionIter)) {
        ss.bitsets_[INVERSION].setBitOn(inversion->getGlobalIndex());
      }

      selectionSets_[mol->getGlobalIndex()] = ss;
    }
  }

  SelectionSet IndexFinder::find(int molIndex) {
    return selectionSets_[molIndex];
  }

  SelectionSet IndexFinder::find(int begMolIndex, int endMolIndex) {
    SelectionSet ss(nObjects_);

    for (int i = begMolIndex; i < endMolIndex; ++i) {
      ss |= selectionSets_[i];
    }
    return ss;
  }
}  // namespace OpenMD
