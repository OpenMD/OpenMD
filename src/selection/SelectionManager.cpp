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

#include "selection/SelectionManager.hpp"

#ifdef IS_MPI
#include <mpi.h>
#endif

#include "brains/SimInfo.hpp"

namespace OpenMD {
  SelectionManager::SelectionManager(SimInfo* info) : info_(info) {
    nObjects_.push_back(info_->getNGlobalAtoms() +
                        info_->getNGlobalRigidBodies());
    nObjects_.push_back(info_->getNGlobalBonds());
    nObjects_.push_back(info_->getNGlobalBends());
    nObjects_.push_back(info_->getNGlobalTorsions());
    nObjects_.push_back(info_->getNGlobalInversions());
    nObjects_.push_back(info_->getNGlobalMolecules());

    stuntdoubles_.resize(nObjects_[STUNTDOUBLE]);
    bonds_.resize(nObjects_[BOND]);
    bends_.resize(nObjects_[BEND]);
    torsions_.resize(nObjects_[TORSION]);
    inversions_.resize(nObjects_[INVERSION]);
    molecules_.resize(nObjects_[MOLECULE]);

    ss_.resize(nObjects_);

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
      molecules_[mol->getGlobalIndex()] = mol;

      for (atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
        stuntdoubles_[atom->getGlobalIndex()] = atom;
      }
      for (rb = mol->beginRigidBody(rbIter); rb != NULL;
           rb = mol->nextRigidBody(rbIter)) {
        stuntdoubles_[rb->getGlobalIndex()] = rb;
      }
      for (bond = mol->beginBond(bondIter); bond != NULL;
           bond = mol->nextBond(bondIter)) {
        bonds_[bond->getGlobalIndex()] = bond;
      }
      for (bend = mol->beginBend(bendIter); bend != NULL;
           bend = mol->nextBend(bendIter)) {
        bends_[bend->getGlobalIndex()] = bend;
      }
      for (torsion = mol->beginTorsion(torsionIter); torsion != NULL;
           torsion = mol->nextTorsion(torsionIter)) {
        torsions_[torsion->getGlobalIndex()] = torsion;
      }
      for (inversion = mol->beginInversion(inversionIter); inversion != NULL;
           inversion = mol->nextInversion(inversionIter)) {
        inversions_[inversion->getGlobalIndex()] = inversion;
      }
    }
  }

  StuntDouble* SelectionManager::beginSelected(int& i) {
#ifdef IS_MPI
    i = 0;
    while (i < static_cast<int>(ss_.bitsets_[STUNTDOUBLE].size())) {
      if (ss_.bitsets_[STUNTDOUBLE][i]) {
        // check that this processor owns this stuntdouble
        if (stuntdoubles_[i] != NULL) return stuntdoubles_[i];
      }
      ++i;
    }
    return NULL;
#else
    i = ss_.bitsets_[STUNTDOUBLE].firstOnBit();
    return i == -1 ? NULL : stuntdoubles_[i];
#endif
  }

  StuntDouble* SelectionManager::nextSelected(int& i) {
#ifdef IS_MPI
    ++i;
    while (i < static_cast<int>(ss_.bitsets_[STUNTDOUBLE].size())) {
      if (ss_.bitsets_[STUNTDOUBLE][i]) {
        // check that this processor owns this stuntdouble
        if (stuntdoubles_[i] != NULL) return stuntdoubles_[i];
      }
      ++i;
    }
    return NULL;
#else
    i = ss_.bitsets_[STUNTDOUBLE].nextOnBit(i);
    return i == -1 ? NULL : stuntdoubles_[i];
#endif
  }

  StuntDouble* SelectionManager::beginUnselected(int& i) {
#ifdef IS_MPI
    i = 0;
    while (i < static_cast<int>(ss_.bitsets_[STUNTDOUBLE].size())) {
      if (!ss_.bitsets_[STUNTDOUBLE][i]) {
        // check that this processor owns this stuntdouble
        if (stuntdoubles_[i] != NULL) return stuntdoubles_[i];
      }
      ++i;
    }
    return NULL;
#else
    i = ss_.bitsets_[STUNTDOUBLE].firstOffBit();
    return i == -1 ? NULL : stuntdoubles_[i];
#endif
  }

  StuntDouble* SelectionManager::nextUnselected(int& i) {
#ifdef IS_MPI
    ++i;
    while (i < static_cast<int>(ss_.bitsets_[STUNTDOUBLE].size())) {
      if (!ss_.bitsets_[STUNTDOUBLE][i]) {
        // check that this processor owns this stuntdouble
        if (stuntdoubles_[i] != NULL) return stuntdoubles_[i];
      }
      ++i;
    }
    return NULL;
#else
    i = ss_.bitsets_[STUNTDOUBLE].nextOffBit(i);
    return i == -1 ? NULL : stuntdoubles_[i];
#endif
  }

  Bond* SelectionManager::beginSelectedBond(int& i) {
#ifdef IS_MPI
    i = 0;
    while (i < static_cast<int>(ss_.bitsets_[BOND].size())) {
      if (ss_.bitsets_[BOND][i]) {
        // check that this processor owns this bond
        if (bonds_[i] != NULL) return bonds_[i];
      }
      ++i;
    }
    return NULL;
#else
    i = ss_.bitsets_[BOND].firstOnBit();
    return i == -1 ? NULL : bonds_[i];
#endif
  }

  Bond* SelectionManager::nextSelectedBond(int& i) {
#ifdef IS_MPI
    ++i;
    while (i < static_cast<int>(ss_.bitsets_[BOND].size())) {
      if (ss_.bitsets_[BOND][i]) {
        // check that this processor owns this bond
        if (bonds_[i] != NULL) return bonds_[i];
      }
      ++i;
    }
    return NULL;
#else
    i = ss_.bitsets_[BOND].nextOnBit(i);
    return i == -1 ? NULL : bonds_[i];
#endif
  }

  Bond* SelectionManager::beginUnselectedBond(int& i) {
#ifdef IS_MPI
    i = 0;
    while (i < static_cast<int>(ss_.bitsets_[BOND].size())) {
      if (!ss_.bitsets_[BOND][i]) {
        // check that this processor owns this bond
        if (bonds_[i] != NULL) return bonds_[i];
      }
      ++i;
    }
    return NULL;
#else
    i = ss_.bitsets_[BOND].firstOffBit();
    return i == -1 ? NULL : bonds_[i];
#endif
  }

  Bond* SelectionManager::nextUnselectedBond(int& i) {
#ifdef IS_MPI
    ++i;
    while (i < static_cast<int>(ss_.bitsets_[BOND].size())) {
      if (!ss_.bitsets_[BOND][i]) {
        // check that this processor owns this bond
        if (bonds_[i] != NULL) return bonds_[i];
      }
      ++i;
    }
    return NULL;
#else
    i = ss_.bitsets_[BOND].nextOffBit(i);
    return i == -1 ? NULL : bonds_[i];
#endif
  }

  Bend* SelectionManager::beginSelectedBend(int& i) {
#ifdef IS_MPI
    i = 0;
    while (i < static_cast<int>(ss_.bitsets_[BEND].size())) {
      if (ss_.bitsets_[BEND][i]) {
        // check that this processor owns this bend
        if (bends_[i] != NULL) return bends_[i];
      }
      ++i;
    }
    return NULL;
#else
    i = ss_.bitsets_[BEND].firstOnBit();
    return i == -1 ? NULL : bends_[i];
#endif
  }

  Bend* SelectionManager::nextSelectedBend(int& i) {
#ifdef IS_MPI
    ++i;
    while (i < static_cast<int>(ss_.bitsets_[BEND].size())) {
      if (ss_.bitsets_[BEND][i]) {
        // check that this processor owns this bend
        if (bends_[i] != NULL) return bends_[i];
      }
      ++i;
    }
    return NULL;
#else
    i = ss_.bitsets_[BEND].nextOnBit(i);
    return i == -1 ? NULL : bends_[i];
#endif
  }

  Bend* SelectionManager::beginUnselectedBend(int& i) {
#ifdef IS_MPI
    i = 0;
    while (i < static_cast<int>(ss_.bitsets_[BEND].size())) {
      if (!ss_.bitsets_[BEND][i]) {
        // check that this processor owns this bend
        if (bends_[i] != NULL) return bends_[i];
      }
      ++i;
    }
    return NULL;
#else
    i = ss_.bitsets_[BEND].firstOffBit();
    return i == -1 ? NULL : bends_[i];
#endif
  }

  Bend* SelectionManager::nextUnselectedBend(int& i) {
#ifdef IS_MPI
    ++i;
    while (i < static_cast<int>(ss_.bitsets_[BEND].size())) {
      if (!ss_.bitsets_[BEND][i]) {
        // check that this processor owns this bend
        if (bends_[i] != NULL) return bends_[i];
      }
      ++i;
    }
    return NULL;
#else
    i = ss_.bitsets_[BEND].nextOffBit(i);
    return i == -1 ? NULL : bends_[i];
#endif
  }

  Torsion* SelectionManager::beginSelectedTorsion(int& i) {
#ifdef IS_MPI
    i = 0;
    while (i < static_cast<int>(ss_.bitsets_[TORSION].size())) {
      if (ss_.bitsets_[TORSION][i]) {
        // check that this processor owns this torsion
        if (torsions_[i] != NULL) return torsions_[i];
      }
      ++i;
    }
    return NULL;
#else
    i = ss_.bitsets_[TORSION].firstOnBit();
    return i == -1 ? NULL : torsions_[i];
#endif
  }

  Torsion* SelectionManager::nextSelectedTorsion(int& i) {
#ifdef IS_MPI
    ++i;
    while (i < static_cast<int>(ss_.bitsets_[TORSION].size())) {
      if (ss_.bitsets_[TORSION][i]) {
        // check that this processor owns this torsion
        if (torsions_[i] != NULL) return torsions_[i];
      }
      ++i;
    }
    return NULL;
#else
    i = ss_.bitsets_[TORSION].nextOnBit(i);
    return i == -1 ? NULL : torsions_[i];
#endif
  }

  Torsion* SelectionManager::beginUnselectedTorsion(int& i) {
#ifdef IS_MPI
    i = 0;
    while (i < static_cast<int>(ss_.bitsets_[TORSION].size())) {
      if (!ss_.bitsets_[TORSION][i]) {
        // check that this processor owns this torsion
        if (torsions_[i] != NULL) return torsions_[i];
      }
      ++i;
    }
    return NULL;
#else
    i = ss_.bitsets_[TORSION].firstOffBit();
    return i == -1 ? NULL : torsions_[i];
#endif
  }

  Torsion* SelectionManager::nextUnselectedTorsion(int& i) {
#ifdef IS_MPI
    ++i;
    while (i < static_cast<int>(ss_.bitsets_[TORSION].size())) {
      if (!ss_.bitsets_[TORSION][i]) {
        // check that this processor owns this torsion
        if (torsions_[i] != NULL) return torsions_[i];
      }
      ++i;
    }
    return NULL;
#else
    i = ss_.bitsets_[TORSION].nextOffBit(i);
    return i == -1 ? NULL : torsions_[i];
#endif
  }

  Inversion* SelectionManager::beginSelectedInversion(int& i) {
#ifdef IS_MPI
    i = 0;
    while (i < static_cast<int>(ss_.bitsets_[INVERSION].size())) {
      if (ss_.bitsets_[INVERSION][i]) {
        // check that this processor owns this inversion
        if (inversions_[i] != NULL) return inversions_[i];
      }
      ++i;
    }
    return NULL;
#else
    i = ss_.bitsets_[INVERSION].firstOnBit();
    return i == -1 ? NULL : inversions_[i];
#endif
  }

  Inversion* SelectionManager::nextSelectedInversion(int& i) {
#ifdef IS_MPI
    ++i;
    while (i < static_cast<int>(ss_.bitsets_[INVERSION].size())) {
      if (ss_.bitsets_[INVERSION][i]) {
        // check that this processor owns this inversion
        if (inversions_[i] != NULL) return inversions_[i];
      }
      ++i;
    }
    return NULL;
#else
    i = ss_.bitsets_[INVERSION].nextOnBit(i);
    return i == -1 ? NULL : inversions_[i];
#endif
  }

  Inversion* SelectionManager::beginUnselectedInversion(int& i) {
#ifdef IS_MPI
    i = 0;
    while (i < static_cast<int>(ss_.bitsets_[INVERSION].size())) {
      if (!ss_.bitsets_[INVERSION][i]) {
        // check that this processor owns this inversion
        if (inversions_[i] != NULL) return inversions_[i];
      }
      ++i;
    }
    return NULL;
#else
    i = ss_.bitsets_[INVERSION].firstOffBit();
    return i == -1 ? NULL : inversions_[i];
#endif
  }

  Inversion* SelectionManager::nextUnselectedInversion(int& i) {
#ifdef IS_MPI
    ++i;
    while (i < static_cast<int>(ss_.bitsets_[INVERSION].size())) {
      if (!ss_.bitsets_[INVERSION][i]) {
        // check that this processor owns this inversion
        if (inversions_[i] != NULL) return inversions_[i];
      }
      ++i;
    }
    return NULL;
#else
    i = ss_.bitsets_[INVERSION].nextOffBit(i);
    return i == -1 ? NULL : inversions_[i];
#endif
  }

  Molecule* SelectionManager::beginSelectedMolecule(int& i) {
#ifdef IS_MPI
    i = 0;
    while (i < static_cast<int>(ss_.bitsets_[MOLECULE].size())) {
      if (ss_.bitsets_[MOLECULE][i]) {
        // check that this processor owns this molecule
        if (molecules_[i] != NULL) return molecules_[i];
      }
      ++i;
    }
    return NULL;
#else
    i = ss_.bitsets_[MOLECULE].firstOnBit();
    return i == -1 ? NULL : molecules_[i];
#endif
  }

  Molecule* SelectionManager::nextSelectedMolecule(int& i) {
#ifdef IS_MPI
    ++i;
    while (i < static_cast<int>(ss_.bitsets_[MOLECULE].size())) {
      if (ss_.bitsets_[MOLECULE][i]) {
        // check that this processor owns this molecule
        if (molecules_[i] != NULL) return molecules_[i];
      }
      ++i;
    }
    return NULL;
#else
    i = ss_.bitsets_[MOLECULE].nextOnBit(i);
    return i == -1 ? NULL : molecules_[i];
#endif
  }

  Molecule* SelectionManager::beginUnselectedMolecule(int& i) {
#ifdef IS_MPI
    i = 0;
    while (i < static_cast<int>(ss_.bitsets_[MOLECULE].size())) {
      if (!ss_.bitsets_[MOLECULE][i]) {
        // check that this processor owns this molecule
        if (molecules_[i] != NULL) return molecules_[i];
      }
      ++i;
    }
    return NULL;
#else
    i = ss_.bitsets_[MOLECULE].firstOffBit();
    return i == -1 ? NULL : molecules_[i];
#endif
  }

  Molecule* SelectionManager::nextUnselectedMolecule(int& i) {
#ifdef IS_MPI
    ++i;
    while (i < static_cast<int>(ss_.bitsets_[MOLECULE].size())) {
      if (!ss_.bitsets_[MOLECULE][i]) {
        // check that this processor owns this molecule
        if (molecules_[i] != NULL) return molecules_[i];
      }
      ++i;
    }
    return NULL;
#else
    i = ss_.bitsets_[MOLECULE].nextOffBit(i);
    return i == -1 ? NULL : molecules_[i];
#endif
  }

  Molecule* SelectionManager::nthSelectedMolecule(int& n) {
    int i;
#ifdef IS_MPI
    i = ss_.bitsets_[MOLECULE].nthOnBit(n);
    if (i == -1) return NULL;
    // check that this processor owns this molecule
    if (molecules_[i] != NULL) return molecules_[i];
    return NULL;
#else
    i = ss_.bitsets_[MOLECULE].nthOnBit(n);
    return i == -1 ? NULL : molecules_[i];
#endif
  }

  /**
   * getSelectedAtomTypes
   *
   * Returns an STL set of AtomType* that are actually selected.
   * Must query all processors to assemble this information.
   *
   */
  AtomTypeSet SelectionManager::getSelectedAtomTypes() {
    AtomTypeSet atomTypes;

    for (size_t i = 0; i < ss_.bitsets_[STUNTDOUBLE].size(); ++i) {
      if (ss_.bitsets_[STUNTDOUBLE][i]) {
        // check that this processor owns this stuntdouble
        if (stuntdoubles_[i] != NULL) {
          if (stuntdoubles_[i]->isAtom()) {
            Atom* atom = static_cast<Atom*>(stuntdoubles_[i]);
            atomTypes.insert(atom->getAtomType());
          }
        }
      }
    }

#ifdef IS_MPI
    // loop over the found atom types on this processor, and add their
    // numerical idents to a vector:

    std::vector<int> foundTypes;
    AtomTypeSet::iterator i;
    for (i = atomTypes.begin(); i != atomTypes.end(); ++i)
      foundTypes.push_back((*i)->getIdent());

    // count_local holds the number of found types on this processor
    int count_local = foundTypes.size();

    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    // we need arrays to hold the counts and displacement vectors for
    // all processors
    std::vector<int> counts(nproc, 0);
    std::vector<int> disps(nproc, 0);

    // fill the counts array
    MPI_Allgather(&count_local, 1, MPI_INT, &counts[0], 1, MPI_INT,
                  MPI_COMM_WORLD);

    // use the processor counts to compute the displacement array
    disps[0]       = 0;
    int totalCount = counts[0];
    for (int iproc = 1; iproc < nproc; iproc++) {
      disps[iproc] = disps[iproc - 1] + counts[iproc - 1];
      totalCount += counts[iproc];
    }

    if (totalCount > 0) {
      // we need a (possibly redundant) set of all found types:
      std::vector<int> ftGlobal(totalCount);

      // now spray out the foundTypes to all the other processors:
      MPI_Allgatherv(&foundTypes[0], count_local, MPI_INT, &ftGlobal[0],
                     &counts[0], &disps[0], MPI_INT, MPI_COMM_WORLD);

      std::vector<int>::iterator j;

      // foundIdents is a stl set, so inserting an already found ident
      // will have no effect.
      std::set<int> foundIdents;

      for (j = ftGlobal.begin(); j != ftGlobal.end(); ++j)
        foundIdents.insert((*j));

      // now iterate over the foundIdents and get the actual atom types
      // that correspond to these:
      ForceField* forceField_ = info_->getForceField();
      std::set<int>::iterator it;
      for (it = foundIdents.begin(); it != foundIdents.end(); ++it)
        atomTypes.insert(forceField_->getAtomType((*it)));
    }
#endif

    return atomTypes;
  }

  MoleculeStampSet SelectionManager::getSelectedMoleculeStamps() {
    MoleculeStampSet moleculeStamps;

    for (size_t i = 0; i < ss_.bitsets_[MOLECULE].size(); ++i) {
      if (ss_.bitsets_[MOLECULE][i]) {
        // check that this processor owns this molecule
        if (molecules_[i] != NULL) {
          Molecule* mol = static_cast<Molecule*>(molecules_[i]);
          moleculeStamps.insert(mol->getMolStamp());
        }
      }
    }

#ifdef IS_MPI
    std::vector<int> foundStamps;
    MoleculeStampSet::iterator i;
    for (i = moleculeStamps.begin(); i != moleculeStamps.end(); ++i)
      foundStamps.push_back((*i)->getIdent());

    // count_local holds the number of found stamps on this processor
    int count_local = foundStamps.size();

    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    // we need arrays to hold the counts and displacement vectors for
    // all processors
    std::vector<int> counts(nproc, 0);
    std::vector<int> disps(nproc, 0);

    // fill the counts array
    MPI_Allgather(&count_local, 1, MPI_INT, &counts[0], 1, MPI_INT,
                  MPI_COMM_WORLD);

    // use the processor counts to compute the displacement array
    disps[0]       = 0;
    int totalCount = counts[0];
    for (int iproc = 1; iproc < nproc; iproc++) {
      disps[iproc] = disps[iproc - 1] + counts[iproc - 1];
      totalCount += counts[iproc];
    }

    if (totalCount > 0) {
      // we need a (possibly redundant) set of all found stamps:
      std::vector<int> fsGlobal(totalCount);

      // now spray out the foundStamps to all the other processors:
      MPI_Allgatherv(&foundStamps[0], count_local, MPI_INT, &fsGlobal[0],
                     &counts[0], &disps[0], MPI_INT, MPI_COMM_WORLD);

      std::vector<int>::iterator j;

      for (j = fsGlobal.begin(); j != fsGlobal.end(); ++j)
        moleculeStamps.insert(info_->getMoleculeStamp(*j));
    }
#endif

    return moleculeStamps;
  }

  SelectionManager SelectionManager::replaceRigidBodiesWithAtoms() const {
    SelectionSet ssAtoms(nObjects_);
    ssAtoms.clearAll();

    SelectionSet ssRBs(nObjects_);
    ssRBs.clearAll();

    SelectionManager tempSeleMan = *this;

    StuntDouble* sd;
    int isd;
    std::vector<Atom*>::iterator ai;
    Atom* atom;
    for (sd = tempSeleMan.beginSelected(isd); sd != NULL;
         sd = tempSeleMan.nextSelected(isd)) {
      if (sd->isRigidBody()) {
        RigidBody* rb = static_cast<RigidBody*>(sd);
        for (atom = rb->beginAtom(ai); atom != NULL; atom = rb->nextAtom(ai)) {
          ssAtoms.bitsets_[STUNTDOUBLE].setBitOn(atom->getGlobalIndex());
        }

        ssRBs.bitsets_[STUNTDOUBLE].setBitOn(rb->getGlobalIndex());
      }
    }

    // Add the atoms in rigid bodies and remove those rigid bodies from our
    // selection set.
    tempSeleMan.ss_ |= ssAtoms;
    tempSeleMan.ss_ -= ssRBs;

    return tempSeleMan;
  }

  SelectionManager SelectionManager::removeAtomsInRigidBodies() const {
    SelectionSet ssAtoms(nObjects_);
    ssAtoms.clearAll();

    SelectionManager tempSeleMan = *this;

    StuntDouble* sd;
    int isd;
    std::vector<Atom*>::iterator ai;
    Atom* atom;
    for (sd = tempSeleMan.beginSelected(isd); sd != NULL;
         sd = tempSeleMan.nextSelected(isd)) {
      if (sd->isRigidBody()) {
        RigidBody* rb = static_cast<RigidBody*>(sd);
        for (atom = rb->beginAtom(ai); atom != NULL; atom = rb->nextAtom(ai)) {
          ssAtoms.bitsets_[STUNTDOUBLE].setBitOn(atom->getGlobalIndex());
        }
      }
    }

    // Remove the atoms in rigid bodies from our selection set.
    tempSeleMan.ss_ -= ssAtoms;

    return tempSeleMan;
  }

  SelectionManager operator|(const SelectionManager& sman1,
                             const SelectionManager& sman2) {
    SelectionManager result(sman1);
    result |= sman2;
    return result;
  }

  SelectionManager operator&(const SelectionManager& sman1,
                             const SelectionManager& sman2) {
    SelectionManager result(sman1);
    result &= sman2;
    return result;
  }

  SelectionManager operator^(const SelectionManager& sman1,
                             const SelectionManager& sman2) {
    SelectionManager result(sman1);
    result ^= sman2;
    return result;
  }

  SelectionManager operator-(const SelectionManager& sman1,
                             const SelectionManager& sman2) {
    SelectionManager result(sman1);
    result -= sman2;
    return result;
  }

}  // namespace OpenMD
