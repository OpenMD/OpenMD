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
 * @file SimInfo.cpp
 * @author    tlin
 * @date  11/02/2004
 * @version 1.0
 */

#include "brains/SimInfo.hpp"

#include <algorithm>
#include <cstdint>
#include <map>
#include <memory>
#include <random>
#include <set>

#ifdef IS_MPI
#include <mpi.h>
#endif

#include "brains/ForceField.hpp"
#include "io/ForceFieldOptions.hpp"
#include "math/Vector3.hpp"
#include "nonbonded/SwitchingFunction.hpp"
#include "primitives/Molecule.hpp"
#include "primitives/StuntDouble.hpp"
#include "selection/SelectionManager.hpp"
#include "utils/MemoryUtils.hpp"
#include "utils/RandNumGen.hpp"
#include "utils/simError.h"

using namespace std;
namespace OpenMD {

  SimInfo::SimInfo(ForceField* ff, Globals* simParams) :
      forceField_(ff), simParams_(simParams), randNumGen_ {nullptr}, nAtoms_(0),
      nBonds_(0), nBends_(0), nTorsions_(0), nInversions_(0), nRigidBodies_(0),
      nIntegrableObjects_(0), nCutoffGroups_(0), nConstraints_(0),
      nFluctuatingCharges_(0), nGlobalMols_(0), nGlobalAtoms_(0),
      nGlobalCutoffGroups_(0), nGlobalIntegrableObjects_(0),
      nGlobalRigidBodies_(0), nGlobalFluctuatingCharges_(0), nGlobalBonds_(0),
      nGlobalBends_(0), nGlobalTorsions_(0), nGlobalInversions_(0),
      nGlobalConstraints_(0), hasNGlobalConstraints_(false), ndf_(0),
      fdf_local(0), ndfRaw_(0), ndfTrans_(0), nZconstraint_(0), sman_(NULL),
      topologyDone_(false), calcBoxDipole_(false), calcBoxQuadrupole_(false),
      useAtomicVirial_(true) {
    MoleculeStamp* molStamp;
    int nMolWithSameStamp;
    int nCutoffAtoms = 0;  // number of atoms belong to cutoff groups
    int nGroups      = 0;  // total cutoff groups defined in meta-data file
    CutoffGroupStamp* cgStamp;
    RigidBodyStamp* rbStamp;
    int nRigidAtoms = 0;

    vector<Component*> components = simParams->getComponents();

    for (vector<Component*>::iterator i = components.begin();
         i != components.end(); ++i) {
      molStamp = (*i)->getMoleculeStamp();
      if ((*i)->haveRegion()) {
        molStamp->setRegion((*i)->getRegion());
      } else {
        // set the region to a disallowed value:
        molStamp->setRegion(-1);
      }

      nMolWithSameStamp = (*i)->getNMol();

      addMoleculeStamp(molStamp, nMolWithSameStamp);

      // calculate atoms in molecules
      nGlobalAtoms_ += molStamp->getNAtoms() * nMolWithSameStamp;
      nGlobalBonds_ += molStamp->getNBonds() * nMolWithSameStamp;
      nGlobalBends_ += molStamp->getNBends() * nMolWithSameStamp;
      nGlobalTorsions_ += molStamp->getNTorsions() * nMolWithSameStamp;
      nGlobalInversions_ += molStamp->getNInversions() * nMolWithSameStamp;

      // calculate atoms in cutoff groups
      int nAtomsInGroups       = 0;
      int nCutoffGroupsInStamp = molStamp->getNCutoffGroups();

      for (int j = 0; j < nCutoffGroupsInStamp; j++) {
        cgStamp = molStamp->getCutoffGroupStamp(j);
        nAtomsInGroups += cgStamp->getNMembers();
      }

      nGroups += nCutoffGroupsInStamp * nMolWithSameStamp;

      nCutoffAtoms += nAtomsInGroups * nMolWithSameStamp;

      // calculate atoms in rigid bodies
      int nAtomsInRigidBodies = 0;
      int nRigidBodiesInStamp = molStamp->getNRigidBodies();

      for (int j = 0; j < nRigidBodiesInStamp; j++) {
        rbStamp = molStamp->getRigidBodyStamp(j);
        nAtomsInRigidBodies += rbStamp->getNMembers();
      }

      nGlobalRigidBodies_ += nRigidBodiesInStamp * nMolWithSameStamp;
      nRigidAtoms += nAtomsInRigidBodies * nMolWithSameStamp;
    }

    // every free atom (atom does not belong to cutoff groups) is a cutoff
    // group therefore the total number of cutoff groups in the system is
    // equal to the total number of atoms minus number of atoms belong to
    // cutoff group defined in meta-data file plus the number of cutoff
    // groups defined in meta-data file

    nGlobalCutoffGroups_ = nGlobalAtoms_ - nCutoffAtoms + nGroups;

    // every free atom (atom does not belong to rigid bodies) is an
    // integrable object therefore the total number of integrable objects
    // in the system is equal to the total number of atoms minus number of
    // atoms belong to rigid body defined in meta-data file plus the number
    // of rigid bodies defined in meta-data file
    nGlobalIntegrableObjects_ =
        nGlobalAtoms_ - nRigidAtoms + nGlobalRigidBodies_;

    nGlobalMols_ = molStampIds_.size();
    molToProcMap_.resize(nGlobalMols_);

    // Initialize the random number generator on each processing
    // element using either the user-supplied seed or a default seed
    std::uint_fast32_t seed;

    if (simParams_->haveSeed())
      seed = static_cast<std::uint_fast32_t>(simParams_->getSeed());
    else
      seed = std::mt19937::default_seed;

    randNumGen_ = std::make_shared<Utils::RandNumGen>(seed);
  }

  SimInfo::~SimInfo() {
    Utils::deletePointers(molecules_);

    delete sman_;
    delete simParams_;
    delete forceField_;
  }

  bool SimInfo::addMolecule(Molecule* mol) {
    MoleculeIterator i;

    i = molecules_.find(mol->getGlobalIndex());
    if (i == molecules_.end()) {
      molecules_.insert(make_pair(mol->getGlobalIndex(), mol));

      nAtoms_ += mol->getNAtoms();
      nBonds_ += mol->getNBonds();
      nBends_ += mol->getNBends();
      nTorsions_ += mol->getNTorsions();
      nInversions_ += mol->getNInversions();
      nRigidBodies_ += mol->getNRigidBodies();
      nIntegrableObjects_ += mol->getNIntegrableObjects();
      nCutoffGroups_ += mol->getNCutoffGroups();
      nConstraints_ += mol->getNConstraintPairs();

      addInteractionPairs(mol);

      return true;
    } else {
      return false;
    }
  }

  bool SimInfo::removeMolecule(Molecule* mol) {
    MoleculeIterator i;
    i = molecules_.find(mol->getGlobalIndex());

    if (i != molecules_.end()) {
      assert(mol == i->second);

      nAtoms_ -= mol->getNAtoms();
      nBonds_ -= mol->getNBonds();
      nBends_ -= mol->getNBends();
      nTorsions_ -= mol->getNTorsions();
      nInversions_ -= mol->getNInversions();
      nRigidBodies_ -= mol->getNRigidBodies();
      nIntegrableObjects_ -= mol->getNIntegrableObjects();
      nCutoffGroups_ -= mol->getNCutoffGroups();
      nConstraints_ -= mol->getNConstraintPairs();

      removeInteractionPairs(mol);
      molecules_.erase(mol->getGlobalIndex());

      delete mol;

      return true;
    } else {
      return false;
    }
  }

  Molecule* SimInfo::beginMolecule(MoleculeIterator& i) {
    i = molecules_.begin();
    return i == molecules_.end() ? NULL : i->second;
  }

  Molecule* SimInfo::nextMolecule(MoleculeIterator& i) {
    ++i;
    return i == molecules_.end() ? NULL : i->second;
  }

  void SimInfo::calcNdf() {
    int ndf_local, nfq_local;
    MoleculeIterator i;
    vector<StuntDouble*>::iterator j;
    vector<Atom*>::iterator k;

    Molecule* mol;
    StuntDouble* sd;
    Atom* atom;

    ndf_local = 0;
    nfq_local = 0;

    for (mol = beginMolecule(i); mol != NULL; mol = nextMolecule(i)) {
      for (sd = mol->beginIntegrableObject(j); sd != NULL;
           sd = mol->nextIntegrableObject(j)) {
        ndf_local += 3;

        if (sd->isDirectional()) {
          if (sd->isLinear()) {
            ndf_local += 2;
          } else {
            ndf_local += 3;
          }
        }
      }

      for (atom = mol->beginFluctuatingCharge(k); atom != NULL;
           atom = mol->nextFluctuatingCharge(k)) {
        if (atom->isFluctuatingCharge()) { nfq_local++; }
      }
    }

    ndfLocal_ = ndf_local;

    // n_constraints is local, so subtract them on each processor
    ndf_local -= nConstraints_;

#ifdef IS_MPI
    MPI_Allreduce(&ndf_local, &ndf_, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&nfq_local, &nGlobalFluctuatingCharges_, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
#else
    ndf_                       = ndf_local;
    nGlobalFluctuatingCharges_ = nfq_local;
#endif

    // nZconstraints_ is global, as are the 3 COM translations for the
    // entire system:
    ndf_ = ndf_ - 3 - nZconstraint_;
  }

  int SimInfo::getFdf() {
#ifdef IS_MPI
    MPI_Allreduce(&fdf_local, &fdf_, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
    fdf_                       = fdf_local;
#endif
    return fdf_;
  }

  unsigned int SimInfo::getNLocalCutoffGroups() {
    int nLocalCutoffAtoms = 0;
    Molecule* mol;
    MoleculeIterator mi;
    CutoffGroup* cg;
    Molecule::CutoffGroupIterator ci;

    for (mol = beginMolecule(mi); mol != NULL; mol = nextMolecule(mi)) {
      for (cg = mol->beginCutoffGroup(ci); cg != NULL;
           cg = mol->nextCutoffGroup(ci)) {
        nLocalCutoffAtoms += cg->getNumAtom();
      }
    }

    return nAtoms_ - nLocalCutoffAtoms + nCutoffGroups_;
  }

  void SimInfo::calcNdfRaw() {
    int ndfRaw_local;

    MoleculeIterator i;
    vector<StuntDouble*>::iterator j;
    Molecule* mol;
    StuntDouble* sd;

    // Raw degrees of freedom that we have to set
    ndfRaw_local = 0;

    for (mol = beginMolecule(i); mol != NULL; mol = nextMolecule(i)) {
      for (sd = mol->beginIntegrableObject(j); sd != NULL;
           sd = mol->nextIntegrableObject(j)) {
        ndfRaw_local += 3;

        if (sd->isDirectional()) {
          if (sd->isLinear()) {
            ndfRaw_local += 2;
          } else {
            ndfRaw_local += 3;
          }
        }
      }
    }

#ifdef IS_MPI
    MPI_Allreduce(&ndfRaw_local, &ndfRaw_, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
    ndfRaw_                    = ndfRaw_local;
#endif
  }

  void SimInfo::calcNdfTrans() {
    int ndfTrans_local;

    ndfTrans_local = 3 * nIntegrableObjects_ - nConstraints_;

#ifdef IS_MPI
    MPI_Allreduce(&ndfTrans_local, &ndfTrans_, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
#else
    ndfTrans_                  = ndfTrans_local;
#endif

    ndfTrans_ = ndfTrans_ - 3 - nZconstraint_;
  }

  void SimInfo::addInteractionPairs(Molecule* mol) {
    ForceFieldOptions& options_ = forceField_->getForceFieldOptions();
    vector<Bond*>::iterator bondIter;
    vector<Bend*>::iterator bendIter;
    vector<Torsion*>::iterator torsionIter;
    vector<Inversion*>::iterator inversionIter;
    Bond* bond;
    Bend* bend;
    Torsion* torsion;
    Inversion* inversion;
    int a;
    int b;
    int c;
    int d;

    // atomGroups can be used to add special interaction maps between
    // groups of atoms that are in two separate rigid bodies.
    // However, most site-site interactions between two rigid bodies
    // are probably not special, just the ones between the physically
    // bonded atoms.  Interactions *within* a single rigid body should
    // always be excluded.  These are done at the bottom of this
    // function.

    map<int, set<int>> atomGroups;
    Molecule::RigidBodyIterator rbIter;
    RigidBody* rb;
    Molecule::IntegrableObjectIterator ii;
    StuntDouble* sd;

    for (sd = mol->beginIntegrableObject(ii); sd != NULL;
         sd = mol->nextIntegrableObject(ii)) {
      if (sd->isRigidBody()) {
        rb                  = static_cast<RigidBody*>(sd);
        vector<Atom*> atoms = rb->getAtoms();
        set<int> rigidAtoms;
        for (int i = 0; i < static_cast<int>(atoms.size()); ++i) {
          rigidAtoms.insert(atoms[i]->getGlobalIndex());
        }
        for (int i = 0; i < static_cast<int>(atoms.size()); ++i) {
          atomGroups.insert(map<int, set<int>>::value_type(
              atoms[i]->getGlobalIndex(), rigidAtoms));
        }
      } else {
        set<int> oneAtomSet;
        oneAtomSet.insert(sd->getGlobalIndex());
        atomGroups.insert(
            map<int, set<int>>::value_type(sd->getGlobalIndex(), oneAtomSet));
      }
    }

    for (bond = mol->beginBond(bondIter); bond != NULL;
         bond = mol->nextBond(bondIter)) {
      a = bond->getAtomA()->getGlobalIndex();
      b = bond->getAtomB()->getGlobalIndex();

      if (options_.havevdw12scale() || options_.haveelectrostatic12scale()) {
        oneTwoInteractions_.addPair(a, b);
      } else {
        excludedInteractions_.addPair(a, b);
      }
    }

    for (bend = mol->beginBend(bendIter); bend != NULL;
         bend = mol->nextBend(bendIter)) {
      a = bend->getAtomA()->getGlobalIndex();
      b = bend->getAtomB()->getGlobalIndex();
      c = bend->getAtomC()->getGlobalIndex();

      if (options_.havevdw12scale() || options_.haveelectrostatic12scale()) {
        oneTwoInteractions_.addPair(a, b);
        oneTwoInteractions_.addPair(b, c);
      } else {
        excludedInteractions_.addPair(a, b);
        excludedInteractions_.addPair(b, c);
      }

      if (options_.havevdw13scale() || options_.haveelectrostatic13scale()) {
        oneThreeInteractions_.addPair(a, c);
      } else {
        excludedInteractions_.addPair(a, c);
      }
    }

    for (torsion = mol->beginTorsion(torsionIter); torsion != NULL;
         torsion = mol->nextTorsion(torsionIter)) {
      a = torsion->getAtomA()->getGlobalIndex();
      b = torsion->getAtomB()->getGlobalIndex();
      c = torsion->getAtomC()->getGlobalIndex();
      d = torsion->getAtomD()->getGlobalIndex();

      if (options_.havevdw12scale() || options_.haveelectrostatic12scale()) {
        oneTwoInteractions_.addPair(a, b);
        oneTwoInteractions_.addPair(b, c);
        oneTwoInteractions_.addPair(c, d);
      } else {
        excludedInteractions_.addPair(a, b);
        excludedInteractions_.addPair(b, c);
        excludedInteractions_.addPair(c, d);
      }

      if (options_.havevdw13scale() || options_.haveelectrostatic13scale()) {
        oneThreeInteractions_.addPair(a, c);
        oneThreeInteractions_.addPair(b, d);
      } else {
        excludedInteractions_.addPair(a, c);
        excludedInteractions_.addPair(b, d);
      }

      if (options_.havevdw14scale() || options_.haveelectrostatic14scale()) {
        oneFourInteractions_.addPair(a, d);
      } else {
        excludedInteractions_.addPair(a, d);
      }
    }

    for (inversion = mol->beginInversion(inversionIter); inversion != NULL;
         inversion = mol->nextInversion(inversionIter)) {
      a = inversion->getAtomA()->getGlobalIndex();
      b = inversion->getAtomB()->getGlobalIndex();
      c = inversion->getAtomC()->getGlobalIndex();
      d = inversion->getAtomD()->getGlobalIndex();

      if (options_.havevdw12scale() || options_.haveelectrostatic12scale()) {
        oneTwoInteractions_.addPair(a, b);
        oneTwoInteractions_.addPair(a, c);
        oneTwoInteractions_.addPair(a, d);
      } else {
        excludedInteractions_.addPair(a, b);
        excludedInteractions_.addPair(a, c);
        excludedInteractions_.addPair(a, d);
      }

      if (options_.havevdw13scale() || options_.haveelectrostatic13scale()) {
        oneThreeInteractions_.addPair(b, c);
        oneThreeInteractions_.addPair(b, d);
        oneThreeInteractions_.addPair(c, d);
      } else {
        excludedInteractions_.addPair(b, c);
        excludedInteractions_.addPair(b, d);
        excludedInteractions_.addPair(c, d);
      }
    }

    for (rb = mol->beginRigidBody(rbIter); rb != NULL;
         rb = mol->nextRigidBody(rbIter)) {
      vector<Atom*> atoms = rb->getAtoms();
      for (int i = 0; i < static_cast<int>(atoms.size()) - 1; ++i) {
        for (int j = i + 1; j < static_cast<int>(atoms.size()); ++j) {
          a = atoms[i]->getGlobalIndex();
          b = atoms[j]->getGlobalIndex();
          excludedInteractions_.addPair(a, b);
        }
      }
    }
  }

  void SimInfo::removeInteractionPairs(Molecule* mol) {
    ForceFieldOptions& options_ = forceField_->getForceFieldOptions();
    vector<Bond*>::iterator bondIter;
    vector<Bend*>::iterator bendIter;
    vector<Torsion*>::iterator torsionIter;
    vector<Inversion*>::iterator inversionIter;
    Bond* bond;
    Bend* bend;
    Torsion* torsion;
    Inversion* inversion;
    int a;
    int b;
    int c;
    int d;

    map<int, set<int>> atomGroups;
    Molecule::RigidBodyIterator rbIter;
    RigidBody* rb;
    Molecule::IntegrableObjectIterator ii;
    StuntDouble* sd;

    for (sd = mol->beginIntegrableObject(ii); sd != NULL;
         sd = mol->nextIntegrableObject(ii)) {
      if (sd->isRigidBody()) {
        rb                  = static_cast<RigidBody*>(sd);
        vector<Atom*> atoms = rb->getAtoms();
        set<int> rigidAtoms;
        for (int i = 0; i < static_cast<int>(atoms.size()); ++i) {
          rigidAtoms.insert(atoms[i]->getGlobalIndex());
        }
        for (int i = 0; i < static_cast<int>(atoms.size()); ++i) {
          atomGroups.insert(map<int, set<int>>::value_type(
              atoms[i]->getGlobalIndex(), rigidAtoms));
        }
      } else {
        set<int> oneAtomSet;
        oneAtomSet.insert(sd->getGlobalIndex());
        atomGroups.insert(
            map<int, set<int>>::value_type(sd->getGlobalIndex(), oneAtomSet));
      }
    }

    for (bond = mol->beginBond(bondIter); bond != NULL;
         bond = mol->nextBond(bondIter)) {
      a = bond->getAtomA()->getGlobalIndex();
      b = bond->getAtomB()->getGlobalIndex();

      if (options_.havevdw12scale() || options_.haveelectrostatic12scale()) {
        oneTwoInteractions_.removePair(a, b);
      } else {
        excludedInteractions_.removePair(a, b);
      }
    }

    for (bend = mol->beginBend(bendIter); bend != NULL;
         bend = mol->nextBend(bendIter)) {
      a = bend->getAtomA()->getGlobalIndex();
      b = bend->getAtomB()->getGlobalIndex();
      c = bend->getAtomC()->getGlobalIndex();

      if (options_.havevdw12scale() || options_.haveelectrostatic12scale()) {
        oneTwoInteractions_.removePair(a, b);
        oneTwoInteractions_.removePair(b, c);
      } else {
        excludedInteractions_.removePair(a, b);
        excludedInteractions_.removePair(b, c);
      }

      if (options_.havevdw13scale() || options_.haveelectrostatic13scale()) {
        oneThreeInteractions_.removePair(a, c);
      } else {
        excludedInteractions_.removePair(a, c);
      }
    }

    for (torsion = mol->beginTorsion(torsionIter); torsion != NULL;
         torsion = mol->nextTorsion(torsionIter)) {
      a = torsion->getAtomA()->getGlobalIndex();
      b = torsion->getAtomB()->getGlobalIndex();
      c = torsion->getAtomC()->getGlobalIndex();
      d = torsion->getAtomD()->getGlobalIndex();

      if (options_.havevdw12scale() || options_.haveelectrostatic12scale()) {
        oneTwoInteractions_.removePair(a, b);
        oneTwoInteractions_.removePair(b, c);
        oneTwoInteractions_.removePair(c, d);
      } else {
        excludedInteractions_.removePair(a, b);
        excludedInteractions_.removePair(b, c);
        excludedInteractions_.removePair(c, d);
      }

      if (options_.havevdw13scale() || options_.haveelectrostatic13scale()) {
        oneThreeInteractions_.removePair(a, c);
        oneThreeInteractions_.removePair(b, d);
      } else {
        excludedInteractions_.removePair(a, c);
        excludedInteractions_.removePair(b, d);
      }

      if (options_.havevdw14scale() || options_.haveelectrostatic14scale()) {
        oneFourInteractions_.removePair(a, d);
      } else {
        excludedInteractions_.removePair(a, d);
      }
    }

    for (inversion = mol->beginInversion(inversionIter); inversion != NULL;
         inversion = mol->nextInversion(inversionIter)) {
      a = inversion->getAtomA()->getGlobalIndex();
      b = inversion->getAtomB()->getGlobalIndex();
      c = inversion->getAtomC()->getGlobalIndex();
      d = inversion->getAtomD()->getGlobalIndex();

      if (options_.havevdw12scale() || options_.haveelectrostatic12scale()) {
        oneTwoInteractions_.removePair(a, b);
        oneTwoInteractions_.removePair(a, c);
        oneTwoInteractions_.removePair(a, d);
      } else {
        excludedInteractions_.removePair(a, b);
        excludedInteractions_.removePair(a, c);
        excludedInteractions_.removePair(a, d);
      }

      if (options_.havevdw13scale() || options_.haveelectrostatic13scale()) {
        oneThreeInteractions_.removePair(b, c);
        oneThreeInteractions_.removePair(b, d);
        oneThreeInteractions_.removePair(c, d);
      } else {
        excludedInteractions_.removePair(b, c);
        excludedInteractions_.removePair(b, d);
        excludedInteractions_.removePair(c, d);
      }
    }

    for (rb = mol->beginRigidBody(rbIter); rb != NULL;
         rb = mol->nextRigidBody(rbIter)) {
      vector<Atom*> atoms = rb->getAtoms();
      for (int i = 0; i < static_cast<int>(atoms.size()) - 1; ++i) {
        for (int j = i + 1; j < static_cast<int>(atoms.size()); ++j) {
          a = atoms[i]->getGlobalIndex();
          b = atoms[j]->getGlobalIndex();
          excludedInteractions_.removePair(a, b);
        }
      }
    }
  }

  void SimInfo::addMoleculeStamp(MoleculeStamp* molStamp, int nmol) {
    int curStampId;

    // index from 0
    curStampId = moleculeStamps_.size();

    moleculeStamps_.push_back(molStamp);
    moleculeStamps_[moleculeStamps_.size() - 1]->setIdent(curStampId);
    molStampIds_.insert(molStampIds_.end(), nmol, curStampId);
  }

  /**
   * update
   *
   *  Performs the global checks and variable settings after the
   *  objects have been created.
   *
   */
  void SimInfo::update() {
    setupSimVariables();
    calcNConstraints();
    calcNdf();
    calcNdfRaw();
    calcNdfTrans();
  }

  /**
   * getSimulatedAtomTypes
   *
   * Returns an STL set of AtomType* that are actually present in this
   * simulation.  Must query all processors to assemble this information.
   *
   */
  AtomTypeSet SimInfo::getSimulatedAtomTypes() {
    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::AtomIterator ai;
    Atom* atom;
    AtomTypeSet atomTypes;

    for (mol = beginMolecule(mi); mol != NULL; mol = nextMolecule(mi)) {
      for (atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
        atomTypes.insert(atom->getAtomType());
      }
    }

#ifdef IS_MPI

    // loop over the found atom types on this processor, and add their
    // numerical idents to a vector:

    vector<int> foundTypes;
    AtomTypeSet::iterator i;
    for (i = atomTypes.begin(); i != atomTypes.end(); ++i)
      foundTypes.push_back((*i)->getIdent());

    // count_local holds the number of found types on this processor
    int count_local = foundTypes.size();

    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    // we need arrays to hold the counts and displacement vectors for
    // all processors
    vector<int> counts(nproc, 0);
    vector<int> disps(nproc, 0);

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

    // we need a (possibly redundant) set of all found types:
    vector<int> ftGlobal(totalCount);

    // now spray out the foundTypes to all the other processors:
    MPI_Allgatherv(&foundTypes[0], count_local, MPI_INT, &ftGlobal[0],
                   &counts[0], &disps[0], MPI_INT, MPI_COMM_WORLD);

    vector<int>::iterator j;

    // foundIdents is a stl set, so inserting an already found ident
    // will have no effect.
    set<int> foundIdents;

    for (j = ftGlobal.begin(); j != ftGlobal.end(); ++j)
      foundIdents.insert((*j));

    // now iterate over the foundIdents and get the actual atom types
    // that correspond to these:
    set<int>::iterator it;
    for (it = foundIdents.begin(); it != foundIdents.end(); ++it)
      atomTypes.insert(forceField_->getAtomType((*it)));

#endif

    return atomTypes;
  }

  int getGlobalCountOfType(AtomType*) {
    /*
    AtomTypeSet atypes = getSimulatedAtomTypes();
    map<AtomType*, int> counts_;

    for(mol = beginMolecule(mi); mol != NULL; mol = nextMolecule(mi)) {
    for(atom = mol->beginAtom(ai); atom != NULL;
    atom = mol->nextAtom(ai)) {
    atom->getAtomType();
    }
    }
  */
    return 0;
  }

  void SimInfo::setupSimVariables() {
    useAtomicVirial_ = simParams_->getUseAtomicVirial();
    // we only call setAccumulateBoxDipole if the accumulateBoxDipole
    // parameter is true
    calcBoxDipole_ = false;
    if (simParams_->haveAccumulateBoxDipole())
      if (simParams_->getAccumulateBoxDipole()) { calcBoxDipole_ = true; }
    // we only call setAccumulateBoxQuadrupole if the accumulateBoxQuadrupole
    // parameter is true
    calcBoxQuadrupole_ = false;
    if (simParams_->haveAccumulateBoxQuadrupole())
      if (simParams_->getAccumulateBoxQuadrupole()) {
        calcBoxQuadrupole_ = true;
      }

    AtomTypeSet::iterator i;
    AtomTypeSet atomTypes;
    atomTypes                   = getSimulatedAtomTypes();
    bool usesElectrostatic      = false;
    bool usesMetallic           = false;
    bool usesDirectional        = false;
    bool usesFluctuatingCharges = false;
    // loop over all of the atom types
    for (i = atomTypes.begin(); i != atomTypes.end(); ++i) {
      usesElectrostatic |= (*i)->isElectrostatic();
      usesMetallic |= (*i)->isMetal();
      usesDirectional |= (*i)->isDirectional();
      usesFluctuatingCharges |= (*i)->isFluctuatingCharge();
    }

#ifdef IS_MPI
    int temp;

    temp = usesDirectional;
    MPI_Allreduce(MPI_IN_PLACE, &temp, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
    usesDirectionalAtoms_ = (temp == 0) ? false : true;

    temp = usesMetallic;
    MPI_Allreduce(MPI_IN_PLACE, &temp, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
    usesMetallicAtoms_ = (temp == 0) ? false : true;

    temp = usesElectrostatic;
    MPI_Allreduce(MPI_IN_PLACE, &temp, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
    usesElectrostaticAtoms_ = (temp == 0) ? false : true;

    temp = usesFluctuatingCharges;
    MPI_Allreduce(MPI_IN_PLACE, &temp, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
    usesFluctuatingCharges_ = (temp == 0) ? false : true;
#else

    usesDirectionalAtoms_   = usesDirectional;
    usesMetallicAtoms_      = usesMetallic;
    usesElectrostaticAtoms_ = usesElectrostatic;
    usesFluctuatingCharges_ = usesFluctuatingCharges;

#endif

    requiresPrepair_        = usesMetallicAtoms_ ? true : false;
    requiresSkipCorrection_ = usesElectrostaticAtoms_ ? true : false;
    requiresSelfCorrection_ = usesElectrostaticAtoms_ ? true : false;
  }

  vector<int> SimInfo::getGlobalAtomIndices() {
    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::AtomIterator ai;
    Atom* atom;

    vector<int> GlobalAtomIndices(getNAtoms(), 0);

    for (mol = beginMolecule(mi); mol != NULL; mol = nextMolecule(mi)) {
      for (atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
        GlobalAtomIndices[atom->getLocalIndex()] = atom->getGlobalIndex();
      }
    }
    return GlobalAtomIndices;
  }

  vector<int> SimInfo::getGlobalGroupIndices() {
    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::CutoffGroupIterator ci;
    CutoffGroup* cg;

    vector<int> GlobalGroupIndices;

    for (mol = beginMolecule(mi); mol != NULL; mol = nextMolecule(mi)) {
      // local index of cutoff group is trivial, it only depends on the
      // order of travesing
      for (cg = mol->beginCutoffGroup(ci); cg != NULL;
           cg = mol->nextCutoffGroup(ci)) {
        GlobalGroupIndices.push_back(cg->getGlobalIndex());
      }
    }
    return GlobalGroupIndices;
  }

  void SimInfo::prepareTopology() {
    // calculate mass ratio of cutoff group
    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::CutoffGroupIterator ci;
    CutoffGroup* cg;
    Molecule::AtomIterator ai;
    Atom* atom;
    RealType totalMass;

    /**
     * The mass factor is the relative mass of an atom to the total
     * mass of the cutoff group it belongs to.  By default, all atoms
     * are their own cutoff groups, and therefore have mass factors of
     * 1.  We need some special handling for massless atoms, which
     * will be treated as carrying the entire mass of the cutoff
     * group.
     */
    massFactors_.clear();
    massFactors_.resize(getNAtoms(), 1.0);

    for (mol = beginMolecule(mi); mol != NULL; mol = nextMolecule(mi)) {
      for (cg = mol->beginCutoffGroup(ci); cg != NULL;
           cg = mol->nextCutoffGroup(ci)) {
        totalMass = cg->getMass();
        for (atom = cg->beginAtom(ai); atom != NULL; atom = cg->nextAtom(ai)) {
          // Check for massless groups - set mfact to 1 if true
          if (totalMass != 0)
            massFactors_[atom->getLocalIndex()] = atom->getMass() / totalMass;
          else
            massFactors_[atom->getLocalIndex()] = 1.0;
        }
      }
    }

    // Build the identArray_ and regions_

    identArray_.clear();
    identArray_.reserve(getNAtoms());
    regions_.clear();
    regions_.reserve(getNAtoms());

    for (mol = beginMolecule(mi); mol != NULL; mol = nextMolecule(mi)) {
      int reg = mol->getRegion();
      for (atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
        identArray_.push_back(atom->getIdent());
        regions_.push_back(reg);
      }
    }

    topologyDone_ = true;
  }

  void SimInfo::addProperty(std::shared_ptr<GenericData> genData) {
    properties_.addProperty(genData);
  }

  void SimInfo::removeProperty(const string& propName) {
    properties_.removeProperty(propName);
  }

  std::vector<string> SimInfo::getPropertyNames() {
    return properties_.getPropertyNames();
  }

  std::vector<std::shared_ptr<GenericData>> SimInfo::getProperties() {
    return properties_.getProperties();
  }

  std::shared_ptr<GenericData> SimInfo::getPropertyByName(
      const string& propName) {
    return properties_.getPropertyByName(propName);
  }

  void SimInfo::setSnapshotManager(SnapshotManager* sman) {
    if (sman_ == sman) { return; }
    delete sman_;
    sman_ = sman;

    SimInfo::MoleculeIterator mi;
    Molecule::AtomIterator ai;
    Molecule::RigidBodyIterator rbIter;
    Molecule::CutoffGroupIterator cgIter;
    Molecule::BondIterator bondIter;
    Molecule::BendIterator bendIter;
    Molecule::TorsionIterator torsionIter;
    Molecule::InversionIterator inversionIter;

    Molecule* mol;
    Atom* atom;
    RigidBody* rb;
    CutoffGroup* cg;
    Bond* bond;
    Bend* bend;
    Torsion* torsion;
    Inversion* inversion;

    for (mol = beginMolecule(mi); mol != NULL; mol = nextMolecule(mi)) {
      for (atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
        atom->setSnapshotManager(sman_);
      }
      for (rb = mol->beginRigidBody(rbIter); rb != NULL;
           rb = mol->nextRigidBody(rbIter)) {
        rb->setSnapshotManager(sman_);
      }
      for (cg = mol->beginCutoffGroup(cgIter); cg != NULL;
           cg = mol->nextCutoffGroup(cgIter)) {
        cg->setSnapshotManager(sman_);
      }
      for (bond = mol->beginBond(bondIter); bond != NULL;
           bond = mol->nextBond(bondIter)) {
        bond->setSnapshotManager(sman_);
      }
      for (bend = mol->beginBend(bendIter); bend != NULL;
           bend = mol->nextBend(bendIter)) {
        bend->setSnapshotManager(sman_);
      }
      for (torsion = mol->beginTorsion(torsionIter); torsion != NULL;
           torsion = mol->nextTorsion(torsionIter)) {
        torsion->setSnapshotManager(sman_);
      }
      for (inversion = mol->beginInversion(inversionIter); inversion != NULL;
           inversion = mol->nextInversion(inversionIter)) {
        inversion->setSnapshotManager(sman_);
      }
    }
  }

  ostream& operator<<(ostream& o, SimInfo&) { return o; }

  StuntDouble* SimInfo::getIOIndexToIntegrableObject(int index) {
    if (index >= int(IOIndexToIntegrableObject.size())) {
      snprintf(
          painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
          "SimInfo::getIOIndexToIntegrableObject Error: Integrable Object\n"
          "\tindex exceeds number of known objects!\n");
      painCave.isFatal = 1;
      simError();
      return NULL;
    } else
      return IOIndexToIntegrableObject.at(index);
  }

  void SimInfo::setIOIndexToIntegrableObject(const vector<StuntDouble*>& v) {
    IOIndexToIntegrableObject = v;
  }

  void SimInfo::calcNConstraints() {
#ifdef IS_MPI
    MPI_Allreduce(&nConstraints_, &nGlobalConstraints_, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
#else
    nGlobalConstraints_     = nConstraints_;
#endif
  }
}  // namespace OpenMD
