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

#include "selection/DistanceFinder.hpp"

#ifdef IS_MPI
#include <mpi.h>
#endif

#include "primitives/Molecule.hpp"

namespace OpenMD {

  DistanceFinder::DistanceFinder(SimInfo* info) : info_(info) {
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

  SelectionSet DistanceFinder::find(const SelectionSet& bs,
                                    RealType distance) {
    return findImpl(bs, distance, -1);
  }

  SelectionSet DistanceFinder::find(const SelectionSet& bs, RealType distance,
                                    int frame) {
    return findImpl(bs, distance, frame);
  }

  SelectionSet DistanceFinder::findImpl(const SelectionSet& bs,
                                        RealType distance, int frame) {
    const bool useCurrent = (frame < 0);
    Snapshot* currSnapshot = useCurrent
      ? info_->getSnapshotManager()->getCurrentSnapshot()
      : info_->getSnapshotManager()->getSnapshot(frame);

    // Position-retrieval helpers that dispatch on frame
    auto getPos = [useCurrent, frame](StuntDouble* sd) -> Vector3d {
      return useCurrent ? sd->getPos() : sd->getPos(frame);
    };
    auto getCom = [useCurrent, frame](Molecule* mol) -> Vector3d {
      return useCurrent ? mol->getCom() : mol->getCom(frame);
    };

    SelectionSet bsResult(nObjects_);
    assert(bsResult.size() == bs.size());

#ifdef IS_MPI
    int mol;
    int proc;
    RealType data[3];
    int worldRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
#endif

    for (unsigned int j = 0; j < stuntdoubles_.size(); ++j) {
      if (stuntdoubles_[j] != NULL) {
        if (stuntdoubles_[j]->isRigidBody()) {
          RigidBody* rb = static_cast<RigidBody*>(stuntdoubles_[j]);
          useCurrent ? rb->updateAtoms() : rb->updateAtoms(frame);
        }
      }
    }

    // Helper: check if pos is within distance of centerPos (min-image)
    // and set the bit if so.
    auto checkDistance = [&](const Vector3d& centerPos, const Vector3d& pos,
			     OpenMDBitSet& bits, size_t idx) {
      Vector3d r = centerPos - pos;
      currSnapshot->wrapVector(r);
      if (r.length() <= distance) { bits.setBitOn(idx); }
    };

    SelectionSet bsTemp = bs;
    bsTemp              = bsTemp.parallelReduce();

    for (size_t i = 0; i < bsTemp.bitsets_[STUNTDOUBLE].size(); ++i) {
      if (!bsTemp.bitsets_[STUNTDOUBLE][i]) continue;

      Vector3d centerPos;

#ifdef IS_MPI
      mol  = info_->getGlobalMolMembership(i);
      proc = info_->getMolToProc(mol);

      if (proc == worldRank) {
        centerPos = getPos(stuntdoubles_[i]);
        data[0]   = centerPos.x();
        data[1]   = centerPos.y();
        data[2]   = centerPos.z();
        MPI_Bcast(data, 3, MPI_REALTYPE, proc, MPI_COMM_WORLD);
      } else {
        MPI_Bcast(data, 3, MPI_REALTYPE, proc, MPI_COMM_WORLD);
        centerPos = Vector3d(data);
      }
#else
      centerPos = getPos(stuntdoubles_[i]);
#endif

      // Molecules
      for (size_t j = 0; j < molecules_.size(); ++j) {
        if (molecules_[j] != NULL)
          checkDistance(centerPos, getCom(molecules_[j]),
			bsResult.bitsets_[MOLECULE], j);
      }

      // StuntDoubles (atoms and rigid bodies)
      for (size_t j = 0; j < stuntdoubles_.size(); ++j) {
        if (stuntdoubles_[j] != NULL)
          checkDistance(centerPos, getPos(stuntdoubles_[j]),
			bsResult.bitsets_[STUNTDOUBLE], j);
      }

      // Bonds (midpoint of two atoms)
      for (size_t j = 0; j < bonds_.size(); ++j) {
        if (bonds_[j] != NULL) {
          Vector3d loc = getPos(bonds_[j]->getAtomA())
	    + getPos(bonds_[j]->getAtomB());
          checkDistance(centerPos, loc / 2.0,
			bsResult.bitsets_[BOND], j);
        }
      }

      // Bends (centroid of three atoms)
      for (size_t j = 0; j < bends_.size(); ++j) {
        if (bends_[j] != NULL) {
          Vector3d loc = getPos(bends_[j]->getAtomA())
	    + getPos(bends_[j]->getAtomB())
	    + getPos(bends_[j]->getAtomC());
          checkDistance(centerPos, loc / 3.0,
			bsResult.bitsets_[BEND], j);
        }
      }

      // Torsions (centroid of four atoms)
      for (size_t j = 0; j < torsions_.size(); ++j) {
        if (torsions_[j] != NULL) {
          Vector3d loc = getPos(torsions_[j]->getAtomA())
	    + getPos(torsions_[j]->getAtomB())
	    + getPos(torsions_[j]->getAtomC())
	    + getPos(torsions_[j]->getAtomD());
          checkDistance(centerPos, loc / 4.0,
			bsResult.bitsets_[TORSION], j);
        }
      }

      // Inversions (centroid of four atoms)
      for (size_t j = 0; j < inversions_.size(); ++j) {
        if (inversions_[j] != NULL) {
          Vector3d loc = getPos(inversions_[j]->getAtomA())
	    + getPos(inversions_[j]->getAtomB())
	    + getPos(inversions_[j]->getAtomC())
	    + getPos(inversions_[j]->getAtomD());
          checkDistance(centerPos, loc / 4.0,
			bsResult.bitsets_[INVERSION], j);
        }
      }
    }
    return bsResult;
  }
}  // namespace OpenMD
