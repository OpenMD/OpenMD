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

#include "FluctuatingChargeConstraints.hpp"

#ifdef IS_MPI
#include <mpi.h>
#endif

#include "primitives/Molecule.hpp"

namespace OpenMD {

  FluctuatingChargeConstraints::FluctuatingChargeConstraints(SimInfo* info) :
      info_(info), initialized_(false), hasFlucQ_(false),
      constrainRegions_(false) {}

  void FluctuatingChargeConstraints::initialize() {
    if (info_->usesFluctuatingCharges()) {
      if (info_->getNFluctuatingCharges() > 0) { hasFlucQ_ = true; }
    }
    initialized_ = true;
  }

  void FluctuatingChargeConstraints::setConstrainRegions(bool cr) {
    constrainRegions_ = cr;

    if (!initialized_) initialize();

    regionKeys_.clear();
    regionForce_.clear();
    regionCMom_.clear();
    regionCharges_.clear();

    if (constrainRegions_) {
      std::vector<int> localRegions = info_->getRegions();

#ifdef IS_MPI
      int size;
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      int mylen = localRegions.size();

      std::vector<int> counts;
      std::vector<int> displs;

      counts.resize(size, 0);
      displs.resize(size, 0);

      MPI_Allgather(&mylen, 1, MPI_INT, &counts[0], 1, MPI_INT, MPI_COMM_WORLD);

      int total = counts[0];

      for (int i = 1; i < size; i++) {
        total += counts[i];
        displs[i] = displs[i - 1] + counts[i - 1];
      }

      std::vector<int> globalRegions(total, 0);

      MPI_Allgatherv(&localRegions[0], mylen, MPI_INT, &globalRegions[0],
                     &counts[0], &displs[0], MPI_INT, MPI_COMM_WORLD);

      localRegions = globalRegions;
#endif

      std::set<int> regions;
      std::vector<int>::iterator iter;
      for (iter = localRegions.begin(); iter != localRegions.end(); ++iter) {
        if (*iter >= 0) regions.insert(*iter);
      }

      // resize the keys vector to the largest found value for regions.
      regionKeys_.resize(*(regions.end()));
      int which = 0;
      for (std::set<int>::iterator r = regions.begin(); r != regions.end();
           ++r) {
        regionKeys_[(*r)] = which;
        which++;
      }
      regionForce_.resize(regionKeys_.size());
      regionCMom_.resize(regionKeys_.size());
      regionCharges_.resize(regionKeys_.size());
      regionChargeMass_.resize(regionKeys_.size());
    }
  }

  void FluctuatingChargeConstraints::applyConstraints() {
    if (!initialized_) initialize();
    if (!hasFlucQ_) return;

    SimInfo::MoleculeIterator i;
    Molecule::FluctuatingChargeIterator j;
    Molecule* mol;
    Atom* atom;

    RealType frc, systemFrc, molFrc, regionFrc;
    int systemCharges;

    // accumulate the system fluctuating charge forces, but we have
    // separate constraints for any charges in defined regions and for
    // molecules with constrained charges:

    systemFrc     = 0.0;
    systemCharges = 0;
    if (constrainRegions_) {
      std::fill(regionForce_.begin(), regionForce_.end(), 0.0);
      std::fill(regionCharges_.begin(), regionCharges_.end(), 0);
    }

    for (mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i)) {
      if (!mol->constrainTotalCharge()) {
        int region = mol->getRegion();

        for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
             atom = mol->nextFluctuatingCharge(j)) {
          frc = atom->getFlucQFrc();
          if (constrainRegions_ && region >= 0) {
            regionForce_[regionKeys_[region]] += frc;
            regionCharges_[regionKeys_[region]] += 1;
          } else {
            systemFrc += frc;
            systemCharges += 1;
          }
        }
      }
    }

#ifdef IS_MPI
    // in parallel, we need to add up the contributions from all
    // processors:
    MPI_Allreduce(MPI_IN_PLACE, &systemFrc, 1, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &systemCharges, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);

    if (constrainRegions_) {
      MPI_Allreduce(MPI_IN_PLACE, &regionForce_[0], regionForce_.size(),
                    MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &regionCharges_[0], regionCharges_.size(),
                    MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }

#endif

    // divide by the total number of fluctuating charges:
    systemFrc /= systemCharges;

    // do the same in the regions:
    if (constrainRegions_) {
      for (unsigned int i = 0; i < regionForce_.size(); ++i) {
        regionForce_[i] /= regionCharges_[i];
      }
    }

    for (mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i)) {
      molFrc = 0.0;

      if (mol->constrainTotalCharge()) {
        for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
             atom = mol->nextFluctuatingCharge(j)) {
          molFrc += atom->getFlucQFrc();
        }
        molFrc /= mol->getNFluctuatingCharges();

        for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
             atom = mol->nextFluctuatingCharge(j)) {
          frc = atom->getFlucQFrc() - molFrc;
          atom->setFlucQFrc(frc);
        }
      } else {
        int region = mol->getRegion();

        regionFrc = 0.0;
        if (constrainRegions_ && region >= 0) {
          regionFrc = regionForce_[regionKeys_[region]];

          for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
               atom = mol->nextFluctuatingCharge(j)) {
            frc = atom->getFlucQFrc() - regionFrc;
            atom->setFlucQFrc(frc);
          }
        } else {
          for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
               atom = mol->nextFluctuatingCharge(j)) {
            frc = atom->getFlucQFrc() - systemFrc;
            atom->setFlucQFrc(frc);
          }
        }
      }
    }
  }

  void FluctuatingChargeConstraints::applyConstraintsOnChargeVelocities() {
    if (!initialized_) initialize();
    if (!hasFlucQ_) return;

    SimInfo::MoleculeIterator i;
    Molecule::FluctuatingChargeIterator j;
    Molecule* mol;
    Atom* atom;

    RealType flucqP, systemCMom, regionCMom, molCMom, molFlucQMass, flucqW;
    RealType systemChargeMass;

    // accumulate the system fluctuating charge velocities, but we have
    // separate constraints for any charges in defined regions and for
    // molecules with constrained charges:

    systemCMom       = 0.0;
    systemChargeMass = 0.0;
    if (constrainRegions_) {
      std::fill(regionCMom_.begin(), regionCMom_.end(), 0.0);
      std::fill(regionChargeMass_.begin(), regionChargeMass_.end(), 0);
    }

    for (mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i)) {
      if (!mol->constrainTotalCharge()) {
        int region = mol->getRegion();

        for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
             atom = mol->nextFluctuatingCharge(j)) {
          flucqP = atom->getFlucQVel() * atom->getChargeMass();
          if (constrainRegions_ && region >= 0) {
            regionCMom_[regionKeys_[region]] += flucqP;
            regionChargeMass_[regionKeys_[region]] += atom->getChargeMass();
          } else {
            systemCMom += flucqP;
            systemChargeMass += atom->getChargeMass();
          }
        }
      }
    }

#ifdef IS_MPI
    // in parallel, we need to add up the contributions from all
    // processors:
    MPI_Allreduce(MPI_IN_PLACE, &systemCMom, 1, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &systemChargeMass, 1, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);

    if (constrainRegions_) {
      MPI_Allreduce(MPI_IN_PLACE, &regionCMom_[0], regionCMom_.size(),
                    MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &regionChargeMass_[0],
                    regionChargeMass_.size(), MPI_REALTYPE, MPI_SUM,
                    MPI_COMM_WORLD);
    }
#endif

    // divide by the total number of fluctuating charges:
    systemCMom /= systemChargeMass;

    // do the same in the regions:
    if (constrainRegions_) {
      for (unsigned int i = 0; i < regionCMom_.size(); ++i) {
        regionCMom_[i] /= regionChargeMass_[i];
      }
    }

    for (mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i)) {
      molCMom      = 0.0;
      molFlucQMass = 0.0;

      if (mol->constrainTotalCharge()) {
        for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
             atom = mol->nextFluctuatingCharge(j)) {
          molCMom += atom->getFlucQVel() * atom->getChargeMass();
          molFlucQMass += atom->getChargeMass();
        }
        molCMom /= molFlucQMass;

        for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
             atom = mol->nextFluctuatingCharge(j)) {
          flucqW = atom->getFlucQVel() - molCMom;
          atom->setFlucQVel(flucqW);
        }
      } else {
        int region = mol->getRegion();

        regionCMom = 0.0;
        if (constrainRegions_ && region >= 0) {
          regionCMom = regionCMom_[regionKeys_[region]];

          for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
               atom = mol->nextFluctuatingCharge(j)) {
            flucqW = atom->getFlucQVel() - regionCMom;
            atom->setFlucQVel(flucqW);
          }
        } else {
          for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
               atom = mol->nextFluctuatingCharge(j)) {
            flucqW = atom->getFlucQVel() - systemCMom;
            atom->setFlucQVel(flucqW);
          }
        }
      }
    }
  }

  int FluctuatingChargeConstraints::getNumberOfFlucQConstraints() {
    int nConstraints = 0;
    if (!initialized_) initialize();
    if (!hasFlucQ_) return 0;
    SimInfo::MoleculeIterator i;
    Molecule* mol;
    int systemConstrain = 0;
    for (mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i)) {
      if (mol->constrainTotalCharge()) {
        nConstraints++;
      } else {
        int region = mol->getRegion();
        if (!constrainRegions_ || region < 0) { systemConstrain = 1; }
      }
    }
    return nConstraints + regionCMom_.size() + systemConstrain;
  }

  int FluctuatingChargeConstraints::getNumberOfFlucQAtoms() {
    int nFlucq = 0;
    if (!initialized_) initialize();
    if (!hasFlucQ_) return 0;
    SimInfo::MoleculeIterator i;
    Molecule::FluctuatingChargeIterator j;
    Molecule* mol;
    Atom* atom;
    for (mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i)) {
      if (mol->constrainTotalCharge()) {
        nFlucq += mol->getNFluctuatingCharges();
      } else {
        int region = mol->getRegion();
        if (constrainRegions_ && region >= 0) {
          for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
               atom = mol->nextFluctuatingCharge(j)) {
            nFlucq++;
          }
        } else {
          for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
               atom = mol->nextFluctuatingCharge(j)) {
            nFlucq++;
          }
        }
      }
    }
    return nFlucq;
  }
}  // namespace OpenMD
