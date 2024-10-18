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

#include "FluctuatingChargeDamped.hpp"

#include "primitives/Molecule.hpp"
#include "utils/Constants.hpp"
#include "utils/simError.h"

namespace OpenMD {

  FluctuatingChargeDamped::FluctuatingChargeDamped(SimInfo* info) :
      FluctuatingChargePropagator(info), maxIterNum_(4), forceTolerance_(1e-6),
      snap(info->getSnapshotManager()->getCurrentSnapshot()) {}

  void FluctuatingChargeDamped::initialize() {
    FluctuatingChargePropagator::initialize();
    if (hasFlucQ_) {
      if (info_->getSimParams()->haveDt()) {
        dt_  = info_->getSimParams()->getDt();
        dt2_ = dt_ * 0.5;
      } else {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "FluctuatingChargeDamped Error: dt is not set\n");
        painCave.isFatal = 1;
        simError();
      }

      if (!fqParams_->haveDragCoefficient()) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "If you use the FluctuatingChargeDamped\n"
                 "\tpropagator, you must set flucQ dragCoefficient .\n");

        painCave.severity = OPENMD_ERROR;
        painCave.isFatal  = 1;
        simError();
      } else {
        drag_ = fqParams_->getDragCoefficient();
      }
    }
  }

  void FluctuatingChargeDamped::moveA() {
    if (!hasFlucQ_) return;

    SimInfo::MoleculeIterator i;
    Molecule::FluctuatingChargeIterator j;
    Molecule* mol;
    Atom* atom;
    RealType cvel, cpos, cfrc, cmass;

    for (mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i)) {
      for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
           atom = mol->nextFluctuatingCharge(j)) {
        cvel  = atom->getFlucQVel();
        cpos  = atom->getFlucQPos();
        cfrc  = atom->getFlucQFrc();
        cmass = atom->getChargeMass();

        // velocity half step
        cvel += dt2_ * cfrc / cmass;
        // position whole step
        cpos += dt_ * cvel;

        atom->setFlucQVel(cvel);
        atom->setFlucQPos(cpos);
      }
    }
  }

  void FluctuatingChargeDamped::applyConstraints() {
    if (!hasFlucQ_) return;

    SimInfo::MoleculeIterator i;
    Molecule::FluctuatingChargeIterator j;
    Molecule* mol;
    Atom* atom;
    RealType cvel, cfrc, cmass, frictionForce;
    RealType velStep, oldFF;  // used to test for convergence

    for (mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i)) {
      for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
           atom = mol->nextFluctuatingCharge(j)) {
        // What remains contains velocity explicitly, but the velocity
        // required is at the full step: v(t + h), while we have
        // initially the velocity at the half step: v(t + h/2).  We
        // need to iterate to converge the friction force vector.

        // this is the velocity at the half-step:

        cvel = atom->getFlucQVel();

        // estimate velocity at full-step using everything but
        // friction forces:

        cfrc    = atom->getFlucQFrc();
        cmass   = atom->getChargeMass();
        velStep = cvel + dt2_ * cfrc / cmass;

        frictionForce = 0.0;

        // iteration starts here:

        for (int k = 0; k < maxIterNum_; k++) {
          oldFF         = frictionForce;
          frictionForce = -drag_ * velStep;
          // re-estimate velocities at full-step using friction forces:

          velStep = cvel + dt2_ * (cfrc + frictionForce) / cmass;

          // check for convergence

          if (fabs(frictionForce - oldFF) <= forceTolerance_)
            break;  // iteration ends here
        }
        atom->addFlucQFrc(frictionForce);
      }
    }
    fqConstraints_->applyConstraints();
  }

  void FluctuatingChargeDamped::moveB() {
    if (!hasFlucQ_) return;
    SimInfo::MoleculeIterator i;
    Molecule::FluctuatingChargeIterator j;
    Molecule* mol;
    Atom* atom;
    RealType cfrc, cvel, cmass;

    for (mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i)) {
      for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
           atom = mol->nextFluctuatingCharge(j)) {
        cvel  = atom->getFlucQVel();
        cfrc  = atom->getFlucQFrc();
        cmass = atom->getChargeMass();

        // velocity half step
        cvel += (dt2_ * cfrc) / cmass;

        atom->setFlucQVel(cvel);
      }
    }
  }

  void FluctuatingChargeDamped::updateSizes() {}

  RealType FluctuatingChargeDamped::calcConservedQuantity() { return 0.0; }
}  // namespace OpenMD
