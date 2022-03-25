/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
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

#include "FluctuatingChargeLangevin.hpp"

#include <cmath>

#include "primitives/Molecule.hpp"
#include "utils/Constants.hpp"
#include "utils/simError.h"

namespace OpenMD {

  FluctuatingChargeLangevin::FluctuatingChargeLangevin(SimInfo* info) :
      FluctuatingChargePropagator(info), maxIterNum_(4), forceTolerance_(1e-6),
      snap_(info->getSnapshotManager()->getCurrentSnapshot()),
      randNumGen_(info->getRandomNumberGenerator()) {}

  void FluctuatingChargeLangevin::initialize() {
    FluctuatingChargePropagator::initialize();
    if (hasFlucQ_) {
      if (info_->getSimParams()->haveDt()) {
        dt_  = info_->getSimParams()->getDt();
        dt2_ = dt_ * 0.5;
      } else {
        sprintf(painCave.errMsg,
                "FluctuatingChargeLangevin Error: dt is not set\n");
        painCave.isFatal = 1;
        simError();
      }

      if (!fqParams_->haveTargetTemp()) {
        sprintf(painCave.errMsg, "You can't use the FluctuatingChargeLangevin "
                                 "propagator without a flucQ.targetTemp!\n");
        painCave.isFatal  = 1;
        painCave.severity = OPENMD_ERROR;
        simError();
      } else {
        targetTemp_ = fqParams_->getTargetTemp();
      }

      // We must set tauThermostat.

      if (!fqParams_->haveDragCoefficient()) {
        sprintf(painCave.errMsg,
                "If you use the FluctuatingChargeLangevin\n"
                "\tpropagator, you must set flucQ.dragCoefficient .\n");

        painCave.severity = OPENMD_ERROR;
        painCave.isFatal  = 1;
        simError();
      } else {
        drag_ = fqParams_->getDragCoefficient();
      }
    }

    RealType stdDev =
        std::sqrt(2.0 * Constants::kb * targetTemp_ * drag_ / dt_);

    forceDistribution_ = std::normal_distribution<RealType>(0.0, stdDev);
  }

  void FluctuatingChargeLangevin::moveA() {
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

  void FluctuatingChargeLangevin::applyConstraints() {
    if (!hasFlucQ_) return;

    SimInfo::MoleculeIterator i;
    Molecule::FluctuatingChargeIterator j;
    Molecule* mol;
    Atom* atom;
    RealType cvel, cfrc, cmass, randomForce, frictionForce;
    RealType velStep, oldFF;  // used to test for convergence

    for (mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i)) {
      for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
           atom = mol->nextFluctuatingCharge(j)) {
        randomForce = forceDistribution_(*randNumGen_);
        atom->addFlucQFrc(randomForce);

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

  void FluctuatingChargeLangevin::moveB() {
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
}  // namespace OpenMD
