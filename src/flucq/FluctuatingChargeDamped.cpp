/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * This software is provided "AS IS," without a warranty of any
 * kind. All express or implied conditions, representations and
 * warranties, including any implied warranty of merchantability,
 * fitness for a particular purpose or non-infringement, are hereby
 * excluded.  The University of Notre Dame and its licensors shall not
 * be liable for any damages suffered by licensee as a result of
 * using, modifying or distributing the software or its
 * derivatives. In no event will the University of Notre Dame or its
 * licensors be liable for any lost revenue, profit or data, or for
 * direct, indirect, special, consequential, incidental or punitive
 * damages, however caused and regardless of the theory of liability,
 * arising out of the use of or inability to use software, even if the
 * University of Notre Dame has been advised of the possibility of
 * such damages.
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

#include "FluctuatingChargeDamped.hpp"

#include "primitives/Molecule.hpp"
#include "utils/Constants.hpp"
#include "utils/simError.h"

namespace OpenMD {

FluctuatingChargeDamped::FluctuatingChargeDamped(SimInfo* info)
    : FluctuatingChargePropagator(info),
      maxIterNum_(4),
      forceTolerance_(1e-6),
      snap(info->getSnapshotManager()->getCurrentSnapshot()) {}

void FluctuatingChargeDamped::initialize() {
  FluctuatingChargePropagator::initialize();
  if (hasFlucQ_) {
    if (info_->getSimParams()->haveDt()) {
      dt_ = info_->getSimParams()->getDt();
      dt2_ = dt_ * 0.5;
    } else {
      sprintf(painCave.errMsg,
              "FluctuatingChargeDamped Error: dt is not set\n");
      painCave.isFatal = 1;
      simError();
    }

    if (!fqParams_->haveDragCoefficient()) {
      sprintf(painCave.errMsg,
              "If you use the FluctuatingChargeDamped\n"
              "\tpropagator, you must set flucQ dragCoefficient .\n");

      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
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
      cvel = atom->getFlucQVel();
      cpos = atom->getFlucQPos();
      cfrc = atom->getFlucQFrc();
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

      cfrc = atom->getFlucQFrc();
      cmass = atom->getChargeMass();
      velStep = cvel + dt2_ * cfrc / cmass;

      frictionForce = 0.0;

      // iteration starts here:

      for (int k = 0; k < maxIterNum_; k++) {
        oldFF = frictionForce;
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
      cvel = atom->getFlucQVel();
      cfrc = atom->getFlucQFrc();
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
