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

#include "FluctuatingChargeNVT.hpp"

#include "primitives/Molecule.hpp"
#include "utils/Constants.hpp"
#include "utils/simError.h"

namespace OpenMD {

  FluctuatingChargeNVT::FluctuatingChargeNVT(SimInfo* info) :
      FluctuatingChargePropagator(info), maxIterNum_(4), chiTolerance_(1e-6),
      snap(info->getSnapshotManager()->getCurrentSnapshot()), thermo(info) {}

  void FluctuatingChargeNVT::initialize() {
    FluctuatingChargePropagator::initialize();
    if (hasFlucQ_) {
      if (info_->getSimParams()->haveDt()) {
        dt_  = info_->getSimParams()->getDt();
        dt2_ = dt_ * 0.5;
      } else {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "FluctuatingChargeNVT Error: dt is not set\n");
        painCave.isFatal = 1;
        simError();
      }

      if (!info_->getSimParams()->getUseIntialExtendedSystemState()) {
        snap->setElectronicThermostat(make_pair(0.0, 0.0));
      }

      if (!fqParams_->haveTargetTemp()) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "You can't use the FluctuatingChargeNVT "
                 "propagator without a flucQ.targetTemp!\n");
        painCave.isFatal  = 1;
        painCave.severity = OPENMD_ERROR;
        simError();
      } else {
        targetTemp_ = fqParams_->getTargetTemp();
      }

      // We must set tauThermostat.

      if (!fqParams_->haveTauThermostat()) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "If you use the FluctuatingChargeNVT\n"
                 "\tpropagator, you must set flucQ.tauThermostat .\n");

        painCave.severity = OPENMD_ERROR;
        painCave.isFatal  = 1;
        simError();
      } else {
        tauThermostat_ = fqParams_->getTauThermostat();
      }
      updateSizes();
    }
  }

  void FluctuatingChargeNVT::moveA() {
    if (!hasFlucQ_) return;

    SimInfo::MoleculeIterator i;
    Molecule::FluctuatingChargeIterator j;
    Molecule* mol;
    Atom* atom;
    RealType cvel, cpos, cfrc, cmass;

    pair<RealType, RealType> thermostat = snap->getElectronicThermostat();
    RealType chi                        = thermostat.first;
    RealType integralOfChidt            = thermostat.second;
    RealType instTemp                   = thermo.getElectronicTemperature();

    for (mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i)) {
      for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
           atom = mol->nextFluctuatingCharge(j)) {
        cvel  = atom->getFlucQVel();
        cpos  = atom->getFlucQPos();
        cfrc  = atom->getFlucQFrc();
        cmass = atom->getChargeMass();

        // velocity half step
        cvel += dt2_ * cfrc / cmass - dt2_ * chi * cvel;
        // position whole step
        cpos += dt_ * cvel;

        atom->setFlucQVel(cvel);
        atom->setFlucQPos(cpos);
      }
    }

    chi += dt2_ * (instTemp / targetTemp_ - 1.0) /
           (tauThermostat_ * tauThermostat_);

    integralOfChidt += chi * dt2_;
    snap->setElectronicThermostat(make_pair(chi, integralOfChidt));
  }

  void FluctuatingChargeNVT::updateSizes() {
    oldVel_.resize(info_->getNFluctuatingCharges());
  }

  void FluctuatingChargeNVT::moveB() {
    if (!hasFlucQ_) return;
    SimInfo::MoleculeIterator i;
    Molecule::FluctuatingChargeIterator j;
    Molecule* mol;
    Atom* atom;
    RealType instTemp;
    pair<RealType, RealType> thermostat = snap->getElectronicThermostat();
    RealType chi                        = thermostat.first;
    RealType oldChi                     = chi;
    RealType prevChi;
    RealType integralOfChidt = thermostat.second;
    int index;
    RealType cfrc, cvel, cmass;

    index = 0;
    for (mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i)) {
      for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
           atom = mol->nextFluctuatingCharge(j)) {
        oldVel_[index] = atom->getFlucQVel();
        ++index;
      }
    }

    // do the iteration:

    for (int k = 0; k < maxIterNum_; k++) {
      index    = 0;
      instTemp = thermo.getElectronicTemperature();
      // evolve chi another half step using the temperature at t + dt/2
      prevChi = chi;
      chi     = oldChi + dt2_ * (instTemp / targetTemp_ - 1.0) /
                         (tauThermostat_ * tauThermostat_);

      for (mol = info_->beginMolecule(i); mol != NULL;
           mol = info_->nextMolecule(i)) {
        for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
             atom = mol->nextFluctuatingCharge(j)) {
          cfrc  = atom->getFlucQFrc();
          cmass = atom->getChargeMass();

          // velocity half step
          cvel = oldVel_[index] + dt2_ * cfrc / cmass -
                 dt2_ * chi * oldVel_[index];
          atom->setFlucQVel(cvel);
          ++index;
        }
      }
      if (fabs(prevChi - chi) <= chiTolerance_) break;
    }
    integralOfChidt += dt2_ * chi;
    snap->setElectronicThermostat(make_pair(chi, integralOfChidt));
  }

  void FluctuatingChargeNVT::resetPropagator() {
    if (!hasFlucQ_) return;
    snap->setElectronicThermostat(make_pair(0.0, 0.0));
  }

  RealType FluctuatingChargeNVT::calcConservedQuantity() {
    if (!hasFlucQ_) return 0.0;
    pair<RealType, RealType> thermostat = snap->getElectronicThermostat();
    RealType chi                        = thermostat.first;
    RealType integralOfChidt            = thermostat.second;
    RealType fkBT =
        info_->getNFluctuatingCharges() * Constants::kB * targetTemp_;

    RealType thermostat_kinetic = fkBT * tauThermostat_ * tauThermostat_ * chi *
                                  chi / (2.0 * Constants::energyConvert);

    RealType thermostat_potential =
        fkBT * integralOfChidt / Constants::energyConvert;

    return thermostat_kinetic + thermostat_potential;
  }
}  // namespace OpenMD
