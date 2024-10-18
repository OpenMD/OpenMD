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

#include "integrators/NVT.hpp"

#include "primitives/Molecule.hpp"
#include "utils/Constants.hpp"
#include "utils/simError.h"

namespace OpenMD {

  NVT::NVT(SimInfo* info) :
      VelocityVerletIntegrator(info), maxIterNum_(4), chiTolerance_(1e-6) {
    Globals* simParams = info_->getSimParams();

    if (!simParams->getUseIntialExtendedSystemState()) {
      Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();
      snap->setThermostat(make_pair(0.0, 0.0));
    }

    if (!simParams->haveTargetTemp()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "You can't use the NVT integrator without a targetTemp_!\n");
      painCave.isFatal  = 1;
      painCave.severity = OPENMD_ERROR;
      simError();
    } else {
      targetTemp_ = simParams->getTargetTemp();
    }

    // We must set tauThermostat.

    if (!simParams->haveTauThermostat()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "If you use the constant temperature\n"
               "\tintegrator, you must set tauThermostat.\n");

      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    } else {
      tauThermostat_ = simParams->getTauThermostat();
    }

    updateSizes();
  }

  void NVT::doUpdateSizes() {
    oldVel_.resize(info_->getNIntegrableObjects());
    oldJi_.resize(info_->getNIntegrableObjects());
  }

  void NVT::moveA() {
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator j;
    Molecule* mol;
    StuntDouble* sd;
    Vector3d Tb;
    Vector3d ji;
    RealType mass;
    Vector3d vel;
    Vector3d pos;
    Vector3d frc;

    pair<RealType, RealType> thermostat = snap->getThermostat();

    // We need the temperature at time = t for the chi update below:

    RealType instTemp = thermo.getTemperature();

    for (mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i)) {
      for (sd = mol->beginIntegrableObject(j); sd != NULL;
           sd = mol->nextIntegrableObject(j)) {
        vel = sd->getVel();
        pos = sd->getPos();
        frc = sd->getFrc();

        mass = sd->getMass();

        // velocity half step (use chi from previous step here):
        vel += dt2 * Constants::energyConvert / mass * frc -
               dt2 * thermostat.first * vel;

        // position whole step
        pos += dt * vel;

        sd->setVel(vel);
        sd->setPos(pos);

        if (sd->isDirectional()) {
          // convert the torque to body frame
          Tb = sd->lab2Body(sd->getTrq());

          // get the angular momentum, and propagate a half step

          ji = sd->getJ();

          ji +=
              dt2 * Constants::energyConvert * Tb - dt2 * thermostat.first * ji;

          rotAlgo_->rotate(sd, ji, dt);

          sd->setJ(ji);
        }
      }
    }

    flucQ_->moveA();
    rattle_->constraintA();

    // Finally, evolve chi a half step (just like a velocity) using
    // temperature at time t, not time t+dt/2

    thermostat.first += dt2 * (instTemp / targetTemp_ - 1.0) /
                        (tauThermostat_ * tauThermostat_);
    thermostat.second += thermostat.first * dt2;

    snap->setThermostat(thermostat);
  }

  void NVT::moveB() {
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator j;
    Molecule* mol;
    StuntDouble* sd;

    Vector3d Tb;
    Vector3d ji;
    Vector3d vel;
    Vector3d frc;
    RealType mass;
    RealType instTemp;
    int index;
    // Set things up for the iteration:

    pair<RealType, RealType> thermostat = snap->getThermostat();
    RealType oldChi                     = thermostat.first;
    RealType prevChi;

    index = 0;
    for (mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i)) {
      for (sd = mol->beginIntegrableObject(j); sd != NULL;
           sd = mol->nextIntegrableObject(j)) {
        oldVel_[index] = sd->getVel();

        if (sd->isDirectional()) oldJi_[index] = sd->getJ();

        ++index;
      }
    }

    // do the iteration:

    for (int k = 0; k < maxIterNum_; k++) {
      index    = 0;
      instTemp = thermo.getTemperature();

      // evolve chi another half step using the temperature at t + dt/2

      prevChi          = thermostat.first;
      thermostat.first = oldChi + dt2 * (instTemp / targetTemp_ - 1.0) /
                                      (tauThermostat_ * tauThermostat_);

      for (mol = info_->beginMolecule(i); mol != NULL;
           mol = info_->nextMolecule(i)) {
        for (sd = mol->beginIntegrableObject(j); sd != NULL;
             sd = mol->nextIntegrableObject(j)) {
          frc  = sd->getFrc();
          mass = sd->getMass();

          // velocity half step

          vel = oldVel_[index] + dt2 / mass * Constants::energyConvert * frc -
                dt2 * thermostat.first * oldVel_[index];

          sd->setVel(vel);

          if (sd->isDirectional()) {
            // get and convert the torque to body frame

            Tb = sd->lab2Body(sd->getTrq());

            ji = oldJi_[index] + dt2 * Constants::energyConvert * Tb -
                 dt2 * thermostat.first * oldJi_[index];

            sd->setJ(ji);
          }

          ++index;
        }
      }

      rattle_->constraintB();

      if (fabs(prevChi - thermostat.first) <= chiTolerance_) break;
    }

    flucQ_->moveB();

    thermostat.second += dt2 * thermostat.first;
    snap->setThermostat(thermostat);
  }

  void NVT::resetIntegrator() { snap->setThermostat(make_pair(0.0, 0.0)); }

  RealType NVT::calcConservedQuantity() {
    pair<RealType, RealType> thermostat = snap->getThermostat();
    RealType conservedQuantity;
    RealType fkBT;
    RealType Energy;
    RealType thermostat_kinetic;
    RealType thermostat_potential;

    fkBT = info_->getNdf() * Constants::kB * targetTemp_;

    Energy = thermo.getTotalEnergy();

    thermostat_kinetic = fkBT * tauThermostat_ * tauThermostat_ *
                         thermostat.first * thermostat.first /
                         (2.0 * Constants::energyConvert);

    thermostat_potential = fkBT * thermostat.second / Constants::energyConvert;

    conservedQuantity = Energy + thermostat_kinetic + thermostat_potential;

    return conservedQuantity;
  }

}  // namespace OpenMD
