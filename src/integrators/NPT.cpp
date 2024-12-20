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

#include "integrators/NPT.hpp"

#include <cmath>

#include "brains/SimInfo.hpp"
#include "brains/Thermo.hpp"
#include "math/SquareMatrix3.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Constants.hpp"
#include "utils/simError.h"

// Basic isotropic thermostating and barostating via the Melchionna
// modification of the Hoover algorithm:
//
//    Melchionna, S., Ciccotti, G., and Holian, B. L., 1993,
//       Molec. Phys., 78, 533.
//
//           and
//
//    Hoover, W. G., 1986, Phys. Rev. A, 34, 2499.

namespace OpenMD {

  NPT::NPT(SimInfo* info) :
      VelocityVerletIntegrator(info), etaTolerance(1e-6), chiTolerance(1e-6),
      maxIterNum_(4) {
    Globals* simParams = info_->getSimParams();

    if (!simParams->getUseIntialExtendedSystemState()) {
      Snapshot* currSnapshot =
          info_->getSnapshotManager()->getCurrentSnapshot();
      currSnapshot->setThermostat(make_pair(0.0, 0.0));
      currSnapshot->setBarostat(Mat3x3d(0.0));
    }

    if (!simParams->haveTargetTemp()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "You can't use the NVT integrator without a targetTemp!\n");
      painCave.isFatal  = 1;
      painCave.severity = OPENMD_ERROR;
      simError();
    } else {
      targetTemp = simParams->getTargetTemp();
    }

    // We must set tauThermostat
    if (!simParams->haveTauThermostat()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "If you use the constant temperature\n"
               "\tintegrator, you must set tauThermostat.\n");

      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    } else {
      tauThermostat = simParams->getTauThermostat();
    }

    if (!simParams->haveTargetPressure()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "NPT error: You can't use the NPT integrator\n"
               "   without a targetPressure!\n");

      painCave.isFatal = 1;
      simError();
    } else {
      targetPressure = simParams->getTargetPressure();
    }

    if (!simParams->haveTauBarostat()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "If you use the NPT integrator, you must set tauBarostat.\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    } else {
      tauBarostat = simParams->getTauBarostat();
    }

    tt2 = tauThermostat * tauThermostat;
    tb2 = tauBarostat * tauBarostat;

    updateSizes();
  }

  void NPT::doUpdateSizes() {
    oldPos.resize(info_->getNIntegrableObjects());
    oldVel.resize(info_->getNIntegrableObjects());
    oldJi.resize(info_->getNIntegrableObjects());
  }

  void NPT::moveA() {
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator j;
    Molecule* mol;
    StuntDouble* sd;
    Vector3d Tb, ji;
    RealType mass;
    Vector3d vel;
    Vector3d pos;
    Vector3d frc;
    Vector3d sc;
    int index;

    thermostat = snap->getThermostat();
    loadEta();

    instaTemp  = thermo.getTemperature();
    press      = thermo.getPressureTensor();
    instaPress = Constants::pressureConvert *
                 (press(0, 0) + press(1, 1) + press(2, 2)) / 3.0;
    instaVol = thermo.getVolume();

    Vector3d COM = thermo.getCom();

    // evolve velocity half step

    calcVelScale();

    for (mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i)) {
      for (sd = mol->beginIntegrableObject(j); sd != NULL;
           sd = mol->nextIntegrableObject(j)) {
        vel = sd->getVel();
        frc = sd->getFrc();

        mass = sd->getMass();

        getVelScaleA(sc, vel);

        // velocity half step  (use chi from previous step here):

        vel += dt2 * Constants::energyConvert / mass * frc - dt2 * sc;
        sd->setVel(vel);

        if (sd->isDirectional()) {
          // get and convert the torque to body frame

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
    // evolve chi and eta  half step

    thermostat.first += dt2 * (instaTemp / targetTemp - 1.0) / tt2;

    evolveEtaA();

    // calculate the integral of chidt
    thermostat.second += dt2 * thermostat.first;

    flucQ_->moveA();

    index = 0;
    for (mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i)) {
      for (sd = mol->beginIntegrableObject(j); sd != NULL;
           sd = mol->nextIntegrableObject(j)) {
        oldPos[index++] = sd->getPos();
      }
    }

    // the first estimation of r(t+dt) is equal to  r(t)

    for (int k = 0; k < maxIterNum_; k++) {
      index = 0;
      for (mol = info_->beginMolecule(i); mol != NULL;
           mol = info_->nextMolecule(i)) {
        for (sd = mol->beginIntegrableObject(j); sd != NULL;
             sd = mol->nextIntegrableObject(j)) {
          vel = sd->getVel();
          pos = sd->getPos();

          this->getPosScale(pos, COM, index, sc);

          pos = oldPos[index] + dt * (vel + sc);
          sd->setPos(pos);

          ++index;
        }
      }

      rattle_->constraintA();
    }

    // Scale the box after all the positions have been moved:

    this->scaleSimBox();

    snap->setThermostat(thermostat);

    saveEta();
  }

  void NPT::moveB(void) {
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator j;
    Molecule* mol;
    StuntDouble* sd;
    int index;
    Vector3d Tb;
    Vector3d ji;
    Vector3d sc;
    Vector3d vel;
    Vector3d frc;
    RealType mass;

    thermostat      = snap->getThermostat();
    RealType oldChi = thermostat.first;
    RealType prevChi;

    loadEta();

    // save velocity and angular momentum
    index = 0;
    for (mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i)) {
      for (sd = mol->beginIntegrableObject(j); sd != NULL;
           sd = mol->nextIntegrableObject(j)) {
        oldVel[index] = sd->getVel();

        if (sd->isDirectional()) oldJi[index] = sd->getJ();

        ++index;
      }
    }

    // do the iteration:
    instaVol = thermo.getVolume();

    for (int k = 0; k < maxIterNum_; k++) {
      instaTemp  = thermo.getTemperature();
      instaPress = thermo.getPressure();

      // evolve chi another half step using the temperature at t + dt/2
      prevChi          = thermostat.first;
      thermostat.first = oldChi + dt2 * (instaTemp / targetTemp - 1.0) / tt2;

      // evolve eta
      this->evolveEtaB();
      this->calcVelScale();

      index = 0;
      for (mol = info_->beginMolecule(i); mol != NULL;
           mol = info_->nextMolecule(i)) {
        for (sd = mol->beginIntegrableObject(j); sd != NULL;
             sd = mol->nextIntegrableObject(j)) {
          frc  = sd->getFrc();
          mass = sd->getMass();

          getVelScaleB(sc, index);

          // velocity half step
          vel = oldVel[index] + dt2 * Constants::energyConvert / mass * frc -
                dt2 * sc;

          sd->setVel(vel);

          if (sd->isDirectional()) {
            // get and convert the torque to body frame
            Tb = sd->lab2Body(sd->getTrq());

            ji = oldJi[index] + dt2 * Constants::energyConvert * Tb -
                 dt2 * thermostat.first * oldJi[index];

            sd->setJ(ji);
          }

          ++index;
        }
      }

      rattle_->constraintB();

      if ((fabs(prevChi - thermostat.first) <= chiTolerance) &&
          this->etaConverged())
        break;
    }

    // calculate integral of chidt
    thermostat.second += dt2 * thermostat.first;

    snap->setThermostat(thermostat);

    flucQ_->moveB();
    saveEta();
  }

  void NPT::resetIntegrator() {
    snap->setThermostat(make_pair(0.0, 0.0));
    resetEta();
  }

  void NPT::resetEta() {
    Mat3x3d etaMat(0.0);
    snap->setBarostat(etaMat);
  }
}  // namespace OpenMD
