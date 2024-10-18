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

#include "integrators/NPA.hpp"

#include "brains/SimInfo.hpp"
#include "brains/Thermo.hpp"
#include "integrators/IntegratorCreator.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Constants.hpp"
#include "utils/simError.h"

namespace OpenMD {

  void NPA::moveA() {
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

          ji += dt2 * Constants::energyConvert * Tb;

          rotAlgo_->rotate(sd, ji, dt);

          sd->setJ(ji);
        }
      }
    }
    // evolve eta a half step

    evolveEtaA();
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

    saveEta();
  }

  void NPA::moveB(void) {
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

    instaVol = thermo.getVolume();
    for (int k = 0; k < maxIterNum_; k++) {
      instaTemp  = thermo.getTemperature();
      instaPress = thermo.getPressure();

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

            ji = oldJi[index] + dt2 * Constants::energyConvert * Tb;

            sd->setJ(ji);
          }

          ++index;
        }
      }

      rattle_->constraintB();

      if (this->etaConverged()) break;
    }

    flucQ_->moveB();
    saveEta();
  }

  void NPA::evolveEtaA() {
    eta(2, 2) += dt2 * instaVol *
                 (press(2, 2) - targetPressure / Constants::pressureConvert) /
                 (NkBT * tb2);
    oldEta = eta;
  }

  void NPA::evolveEtaB() {
    prevEta = eta;
    eta(2, 2) =
        oldEta(2, 2) +
        dt2 * instaVol *
            (press(2, 2) - targetPressure / Constants::pressureConvert) /
            (NkBT * tb2);
  }

  void NPA::calcVelScale() {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        vScale(i, j) = eta(i, j);
      }
    }
  }

  void NPA::getVelScaleA(Vector3d& sc, const Vector3d& vel) {
    sc = vScale * vel;
  }

  void NPA::getVelScaleB(Vector3d& sc, int index) {
    sc = vScale * oldVel[index];
  }

  void NPA::getPosScale(const Vector3d& pos, const Vector3d& COM, int index,
                        Vector3d& sc) {
    Vector3d rj = (oldPos[index] + pos) / (RealType)2.0 - COM;
    sc          = eta * rj;
  }

  void NPA::scaleSimBox() {
    Mat3x3d scaleMat;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        scaleMat(i, j) = 0.0;
        if (i == j) { scaleMat(i, j) = 1.0; }
      }
    }

    scaleMat(2, 2) = exp(dt * eta(2, 2));
    Mat3x3d hmat   = snap->getHmat();
    hmat           = hmat * scaleMat;
    snap->setHmat(hmat);
  }

  bool NPA::etaConverged() {
    int i;
    RealType diffEta, sumEta;

    sumEta = 0;
    for (i = 0; i < 3; i++) {
      sumEta += pow(prevEta(i, i) - eta(i, i), 2);
    }

    diffEta = sqrt(sumEta / 3.0);

    return (diffEta <= etaTolerance);
  }

  RealType NPA::calcConservedQuantity() {
    loadEta();

    // We need NkBT a lot, so just set it here: This is the RAW number
    // of integrableObjects, so no subtraction or addition of constraints or
    // orientational degrees of freedom:
    NkBT = info_->getNGlobalIntegrableObjects() * Constants::kB * targetTemp;

    RealType conservedQuantity;
    RealType totalEnergy;
    RealType barostat_kinetic;
    RealType barostat_potential;
    RealType trEta;

    totalEnergy = thermo.getTotalEnergy();

    SquareMatrix<RealType, 3> tmp = eta.transpose() * eta;
    trEta                         = tmp.trace();

    barostat_kinetic = NkBT * tb2 * trEta / (2.0 * Constants::energyConvert);

    barostat_potential =
        (targetPressure * thermo.getVolume() / Constants::pressureConvert) /
        Constants::energyConvert;

    conservedQuantity = totalEnergy + barostat_kinetic + barostat_potential;

    return conservedQuantity;
  }

  void NPA::loadEta() { eta = snap->getBarostat(); }

  void NPA::saveEta() { snap->setBarostat(eta); }
}  // namespace OpenMD
