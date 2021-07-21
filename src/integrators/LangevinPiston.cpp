/*
 * Copyright (c) 2004-2021 The University of Notre Dame. All Rights Reserved.
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
 * warranties, including any implied warranty of merhantability,
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

#include "integrators/LangevinPiston.hpp"

#include <cmath>
#include <random>

#ifdef IS_MPI
#include <mpi.h>
#endif

#include "brains/SimInfo.hpp"
#include "brains/Thermo.hpp"
#include "integrators/IntegratorCreator.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Constants.hpp"
#include "utils/simError.h"

namespace OpenMD {

  LangevinPiston::LangevinPiston(SimInfo* info) : NPT(info) {
    // NkBT has units of amu Ang^2 fs^-2 :
    NkBT = info_->getNGlobalIntegrableObjects() * Constants::kB * targetTemp;

    // W_ has units of amu Ang^2
    // W_ = 3.0 * NkBT * tb2;
    W_ = NkBT * tb2;  // our eta scales all three box directions

    // gamma_ has units of fs^-1
    if (!simParams->haveLangevinPistonDrag()) {
      sprintf(painCave.errMsg, "To use the LangevinPiston integrator, you must "
                               "set langevinPistonDrag "
                               "(fs^-1).\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    } else {
      gamma_ = simParams->getLangevinPistonDrag();
    }

#ifdef IS_MPI
    if (worldRank == 0) {
#endif
      randNumGen_ = info->getRandomNumberGenerator();

      // standard deviation units: amu Angs^2 fs^-2
      RealType stdDev =
          std::sqrt(2.0 * W_ * gamma_ * Constants::kB * targetTemp / dt);

      forceDistribution_ = std::normal_distribution<RealType>(0.0, stdDev);
#ifdef IS_MPI
    }
#endif

    // randomForce will have units amu Ang^2 fs^-2:
    genRandomForce(randomForce_);
  }

  void LangevinPiston::moveA() {
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

        // velocity half step

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

  void LangevinPiston::moveB(void) {
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

  void LangevinPiston::evolveEtaA() {
    // volume is Angs^3
    // pressures are in atm
    // pressureConvert takes amu*fs^-2*Ang^-1 -> atm
    eta += dt2 * (instaVol * (instaPress - targetPressure) /
                      (Constants::pressureConvert * W_) -
                  gamma_ * eta + randomForce_ / W_);
    oldEta = eta;
  }

  void LangevinPiston::evolveEtaB() {
    prevEta = eta;

    genRandomForce(randomForce_);

    eta = oldEta + dt2 * (instaVol * (instaPress - targetPressure) /
                              (Constants::pressureConvert * W_) -
                          gamma_ * eta + randomForce_ / W_);
  }

  void LangevinPiston::calcVelScale() { vScale = eta; }

  void LangevinPiston::getVelScaleA(Vector3d& sc, const Vector3d& vel) {
    sc = vScale * vel;
  }

  void LangevinPiston::getVelScaleB(Vector3d& sc, int index) {
    sc = vScale * oldVel[index];
  }

  void LangevinPiston::getPosScale(const Vector3d& pos, const Vector3d& COM,
                                   int index, Vector3d& sc) {
    Vector3d rj = (oldPos[index] + pos) / (RealType)2.0 - COM;
    sc          = eta * rj;
  }

  void LangevinPiston::scaleSimBox() {
    RealType scaleFactor;

    // This is from solving the first order equation that defines eta
    scaleFactor = exp(dt * eta);

    if ((scaleFactor > 1.1) || (scaleFactor < 0.9)) {
      sprintf(painCave.errMsg,
              "LangevinPiston error: Attempting a Box scaling of more than 10 "
              "percent\n"
              " check your tauBarostat, as it is probably too small!\n"
              " eta = %lf, scaleFactor = %lf\n",
              eta, scaleFactor);
      painCave.isFatal = 1;
      simError();
    } else {
      Mat3x3d hmat = snap->getHmat();
      hmat *= scaleFactor;
      snap->setHmat(hmat);
    }
  }

  bool LangevinPiston::etaConverged() {
    return (fabs(prevEta - eta) <= etaTolerance);
  }

  void LangevinPiston::loadEta() {
    Mat3x3d etaMat = snap->getBarostat();
    eta            = etaMat(0, 0);
  }

  void LangevinPiston::saveEta() {
    Mat3x3d etaMat(0.0);
    etaMat(0, 0) = eta;
    etaMat(1, 1) = eta;
    etaMat(2, 2) = eta;
    snap->setBarostat(etaMat);
  }

  void LangevinPiston::genRandomForce(RealType& randomForce) {
#ifdef IS_MPI
    if (worldRank == 0) {
#endif
      randomForce = forceDistribution_(*randNumGen_);
#ifdef IS_MPI
    }
    // push this out to the other processors
    // Same command on all nodes:
    MPI_Bcast(&randomForce, 1, MPI_REALTYPE, 0, MPI_COMM_WORLD);
#endif

    return;
  }
}  // namespace OpenMD
