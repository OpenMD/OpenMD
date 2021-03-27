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

#include "integrators/NPAT.hpp"

#include "brains/SimInfo.hpp"
#include "brains/Thermo.hpp"
#include "integrators/IntegratorCreator.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Constants.hpp"
#include "utils/simError.h"

namespace OpenMD {

NPAT::NPAT(SimInfo* info) : NPT(info) {
  Globals* simParams = info_->getSimParams();

  // Default value of privilegedAxis is "z"
  if (simParams->getPrivilegedAxis() == "x")
    axis_ = 0;
  else if (simParams->getPrivilegedAxis() == "y")
    axis_ = 1;
  else if (simParams->getPrivilegedAxis() == "z")
    axis_ = 2;
}

void NPAT::evolveEtaA() {
  eta(axis_, axis_) +=
      dt2 * instaVol *
      (press(axis_, axis_) - targetPressure / Constants::pressureConvert) /
      (NkBT * tb2);
  oldEta_ = eta;
}

void NPAT::evolveEtaB() {
  prevEta_ = eta;
  eta(axis_, axis_) =
      oldEta_(axis_, axis_) +
      dt2 * instaVol *
          (press(axis_, axis_) - targetPressure / Constants::pressureConvert) /
          (NkBT * tb2);
}

void NPAT::calcVelScale() {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      vScale_(i, j) = eta(i, j);

      if (i == j) {
        vScale_(i, j) += thermostat.first;
      }
    }
  }
}

void NPAT::getVelScaleA(Vector3d& sc, const Vector3d& vel) {
  sc = vScale_ * vel;
}

void NPAT::getVelScaleB(Vector3d& sc, int index) {
  sc = vScale_ * oldVel[index];
}

void NPAT::getPosScale(const Vector3d& pos, const Vector3d& COM, int index,
                       Vector3d& sc) {
  /**@todo */
  Vector3d rj = (oldPos[index] + pos) / (RealType)2.0 - COM;
  sc = eta * rj;
}

void NPAT::scaleSimBox() {
  Mat3x3d scaleMat;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      scaleMat(i, j) = 0.0;
      if (i == j) {
        scaleMat(i, j) = 1.0;
      }
    }
  }

  scaleMat(axis_, axis_) = exp(dt * eta(axis_, axis_));
  Mat3x3d hmat = snap->getHmat();
  hmat = hmat * scaleMat;
  snap->setHmat(hmat);
}

bool NPAT::etaConverged() {
  int i;
  RealType diffEta, sumEta;

  sumEta = 0;
  for (i = 0; i < 3; i++) {
    sumEta += pow(prevEta_(i, i) - eta(i, i), 2);
  }

  diffEta = sqrt(sumEta / 3.0);

  return (diffEta <= etaTolerance);
}

RealType NPAT::calcConservedQuantity() {
  thermostat = snap->getThermostat();
  loadEta();

  // We need NkBT a lot, so just set it here: This is the RAW number
  // of integrableObjects, so no subtraction or addition of constraints or
  // orientational degrees of freedom:
  NkBT = info_->getNGlobalIntegrableObjects() * Constants::kB * targetTemp;

  // fkBT is used because the thermostat operates on more degrees of freedom
  // than the barostat (when there are particles with orientational degrees
  // of freedom).
  fkBT = info_->getNdf() * Constants::kB * targetTemp;

  RealType conservedQuantity;
  RealType totalEnergy;
  RealType thermostat_kinetic;
  RealType thermostat_potential;
  RealType barostat_kinetic;
  RealType barostat_potential;
  RealType trEta;

  totalEnergy = thermo.getTotalEnergy();

  thermostat_kinetic = fkBT * tt2 * thermostat.first * thermostat.first /
                       (2.0 * Constants::energyConvert);

  thermostat_potential = fkBT * thermostat.second / Constants::energyConvert;

  SquareMatrix<RealType, 3> tmp = eta.transpose() * eta;
  trEta = tmp.trace();

  barostat_kinetic = NkBT * tb2 * trEta / (2.0 * Constants::energyConvert);

  barostat_potential =
      (targetPressure * thermo.getVolume() / Constants::pressureConvert) /
      Constants::energyConvert;

  conservedQuantity = totalEnergy + thermostat_kinetic + thermostat_potential +
                      barostat_kinetic + barostat_potential;

  return conservedQuantity;
}

void NPAT::loadEta() {
  eta = snap->getBarostat();

  // if (!eta.isDiagonal()) {
  //    sprintf( painCave.errMsg,
  //             "NPAT error: the diagonal elements of eta matrix are not the
  //             same or etaMat is not a diagonal matrix");
  //    painCave.isFatal = 1;
  //    simError();
  //}
}

void NPAT::saveEta() { snap->setBarostat(eta); }

}  // namespace OpenMD
