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

#include "NPTi.hpp"

#include "brains/SimInfo.hpp"
#include "brains/Thermo.hpp"
#include "integrators/NPT.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Constants.hpp"
#include "utils/simError.h"

namespace OpenMD {

  // Basic isotropic thermostating and barostating via the Melchionna
  // modification of the Hoover algorithm:
  //
  //    Melchionna, S., Ciccotti, G., and Holian, B. L., 1993,
  //       Molec. Phys., 78, 533.
  //
  //           and
  //
  //    Hoover, W. G., 1986, Phys. Rev. A, 34, 2499.

  NPTi::NPTi(SimInfo* info) : NPT(info) {}

  void NPTi::evolveEtaA() {
    eta += dt2 * (instaVol * (instaPress - targetPressure) /
                  (Constants::pressureConvert * NkBT * tb2));
    oldEta = eta;
  }

  void NPTi::evolveEtaB() {
    prevEta = eta;
    eta     = oldEta + dt2 * (instaVol * (instaPress - targetPressure) /
                          (Constants::pressureConvert * NkBT * tb2));
  }

  void NPTi::calcVelScale() { vScale = thermostat.first + eta; }

  void NPTi::getVelScaleA(Vector3d& sc, const Vector3d& vel) {
    sc = vel * vScale;
  }

  void NPTi::getVelScaleB(Vector3d& sc, int index) {
    sc = oldVel[index] * vScale;
  }

  void NPTi::getPosScale(const Vector3d& pos, const Vector3d& COM, int index,
                         Vector3d& sc) {
    /**@todo*/
    sc = (oldPos[index] + pos) / (RealType)2.0 - COM;
    sc *= eta;
  }

  void NPTi::scaleSimBox() {
    RealType scaleFactor;

    scaleFactor = exp(dt * eta);

    if ((scaleFactor > 1.1) || (scaleFactor < 0.9)) {
      sprintf(painCave.errMsg,
              "NPTi error: Attempting a Box scaling of more than 10 percent"
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

  bool NPTi::etaConverged() { return (fabs(prevEta - eta) <= etaTolerance); }

  RealType NPTi::calcConservedQuantity() {
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
    RealType Energy;
    RealType thermostat_kinetic;
    RealType thermostat_potential;
    RealType barostat_kinetic;
    RealType barostat_potential;

    Energy = thermo.getTotalEnergy();

    thermostat_kinetic = fkBT * tt2 * thermostat.first * thermostat.first /
                         (2.0 * Constants::energyConvert);

    thermostat_potential = fkBT * thermostat.second / Constants::energyConvert;

    barostat_kinetic =
        3.0 * NkBT * tb2 * eta * eta / (2.0 * Constants::energyConvert);

    barostat_potential =
        (targetPressure * thermo.getVolume() / Constants::pressureConvert) /
        Constants::energyConvert;

    conservedQuantity = Energy + thermostat_kinetic + thermostat_potential +
                        barostat_kinetic + barostat_potential;

    return conservedQuantity;
  }

  void NPTi::loadEta() {
    Mat3x3d etaMat = snap->getBarostat();
    eta            = etaMat(0, 0);
    // if (fabs(etaMat(1,1) - eta) >= OpenMD::epsilon || fabs(etaMat(1,1) - eta)
    // >= OpenMD::epsilon || !etaMat.isDiagonal()) {
    //    sprintf( painCave.errMsg,
    //             "NPTi error: the diagonal elements of  eta matrix are not the
    //             same or etaMat is not a diagonal matrix");
    //    painCave.isFatal = 1;
    //    simError();
    //}
  }

  void NPTi::saveEta() {
    Mat3x3d etaMat(0.0);
    etaMat(0, 0) = eta;
    etaMat(1, 1) = eta;
    etaMat(2, 2) = eta;
    snap->setBarostat(etaMat);
  }
}  // namespace OpenMD
