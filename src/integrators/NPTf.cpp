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

#include "integrators/NPTf.hpp"

#include "brains/SimInfo.hpp"
#include "brains/Thermo.hpp"
#include "integrators/IntegratorCreator.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Constants.hpp"
#include "utils/simError.h"

namespace OpenMD {

  // Basic non-isotropic thermostating and barostating via the Melchionna
  // modification of the Hoover algorithm:
  //
  //    Melchionna, S., Ciccotti, G., and Holian, B. L., 1993,
  //       Molec. Phys., 78, 533.
  //
  //           and
  //
  //    Hoover, W. G., 1986, Phys. Rev. A, 34, 2499.

  void NPTf::evolveEtaA() {
    int i, j;

    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        if (i == j) {
          eta(i, j) +=
              dt2 * instaVol *
              (press(i, j) - targetPressure / Constants::pressureConvert) /
              (NkBT * tb2);
        } else {
          eta(i, j) += dt2 * instaVol * press(i, j) / (NkBT * tb2);
        }
      }
    }

    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        oldEta(i, j) = eta(i, j);
      }
    }
  }

  void NPTf::evolveEtaB() {
    int i;
    int j;

    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        prevEta(i, j) = eta(i, j);
      }
    }

    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        if (i == j) {
          eta(i, j) =
              oldEta(i, j) +
              dt2 * instaVol *
                  (press(i, j) - targetPressure / Constants::pressureConvert) /
                  (NkBT * tb2);
        } else {
          eta(i, j) =
              oldEta(i, j) + dt2 * instaVol * press(i, j) / (NkBT * tb2);
        }
      }
    }
  }

  void NPTf::calcVelScale() {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        vScale(i, j) = eta(i, j);

        if (i == j) { vScale(i, j) += thermostat.first; }
      }
    }
  }

  void NPTf::getVelScaleA(Vector3d& sc, const Vector3d& vel) {
    sc = vScale * vel;
  }

  void NPTf::getVelScaleB(Vector3d& sc, int index) {
    sc = vScale * oldVel[index];
  }

  void NPTf::getPosScale(const Vector3d& pos, const Vector3d& COM, int index,
                         Vector3d& sc) {
    /**@todo */
    Vector3d rj = (oldPos[index] + pos) / (RealType)2.0 - COM;
    sc          = eta * rj;
  }

  void NPTf::scaleSimBox() {
    int i;
    int j;
    int k;
    Mat3x3d scaleMat;
    RealType eta2ij;
    RealType bigScale, smallScale, offDiagMax;
    Mat3x3d hm;
    Mat3x3d hmnew;

    // Scale the box after all the positions have been moved:

    // Use a taylor expansion for eta products:  Hmat = Hmat . exp(dt * etaMat)
    //  Hmat = Hmat . ( Ident + dt * etaMat  + dt^2 * etaMat*etaMat / 2)

    bigScale   = 1.0;
    smallScale = 1.0;
    offDiagMax = 0.0;

    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        // Calculate the matrix Product of the eta array (we only need
        // the ij element right now):

        eta2ij = 0.0;
        for (k = 0; k < 3; k++) {
          eta2ij += eta(i, k) * eta(k, j);
        }

        scaleMat(i, j) = 0.0;
        // identity matrix (see above):
        if (i == j) scaleMat(i, j) = 1.0;
        // Taylor expansion for the exponential truncated at second order:
        scaleMat(i, j) += dt * eta(i, j) + 0.5 * dt * dt * eta2ij;

        if (i != j)
          if (fabs(scaleMat(i, j)) > offDiagMax)
            offDiagMax = fabs(scaleMat(i, j));
      }

      if (scaleMat(i, i) > bigScale) bigScale = scaleMat(i, i);
      if (scaleMat(i, i) < smallScale) smallScale = scaleMat(i, i);
    }

    if ((bigScale > 1.01) || (smallScale < 0.99)) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "NPTf error: Attempting a Box scaling of more than 1 percent.\n"
               " Check your tauBarostat, as it is probably too small!\n\n"
               " scaleMat = [%lf\t%lf\t%lf]\n"
               "            [%lf\t%lf\t%lf]\n"
               "            [%lf\t%lf\t%lf]\n"
               "      eta = [%lf\t%lf\t%lf]\n"
               "            [%lf\t%lf\t%lf]\n"
               "            [%lf\t%lf\t%lf]\n",
               scaleMat(0, 0), scaleMat(0, 1), scaleMat(0, 2), scaleMat(1, 0),
               scaleMat(1, 1), scaleMat(1, 2), scaleMat(2, 0), scaleMat(2, 1),
               scaleMat(2, 2), eta(0, 0), eta(0, 1), eta(0, 2), eta(1, 0),
               eta(1, 1), eta(1, 2), eta(2, 0), eta(2, 1), eta(2, 2));
      painCave.isFatal = 1;
      simError();
    } else if (offDiagMax > 0.01) {
      snprintf(
          painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
          "NPTf error: Attempting an off-diagonal Box scaling of more than 1 "
          "percent.\n"
          " Check your tauBarostat, as it is probably too small!\n\n"
          " scaleMat = [%lf\t%lf\t%lf]\n"
          "            [%lf\t%lf\t%lf]\n"
          "            [%lf\t%lf\t%lf]\n"
          "      eta = [%lf\t%lf\t%lf]\n"
          "            [%lf\t%lf\t%lf]\n"
          "            [%lf\t%lf\t%lf]\n",
          scaleMat(0, 0), scaleMat(0, 1), scaleMat(0, 2), scaleMat(1, 0),
          scaleMat(1, 1), scaleMat(1, 2), scaleMat(2, 0), scaleMat(2, 1),
          scaleMat(2, 2), eta(0, 0), eta(0, 1), eta(0, 2), eta(1, 0), eta(1, 1),
          eta(1, 2), eta(2, 0), eta(2, 1), eta(2, 2));
      painCave.isFatal = 1;
      simError();
    } else {
      Mat3x3d hmat = snap->getHmat();
      hmat         = hmat * scaleMat;
      snap->setHmat(hmat);
    }
  }

  bool NPTf::etaConverged() {
    int i;
    RealType diffEta, sumEta;

    sumEta = 0;
    for (i = 0; i < 3; i++) {
      sumEta += pow(prevEta(i, i) - eta(i, i), 2);
    }

    diffEta = sqrt(sumEta / 3.0);

    return (diffEta <= etaTolerance);
  }

  RealType NPTf::calcConservedQuantity() {
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
    trEta                         = tmp.trace();

    barostat_kinetic = NkBT * tb2 * trEta / (2.0 * Constants::energyConvert);

    barostat_potential =
        (targetPressure * thermo.getVolume() / Constants::pressureConvert) /
        Constants::energyConvert;

    conservedQuantity = totalEnergy + thermostat_kinetic +
                        thermostat_potential + barostat_kinetic +
                        barostat_potential;

    return conservedQuantity;
  }

  void NPTf::loadEta() {
    eta = snap->getBarostat();

    // if (!eta.isDiagonal()) {
    //    snprintf( painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
    //             "NPTf error: the diagonal elements of eta matrix are not the
    //             same or etaMat is not a diagonal matrix");
    //    painCave.isFatal = 1;
    //    simError();
    //}
  }

  void NPTf::saveEta() { snap->setBarostat(eta); }

}  // namespace OpenMD
