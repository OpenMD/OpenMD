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

#include "integrators/NPTsz.hpp"

#include "brains/SimInfo.hpp"
#include "brains/Thermo.hpp"
#include "integrators/IntegratorCreator.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Constants.hpp"
#include "utils/simError.h"

namespace OpenMD {

  /**
   * There is no known conserved quantity for the NPTsz integrator,
   * but we still compute the equivalent quantity from a fully
   * flexible constant pressure integrator.
   */

  NPTsz::NPTsz(SimInfo* info) : NPTf(info) {
    Globals* simParams = info_->getSimParams();

    // Default value of privilegedAxis is "z"
    if (simParams->getPrivilegedAxis() == "x")
      axis_ = 0;
    else if (simParams->getPrivilegedAxis() == "y")
      axis_ = 1;
    else if (simParams->getPrivilegedAxis() == "z")
      axis_ = 2;

    // Compute complementary axes to the privileged axis
    axis1_ = (axis_ + 1) % 3;
    axis2_ = (axis_ + 2) % 3;
  }

  RealType NPTsz::calcConservedQuantity() {
    thermostat = snap->getThermostat();
    loadEta();

    // We need NkBT a lot, so just set it here: This is the RAW number
    // of integrableObjects, so no subtraction or addition of
    // constraints or orientational degrees of freedom:
    NkBT = info_->getNGlobalIntegrableObjects() * Constants::kB * targetTemp;

    // fkBT is used because the thermostat operates on more degrees of
    // freedom than the barostat (when there are particles with
    // orientational degrees of freedom).
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

  void NPTsz::scaleSimBox() {
    int i, j;
    Mat3x3d scaleMat;
    RealType scaleFactor;
    RealType bigScale, smallScale;
    Mat3x3d hm;
    Mat3x3d hmnew;

    // Scale the box after all the positions have been moved:

    // Use a taylor expansion for eta products:
    //  Hmat = Hmat . exp(dt * etaMat)
    //  Hmat = Hmat . ( Ident + dt * etaMat  + dt^2 * etaMat*etaMat / 2)

    bigScale   = 1.0;
    smallScale = 1.0;

    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        scaleMat(i, j) = 0.0;
        if (i == j) { scaleMat(i, j) = 1.0; }
      }
    }

    // scale x & y together:
    scaleFactor =
        0.5 * (exp(dt * eta(axis1_, axis1_)) + exp(dt * eta(axis2_, axis2_)));
    scaleMat(axis1_, axis1_) = scaleFactor;
    scaleMat(axis2_, axis2_) = scaleFactor;

    bigScale   = scaleFactor;
    smallScale = scaleFactor;

    // scale z separately
    scaleFactor            = exp(dt * eta(axis_, axis_));
    scaleMat(axis_, axis_) = scaleFactor;
    if (scaleFactor > bigScale) bigScale = scaleFactor;
    if (scaleFactor < smallScale) smallScale = scaleFactor;

    if ((bigScale > 1.1) || (smallScale < 0.9)) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "NPTsz: Attempting a Box scaling of more than 10 percent.\n"
               "\tCheck your tauBarostat, as it is probably too small!\n\n"
               "\tscaleMat = [%lf\t%lf\t%lf]\n"
               "\t           [%lf\t%lf\t%lf]\n"
               "\t           [%lf\t%lf\t%lf]\n",
               scaleMat(0, 0), scaleMat(0, 1), scaleMat(0, 2), scaleMat(1, 0),
               scaleMat(1, 1), scaleMat(1, 2), scaleMat(2, 0), scaleMat(2, 1),
               scaleMat(2, 2));
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    } else {
      Mat3x3d hmat = snap->getHmat();
      hmat         = hmat * scaleMat;
      snap->setHmat(hmat);
    }
  }

  void NPTsz::loadEta() { eta = snap->getBarostat(); }
}  // namespace OpenMD
