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

#include "types/FluctuatingChargeAdapter.hpp"

#include <config.h>

#include <cstdio>
#include <memory>

#include "nonbonded/SlaterIntegrals.hpp"
#include "utils/simError.h"

namespace OpenMD {

  bool FluctuatingChargeAdapter::isFluctuatingCharge() {
    return at_->hasProperty(FQtypeID);
  }

  FluctuatingAtypeParameters
      FluctuatingChargeAdapter::getFluctuatingChargeParam() {
    if (!isFluctuatingCharge()) {
      snprintf(
          painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
          "FluctuatingChargeAdapter::getFluctuatingChargeParam was passed an "
          "atomType "
          "(%s)\n"
          "\tthat does not appear to be a fluctuating charge atom.\n",
          at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    std::shared_ptr<GenericData> data = at_->getPropertyByName(FQtypeID);
    if (data == nullptr) {
      snprintf(
          painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
          "FluctuatingChargeAdapter::getFluctuatingChargeParam could not find "
          "fluctuating charge\n"
          "\tparameters for atomType %s.\n",
          at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    std::shared_ptr<FluctuatingAtypeData> fqData =
        std::dynamic_pointer_cast<FluctuatingAtypeData>(data);
    if (fqData == nullptr) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "FluctuatingChargeAdapter::getFluctuatingChargeParam could not "
               "convert\n"
               "\tGenericData to FluctuatingAtypeData for atom type %s\n",
               at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    return fqData->getData();
  }

  bool FluctuatingChargeAdapter::isMetallic() {
    FluctuatingAtypeParameters fqParam = getFluctuatingChargeParam();
    return fqParam.isMetallic;
  }
  bool FluctuatingChargeAdapter::usesSlaterIntramolecular() {
    FluctuatingAtypeParameters fqParam = getFluctuatingChargeParam();
    return fqParam.usesSlaterIntramolecular;
  }
  RealType FluctuatingChargeAdapter::getChargeMass() {
    FluctuatingAtypeParameters fqParam = getFluctuatingChargeParam();
    return fqParam.chargeMass;
  }
  RealType FluctuatingChargeAdapter::getElectronegativity() {
    FluctuatingAtypeParameters fqParam = getFluctuatingChargeParam();
    return fqParam.electronegativity;
  }
  RealType FluctuatingChargeAdapter::getHardness() {
    FluctuatingAtypeParameters fqParam = getFluctuatingChargeParam();
    return fqParam.hardness;
  }
  int FluctuatingChargeAdapter::getSlaterN() {
    FluctuatingAtypeParameters fqParam = getFluctuatingChargeParam();
    return fqParam.slaterN;
  }
  RealType FluctuatingChargeAdapter::getNValence() {
    FluctuatingAtypeParameters fqParam = getFluctuatingChargeParam();
    return fqParam.nValence;
  }
  RealType FluctuatingChargeAdapter::getNMobile() {
    FluctuatingAtypeParameters fqParam = getFluctuatingChargeParam();
    return fqParam.nMobile;
  }
  RealType FluctuatingChargeAdapter::getSlaterZeta() {
    FluctuatingAtypeParameters fqParam = getFluctuatingChargeParam();
    return fqParam.slaterZeta;
  }
  DoublePolynomial FluctuatingChargeAdapter::getSelfPolynomial() {
    FluctuatingAtypeParameters fqParam = getFluctuatingChargeParam();
    return fqParam.vself;
  }

  void FluctuatingChargeAdapter::makeFluctuatingCharge(
      RealType chargeMass, RealType electronegativity, RealType hardness,
      int slaterN, RealType slaterZeta) {
    if (isFluctuatingCharge()) { at_->removeProperty(FQtypeID); }

    FluctuatingAtypeParameters fqParam {};
    fqParam.chargeMass               = chargeMass;
    fqParam.usesSlaterIntramolecular = true;

    fqParam.electronegativity = electronegativity;
    fqParam.hardness          = hardness;
    fqParam.slaterN           = slaterN;
    fqParam.slaterZeta        = slaterZeta;

    fqParam.vself.setCoefficient(1, electronegativity);
    fqParam.vself.setCoefficient(2, 0.5 * hardness);

    at_->addProperty(std::make_shared<FluctuatingAtypeData>(FQtypeID, fqParam));
  }

  void FluctuatingChargeAdapter::makeFluctuatingCharge(RealType chargeMass,
                                                       RealType nValence,
                                                       DoublePolynomial vs) {
    if (isFluctuatingCharge()) { at_->removeProperty(FQtypeID); }

    FluctuatingAtypeParameters fqParam {};

    fqParam.chargeMass               = chargeMass;
    fqParam.usesSlaterIntramolecular = false;

    // old-style EAMPoly has nV = nM
    fqParam.isMetallic = true;
    fqParam.nValence   = nValence;
    fqParam.nMobile    = nValence;

    fqParam.vself = vs;

    at_->addProperty(std::make_shared<FluctuatingAtypeData>(FQtypeID, fqParam));
  }

  void FluctuatingChargeAdapter::makeFluctuatingCharge(RealType chargeMass,
                                                       RealType nValence,
                                                       RealType nMobile,
                                                       DoublePolynomial vs) {
    if (isFluctuatingCharge()) { at_->removeProperty(FQtypeID); }

    FluctuatingAtypeParameters fqParam {};

    fqParam.chargeMass               = chargeMass;
    fqParam.usesSlaterIntramolecular = false;
    fqParam.isMetallic               = true;
    fqParam.nValence                 = nValence;
    fqParam.nMobile                  = nMobile;

    fqParam.vself = vs;

    at_->addProperty(std::make_shared<FluctuatingAtypeData>(FQtypeID, fqParam));
  }
}  // namespace OpenMD
