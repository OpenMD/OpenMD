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

#include "config.h"
#include <cstdio>
#include <memory>

#include "types/FluctuatingChargeAdapter.hpp"
#include "utils/simError.h"
#include "nonbonded/SlaterIntegrals.hpp"

namespace OpenMD {

  bool FluctuatingChargeAdapter::isFluctuatingCharge() {
    return at_->hasProperty(FQtypeID);
  }

  FluctuatingAtypeParameters FluctuatingChargeAdapter::getFluctuatingChargeParam() {

    if (!isFluctuatingCharge()) {
      sprintf( painCave.errMsg,
               "FluctuatingChargeAdapter::getFluctuatingChargeParam was passed an atomType (%s)\n"
               "\tthat does not appear to be a fluctuating charge atom.\n",
               at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }

    std::shared_ptr<GenericData> data = at_->getPropertyByName(FQtypeID);
    if (data == nullptr) {
      sprintf( painCave.errMsg,
               "FluctuatingChargeAdapter::getFluctuatingChargeParam could not find fluctuating charge\n"
               "\tparameters for atomType %s.\n", at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }

    std::shared_ptr<FluctuatingAtypeData> fqData = std::dynamic_pointer_cast<FluctuatingAtypeData>(data);
    if (fqData == nullptr) {
      sprintf( painCave.errMsg,
               "FluctuatingChargeAdapter::getFluctuatingChargeParam could not convert\n"
               "\tGenericData to FluctuatingAtypeData for atom type %s\n",
               at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
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


  void FluctuatingChargeAdapter::makeFluctuatingCharge(RealType chargeMass,
                                                       RealType electronegativity,
                                                       RealType hardness,
                                                       int slaterN,
                                                       RealType slaterZeta) {

    if (isFluctuatingCharge()){
      at_->removeProperty(FQtypeID);
    }

    FluctuatingAtypeParameters fqParam {};
    fqParam.chargeMass = chargeMass;
    fqParam.usesSlaterIntramolecular = true;

    fqParam.electronegativity = electronegativity;
    fqParam.hardness = hardness;
    fqParam.slaterN = slaterN;
    fqParam.slaterZeta = slaterZeta;

    fqParam.vself.setCoefficient(1, electronegativity);
    fqParam.vself.setCoefficient(2, 0.5 * hardness);
    
    at_->addProperty(std::make_shared<FluctuatingAtypeData>(FQtypeID, fqParam));
  }
  
  void FluctuatingChargeAdapter::makeFluctuatingCharge(RealType chargeMass,
                                                       RealType nValence,
                                                       DoublePolynomial vs) {
    
    if (isFluctuatingCharge()){
      at_->removeProperty(FQtypeID);
    }
    
    FluctuatingAtypeParameters fqParam {};

    fqParam.chargeMass = chargeMass;
    fqParam.usesSlaterIntramolecular = false;

    // old-style EAMPoly has nV = nM
    fqParam.isMetallic = true;   
    fqParam.nValence = nValence;
    fqParam.nMobile = nValence;

    fqParam.vself = vs;

    at_->addProperty(std::make_shared<FluctuatingAtypeData>(FQtypeID, fqParam));
  }
  
  void FluctuatingChargeAdapter::makeFluctuatingCharge(RealType chargeMass,
                                                       RealType nValence,
                                                       RealType nMobile,
                                                       DoublePolynomial vs) {
    
    if (isFluctuatingCharge()){
      at_->removeProperty(FQtypeID);
    }
    
    FluctuatingAtypeParameters fqParam {};

    fqParam.chargeMass = chargeMass;
    fqParam.usesSlaterIntramolecular = false;
    fqParam.isMetallic = true;   
    fqParam.nValence = nValence;
    fqParam.nMobile = nMobile;

    fqParam.vself = vs;

    at_->addProperty(std::make_shared<FluctuatingAtypeData>(FQtypeID, fqParam));
  }
}
