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

#include "types/EAMAdapter.hpp"

#include <cmath>
#include <cstdio>
#include <memory>

#include "utils/CaseConversion.hpp"
#include "utils/simError.h"

namespace OpenMD {

  bool EAMAdapter::isEAM() { return at_->hasProperty(EAMtypeID); }

  EAMParameters EAMAdapter::getEAMParam() {
    if (!isEAM()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "EAMAdapter::getEAMParam was passed an atomType (%s)\n"
               "\tthat does not appear to be an EAM atom.\n",
               at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    std::shared_ptr<GenericData> data = at_->getPropertyByName(EAMtypeID);
    if (data == nullptr) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "EAMAdapter::getEAMParam could not find EAM\n"
               "\tparameters for atomType %s.\n",
               at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    std::shared_ptr<EAMData> eamData = std::dynamic_pointer_cast<EAMData>(data);
    if (eamData == nullptr) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "EAMAdapter::getEAMParam could not convert\n"
               "\tGenericData to EAMData for atom type %s\n",
               at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    return eamData->getData();
  }

  FuncflParameters EAMAdapter::getFuncflParam() {
    if (!isEAM()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "EAMAdapter::getFuncflParam was passed an atomType (%s)\n"
               "\tthat does not appear to be an EAM atom.\n",
               at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    std::shared_ptr<GenericData> data = at_->getPropertyByName(FuncflTypeID);
    if (data == nullptr) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "EAMAdapter::getFuncflParam could not find Funcfl\n"
               "\tparameters for atomType %s.\n",
               at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    std::shared_ptr<FuncflData> funcflData =
        std::dynamic_pointer_cast<FuncflData>(data);
    if (funcflData == nullptr) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "EAMAdapter::getFuncflParam could not convert\n"
               "\tGenericData to FuncflData for atom type %s\n",
               at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    return funcflData->getData();
  }

  ZhouParameters EAMAdapter::getZhouParam() {
    if (!isEAM()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "EAMAdapter::getZhouParam was passed an atomType (%s)\n"
               "\tthat does not appear to be an EAM atom.\n",
               at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    std::shared_ptr<GenericData> data = at_->getPropertyByName(ZhouTypeID);
    if (data == nullptr) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "EAMAdapter::getZhouParam could not find Zhou\n"
               "\tparameters for atomType %s.\n",
               at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    std::shared_ptr<ZhouData> zhouData =
        std::dynamic_pointer_cast<ZhouData>(data);
    if (zhouData == nullptr) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "EAMAdapter::getZhouParam could not convert\n"
               "\tGenericData to ZhouData for atom type %s\n",
               at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    return zhouData->getData();
  }

  EAMType EAMAdapter::getEAMType() {
    EAMParameters eamParam = getEAMParam();
    return eamParam.eamType;
  }

  std::string EAMAdapter::getLatticeType() {
    EAMParameters eamParam = getEAMParam();
    return eamParam.latticeType;
  }

  RealType EAMAdapter::getLatticeConstant() {
    EAMParameters eamParam = getEAMParam();
    return eamParam.latticeConstant;
  }

  CubicSplinePtr EAMAdapter::getZSpline() {
    FuncflParameters funcflParam = getFuncflParam();
    int nr                       = funcflParam.nr;
    RealType dr                  = funcflParam.dr;
    std::vector<RealType> rvals;

    for (int i = 0; i < nr; i++)
      rvals.push_back(RealType(i) * dr);

    CubicSplinePtr cs {std::make_shared<CubicSpline>()};
    cs->addPoints(rvals, funcflParam.Z);
    return cs;
  }

  CubicSplinePtr EAMAdapter::getRhoSpline() {
    FuncflParameters funcflParam = getFuncflParam();
    int nr                       = funcflParam.nr;
    RealType dr                  = funcflParam.dr;
    std::vector<RealType> rvals;

    for (int i = 0; i < nr; i++)
      rvals.push_back(RealType(i) * dr);

    CubicSplinePtr cs {std::make_shared<CubicSpline>()};
    cs->addPoints(rvals, funcflParam.rho);
    return cs;
  }

  CubicSplinePtr EAMAdapter::getFSpline() {
    FuncflParameters funcflParam = getFuncflParam();
    int nrho                     = funcflParam.nrho;
    RealType drho                = funcflParam.drho;
    std::vector<RealType> rhovals;

    for (int i = 0; i < nrho; i++) {
      rhovals.push_back(RealType(i) * drho);
    }

    CubicSplinePtr cs {std::make_shared<CubicSpline>()};
    cs->addPoints(rhovals, funcflParam.F);
    return cs;
  }

  void EAMAdapter::makeFuncfl(RealType latticeConstant, std::string latticeType,
                              int nrho, RealType drho, int nr, RealType dr,
                              RealType rcut, vector<RealType> Z,
                              vector<RealType> rho, vector<RealType> F) {
    if (isEAM()) {
      at_->removeProperty(EAMtypeID);
      at_->removeProperty(FuncflTypeID);
    }

    EAMParameters eamParam {};
    FuncflParameters funcflParam {};

    eamParam.eamType         = eamFuncfl;
    eamParam.latticeConstant = latticeConstant;
    eamParam.latticeType     = latticeType;

    funcflParam.nrho = nrho;
    funcflParam.drho = drho;
    funcflParam.nr   = nr;
    funcflParam.dr   = dr;
    funcflParam.rcut = rcut;
    funcflParam.Z    = Z;
    funcflParam.rho  = rho;
    funcflParam.F    = F;

    at_->addProperty(std::make_shared<EAMData>(EAMtypeID, eamParam));
    at_->addProperty(std::make_shared<FuncflData>(FuncflTypeID, funcflParam));
  }

  void EAMAdapter::makeZhou2001(std::string latticeType, RealType re,
                                RealType fe, RealType rhoe, RealType alpha,
                                RealType beta, RealType A, RealType B,
                                RealType kappa, RealType lambda,
                                std::vector<RealType> Fn,
                                std::vector<RealType> F, RealType eta,
                                RealType Fe) {
    if (isEAM()) {
      at_->removeProperty(EAMtypeID);
      at_->removeProperty(ZhouTypeID);
    }

    EAMParameters eamParam {};
    ZhouParameters zhouParam {};

    eamParam.eamType = eamZhou2001;

    toUpper(latticeType);
    eamParam.latticeType = latticeType;
    // default to FCC if we don't specify HCP or BCC:
    if (latticeType == "HCP")
      eamParam.latticeConstant = re;
    else if (latticeType == "BCC")
      eamParam.latticeConstant = 2.0 * re / sqrt(3.0);
    else
      eamParam.latticeConstant = 2.0 * re / sqrt(2.0);

    zhouParam.re     = re;
    zhouParam.fe     = fe;
    zhouParam.rhoe   = rhoe;
    zhouParam.alpha  = alpha;
    zhouParam.beta   = beta;
    zhouParam.A      = A;
    zhouParam.B      = B;
    zhouParam.kappa  = kappa;
    zhouParam.lambda = lambda;
    zhouParam.Fn     = Fn;
    zhouParam.F      = F;
    zhouParam.eta    = eta;
    zhouParam.Fe     = Fe;

    at_->addProperty(std::make_shared<EAMData>(EAMtypeID, eamParam));
    at_->addProperty(std::make_shared<ZhouData>(ZhouTypeID, zhouParam));
  }

  void EAMAdapter::makeZhou2004(std::string latticeType, RealType re,
                                RealType fe, RealType rhoe, RealType rhos,
                                RealType alpha, RealType beta, RealType A,
                                RealType B, RealType kappa, RealType lambda,
                                std::vector<RealType> Fn,
                                std::vector<RealType> F, RealType eta,
                                RealType Fe, RealType rhol, RealType rhoh) {
    if (isEAM()) {
      at_->removeProperty(EAMtypeID);
      at_->removeProperty(ZhouTypeID);
    }

    EAMParameters eamParam {};
    ZhouParameters zhouParam {};

    eamParam.eamType = eamZhou2004;

    toUpper(latticeType);
    eamParam.latticeType = latticeType;
    // default to FCC if we don't specify HCP or BCC:
    if (latticeType == "HCP")
      eamParam.latticeConstant = re;
    else if (latticeType == "BCC")
      eamParam.latticeConstant = 2.0 * re / sqrt(3.0);
    else
      eamParam.latticeConstant = 2.0 * re / sqrt(2.0);

    zhouParam.re     = re;
    zhouParam.fe     = fe;
    zhouParam.rhoe   = rhoe;
    zhouParam.alpha  = alpha;
    zhouParam.beta   = beta;
    zhouParam.A      = A;
    zhouParam.B      = B;
    zhouParam.kappa  = kappa;
    zhouParam.lambda = lambda;
    zhouParam.Fn     = Fn;
    zhouParam.F      = F;
    zhouParam.eta    = eta;
    zhouParam.Fe     = Fe;
    zhouParam.rhos   = rhos;
    zhouParam.rhol   = rhol;
    zhouParam.rhoh   = rhoh;

    at_->addProperty(std::make_shared<EAMData>(EAMtypeID, eamParam));
    at_->addProperty(std::make_shared<ZhouData>(ZhouTypeID, zhouParam));
  }

  void EAMAdapter::makeZhou2005(std::string latticeType, RealType re,
                                RealType fe, RealType rhoe, RealType rhos,
                                RealType alpha, RealType beta, RealType A,
                                RealType B, RealType kappa, RealType lambda,
                                std::vector<RealType> Fn,
                                std::vector<RealType> F, RealType F3plus,
                                RealType F3minus, RealType eta, RealType Fe) {
    if (isEAM()) {
      at_->removeProperty(EAMtypeID);
      at_->removeProperty(ZhouTypeID);
    }

    EAMParameters eamParam {};
    ZhouParameters zhouParam {};

    eamParam.eamType = eamZhou2005;

    toUpper(latticeType);
    eamParam.latticeType = latticeType;
    // default to FCC if we don't specify HCP or BCC:
    if (latticeType == "HCP")
      eamParam.latticeConstant = re;
    else if (latticeType == "BCC")
      eamParam.latticeConstant = 2.0 * re / sqrt(3.0);
    else
      eamParam.latticeConstant = 2.0 * re / sqrt(2.0);

    zhouParam.re      = re;
    zhouParam.fe      = fe;
    zhouParam.rhoe    = rhoe;
    zhouParam.alpha   = alpha;
    zhouParam.beta    = beta;
    zhouParam.A       = A;
    zhouParam.B       = B;
    zhouParam.kappa   = kappa;
    zhouParam.lambda  = lambda;
    zhouParam.Fn      = Fn;
    zhouParam.F       = F;
    zhouParam.F3plus  = F3plus;
    zhouParam.F3minus = F3minus;
    zhouParam.eta     = eta;
    zhouParam.Fe      = Fe;
    zhouParam.rhos    = rhos;

    at_->addProperty(std::make_shared<EAMData>(EAMtypeID, eamParam));
    at_->addProperty(std::make_shared<ZhouData>(ZhouTypeID, zhouParam));
  }

  void EAMAdapter::makeZhou2005Oxygen(RealType re, RealType fe, RealType alpha,
                                      RealType beta, RealType A, RealType B,
                                      RealType kappa, RealType lambda,
                                      RealType gamma, RealType nu,
                                      std::vector<RealType> OrhoLimits,
                                      std::vector<RealType> OrhoE,
                                      std::vector<std::vector<RealType>> OF) {
    if (isEAM()) {
      at_->removeProperty(EAMtypeID);
      at_->removeProperty(ZhouTypeID);
    }

    EAMParameters eamParam {};
    ZhouParameters zhouParam {};

    eamParam.eamType = eamZhou2005Oxygen;

    eamParam.latticeConstant = re;

    zhouParam.re         = re;
    zhouParam.fe         = fe;
    zhouParam.alpha      = alpha;
    zhouParam.beta       = beta;
    zhouParam.A          = A;
    zhouParam.B          = B;
    zhouParam.kappa      = kappa;
    zhouParam.lambda     = lambda;
    zhouParam.gamma      = gamma;
    zhouParam.nu         = nu;
    zhouParam.OrhoLimits = OrhoLimits;
    zhouParam.OrhoE      = OrhoE;
    zhouParam.OF         = OF;

    at_->addProperty(std::make_shared<EAMData>(EAMtypeID, eamParam));
    at_->addProperty(std::make_shared<ZhouData>(ZhouTypeID, zhouParam));
  }

  void EAMAdapter::makeZhouRose(RealType re, RealType fe, RealType rhoe,
                                RealType alpha, RealType beta, RealType A,
                                RealType B, RealType kappa, RealType lambda,
                                RealType F0) {
    if (isEAM()) {
      at_->removeProperty(EAMtypeID);
      at_->removeProperty(ZhouTypeID);
    }

    EAMParameters eamParam {};
    ZhouParameters zhouParam {};

    eamParam.eamType         = eamZhouRose;
    eamParam.latticeConstant = re;
    zhouParam.re             = re;
    zhouParam.fe             = fe;
    zhouParam.rhoe           = rhoe;
    zhouParam.alpha          = alpha;
    zhouParam.beta           = beta;
    zhouParam.A              = A;
    zhouParam.B              = B;
    zhouParam.kappa          = kappa;
    zhouParam.lambda         = lambda;
    zhouParam.F0             = F0;

    at_->addProperty(std::make_shared<EAMData>(EAMtypeID, eamParam));
    at_->addProperty(std::make_shared<ZhouData>(ZhouTypeID, zhouParam));
  }

  void EAMAdapter::makeOxygenFuncfl(RealType re, RealType fe, RealType alpha,
                                    RealType beta, RealType A, RealType B,
                                    RealType kappa, RealType lambda,
                                    RealType drho, RealType nrho,
                                    std::vector<RealType> F) {
    if (isEAM()) {
      at_->removeProperty(EAMtypeID);
      at_->removeProperty(ZhouTypeID);
      at_->removeProperty(FuncflTypeID);
    }

    EAMParameters eamParam {};
    FuncflParameters funcflParam {};
    ZhouParameters zhouParam {};

    eamParam.eamType         = eamOxygenFuncfl;
    eamParam.latticeConstant = re;
    zhouParam.re             = re;
    zhouParam.fe             = fe;
    zhouParam.alpha          = alpha;
    zhouParam.beta           = beta;
    zhouParam.A              = A;
    zhouParam.B              = B;
    zhouParam.kappa          = kappa;
    zhouParam.lambda         = lambda;

    funcflParam.F    = F;
    funcflParam.drho = drho;
    funcflParam.nrho = nrho;

    at_->addProperty(std::make_shared<EAMData>(EAMtypeID, eamParam));
    at_->addProperty(std::make_shared<ZhouData>(ZhouTypeID, zhouParam));
    at_->addProperty(std::make_shared<FuncflData>(FuncflTypeID, funcflParam));
  }

  int EAMAdapter::getNrho() {
    FuncflParameters funcflParam = getFuncflParam();
    return funcflParam.nrho;
  }

  RealType EAMAdapter::getDrho() {
    FuncflParameters funcflParam = getFuncflParam();
    return funcflParam.drho;
  }

  int EAMAdapter::getNr() {
    FuncflParameters funcflParam = getFuncflParam();
    return funcflParam.nr;
  }

  RealType EAMAdapter::getDr() {
    FuncflParameters funcflParam = getFuncflParam();
    return funcflParam.dr;
  }

  RealType EAMAdapter::getRcut() {
    FuncflParameters funcflParam = getFuncflParam();
    return funcflParam.rcut;
  }

  RealType EAMAdapter::getRe() {
    ZhouParameters zhouParam = getZhouParam();
    return zhouParam.re;
  }
  RealType EAMAdapter::get_fe() {
    ZhouParameters zhouParam = getZhouParam();
    return zhouParam.fe;
  }
  RealType EAMAdapter::getRhoe() {
    ZhouParameters zhouParam = getZhouParam();
    return zhouParam.rhoe;
  }
  RealType EAMAdapter::getAlpha() {
    ZhouParameters zhouParam = getZhouParam();
    return zhouParam.alpha;
  }
  RealType EAMAdapter::getBeta() {
    ZhouParameters zhouParam = getZhouParam();
    return zhouParam.beta;
  }
  RealType EAMAdapter::getA() {
    ZhouParameters zhouParam = getZhouParam();
    return zhouParam.A;
  }
  RealType EAMAdapter::getB() {
    ZhouParameters zhouParam = getZhouParam();
    return zhouParam.B;
  }
  RealType EAMAdapter::getKappa() {
    ZhouParameters zhouParam = getZhouParam();
    return zhouParam.kappa;
  }
  RealType EAMAdapter::getLambda() {
    ZhouParameters zhouParam = getZhouParam();
    return zhouParam.lambda;
  }
  RealType EAMAdapter::getGamma() {
    ZhouParameters zhouParam = getZhouParam();
    return zhouParam.gamma;
  }
  RealType EAMAdapter::getNu() {
    ZhouParameters zhouParam = getZhouParam();
    return zhouParam.nu;
  }
  std::vector<RealType> EAMAdapter::getFn() {
    ZhouParameters zhouParam = getZhouParam();
    return zhouParam.Fn;
  }
  std::vector<RealType> EAMAdapter::getF() {
    ZhouParameters zhouParam = getZhouParam();
    return zhouParam.F;
  }
  RealType EAMAdapter::getF3plus() {
    ZhouParameters zhouParam = getZhouParam();
    return zhouParam.F3plus;
  }
  RealType EAMAdapter::getF3minus() {
    ZhouParameters zhouParam = getZhouParam();
    return zhouParam.F3minus;
  }

  RealType EAMAdapter::getEta() {
    ZhouParameters zhouParam = getZhouParam();
    return zhouParam.eta;
  }

  RealType EAMAdapter::getFe() {
    ZhouParameters zhouParam = getZhouParam();
    return zhouParam.Fe;
  }

  RealType EAMAdapter::getRhos() {
    ZhouParameters zhouParam = getZhouParam();
    return zhouParam.rhos;
  }

  RealType EAMAdapter::getRhol() {
    ZhouParameters zhouParam = getZhouParam();
    return zhouParam.rhol;
  }

  RealType EAMAdapter::getRhoh() {
    ZhouParameters zhouParam = getZhouParam();
    return zhouParam.rhoh;
  }

  std::vector<RealType> EAMAdapter::getOrhoLimits() {
    ZhouParameters zhouParam = getZhouParam();
    return zhouParam.OrhoLimits;
  }
  std::vector<RealType> EAMAdapter::getOrhoE() {
    ZhouParameters zhouParam = getZhouParam();
    return zhouParam.OrhoE;
  }

  std::vector<std::vector<RealType>> EAMAdapter::getOF() {
    ZhouParameters zhouParam = getZhouParam();
    return zhouParam.OF;
  }
  RealType EAMAdapter::getF0() {
    ZhouParameters zhouParam = getZhouParam();
    return zhouParam.F0;
  }
}  // namespace OpenMD
