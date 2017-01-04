/*
 * Copyright (c) 2012 The University of Notre Dame. All Rights Reserved.
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
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
#include "types/EAMAdapter.hpp"
#include "utils/simError.h"
#include <cstdio>

namespace OpenMD {

  bool EAMAdapter::isEAM() {
    return at_->hasProperty(EAMtypeID);
  }
  
  EAMAtypeParameters* EAMAdapter::getEAMParam() {
    
    if (!isEAM()) {
      sprintf( painCave.errMsg,               
               "EAMAdapter::getEAMParam was passed an atomType (%s)\n"
               "\tthat does not appear to be an EAM atom.\n",
               at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
    
    GenericData* data = at_->getPropertyByName(EAMtypeID);
    if (data == NULL) {
      sprintf( painCave.errMsg, 
               "EAMAdapter::getEAMParam could not find EAM\n"
               "\tparameters for atomType %s.\n", at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError(); 
    }
    
    EAMAtypeData* eamData = dynamic_cast<EAMAtypeData*>(data);
    if (eamData == NULL) {
      sprintf( painCave.errMsg,
               "EAMAdapter::getEAMParam could not convert\n"
               "\tGenericData to EAMAtypeData for atom type %s\n", 
               at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();          
    }
    
    return eamData->getData();
  }
  
  RealType EAMAdapter::getLatticeConstant() {    
    EAMAtypeParameters* eamParam = getEAMParam();
    return eamParam->latticeConstant;
  }

  std::string EAMAdapter::getLatticeType() {    
    EAMAtypeParameters* eamParam = getEAMParam();
    return eamParam->latticeType;
  }

  int EAMAdapter::getNrho() {    
    EAMAtypeParameters* eamParam = getEAMParam();
    return eamParam->nrho;
  }

  RealType EAMAdapter::getDrho() {    
    EAMAtypeParameters* eamParam = getEAMParam();
    return eamParam->drho;
  }

  int EAMAdapter::getNr() {    
    EAMAtypeParameters* eamParam = getEAMParam();
    return eamParam->nr;
  }

  RealType EAMAdapter::getDr() {    
    EAMAtypeParameters* eamParam = getEAMParam();
    return eamParam->dr;
  }

  RealType EAMAdapter::getRcut() {    
    EAMAtypeParameters* eamParam = getEAMParam();
    return eamParam->rcut;
  }

  RealType EAMAdapter::getRe() {    
    EAMAtypeParameters* eamParam = getEAMParam();
    return eamParam->re;
  }
  RealType EAMAdapter::get_fe() {    
    EAMAtypeParameters* eamParam = getEAMParam();
    return eamParam->fe;
  }
  RealType EAMAdapter::getRhoe() {    
    EAMAtypeParameters* eamParam = getEAMParam();
    return eamParam->rhoe;
  }
  RealType EAMAdapter::getAlpha() {    
    EAMAtypeParameters* eamParam = getEAMParam();
    return eamParam->alpha;
  }
  RealType EAMAdapter::getBeta() {    
    EAMAtypeParameters* eamParam = getEAMParam();
    return eamParam->beta;
  }
  RealType EAMAdapter::getA() {    
    EAMAtypeParameters* eamParam = getEAMParam();
    return eamParam->A;
  }
  RealType EAMAdapter::getB() {    
    EAMAtypeParameters* eamParam = getEAMParam();
    return eamParam->B;
  }
  RealType EAMAdapter::getKappa() {    
    EAMAtypeParameters* eamParam = getEAMParam();
    return eamParam->kappa;
  }
  RealType EAMAdapter::getLambda() {    
    EAMAtypeParameters* eamParam = getEAMParam();
    return eamParam->lambda;
  }
  std::vector<RealType> EAMAdapter::getFn() {    
    EAMAtypeParameters* eamParam = getEAMParam();
    return eamParam->Fn;
  }
  std::vector<RealType> EAMAdapter::getF() {    
    EAMAtypeParameters* eamParam = getEAMParam();
    return eamParam->F;
  }
  RealType EAMAdapter::getEta() {    
    EAMAtypeParameters* eamParam = getEAMParam();
    return eamParam->eta;
  }
  RealType EAMAdapter::getFe() {    
    EAMAtypeParameters* eamParam = getEAMParam();
    return eamParam->Fe;
  }
  bool EAMAdapter::hasSplines() {    
    EAMAtypeParameters* eamParam = getEAMParam();
    return eamParam->hasSplines;
  }
  
  CubicSpline* EAMAdapter::getZSpline() {    
    EAMAtypeParameters* eamParam = getEAMParam();
    int nr = eamParam->nr;
    RealType dr = eamParam->dr;
    vector<RealType> rvals;
    
    for (int i = 0; i < nr; i++) rvals.push_back(RealType(i) * dr);
      
    CubicSpline* cs = new CubicSpline();
    cs->addPoints(rvals, eamParam->Z);
    return cs;
  }

  CubicSpline* EAMAdapter::getRhoSpline() {    
    EAMAtypeParameters* eamParam = getEAMParam();
    int nr = eamParam->nr;
    RealType dr = eamParam->dr;
    vector<RealType> rvals;
    
    for (int i = 0; i < nr; i++) rvals.push_back(RealType(i) * dr);
      
    CubicSpline* cs = new CubicSpline();
    cs->addPoints(rvals, eamParam->rho);
    return cs;
  }

  CubicSpline* EAMAdapter::getFSpline() {    
    EAMAtypeParameters* eamParam = getEAMParam();
    int nrho = eamParam->nrho;
    RealType drho = eamParam->drho;
    vector<RealType> rhovals;
    
    for (int i = 0; i < nrho; i++) {
      rhovals.push_back(RealType(i) * drho);
    }
      
    CubicSpline* cs = new CubicSpline();
    cs->addPoints(rhovals, eamParam->F);
    return cs;
  }

  void EAMAdapter::makeEAM(RealType latticeConstant, std::string latticeType,
                           int nrho, RealType drho, 
                           int nr, RealType dr, RealType rcut,
                           vector<RealType> Z, vector<RealType> rho, 
                           vector<RealType> F, bool hasSplines) {

    if (isEAM()){
      at_->removeProperty(EAMtypeID);
    }

    EAMAtypeParameters* eamParam = new EAMAtypeParameters();
    eamParam->latticeConstant = latticeConstant;
    eamParam->latticeType = latticeType;
    eamParam->nrho = nrho;
    eamParam->drho = drho;
    eamParam->nr = nr;
    eamParam->dr = dr;
    eamParam->rcut = rcut;
    eamParam->Z = Z;
    eamParam->rho = rho;
    eamParam->F = F;
    eamParam->hasSplines = hasSplines;
    
    at_->addProperty(new EAMAtypeData(EAMtypeID, eamParam));
  }

  void EAMAdapter::makeEAM(RealType latticeConstant, std::string latticeType,
                           RealType re,
                           RealType fe,
                           RealType rhoe,
                           RealType alpha,
                           RealType beta,
                           RealType A,
                           RealType B,
                           RealType kappa,
                           RealType lambda,
                           std::vector<RealType> Fn,
                           std::vector<RealType> F,
                           RealType eta,
                           RealType Fe,
                           bool hasSplines){ 
    
    if (isEAM()){
      at_->removeProperty(EAMtypeID);
    }
    
    EAMAtypeParameters* eamParam = new EAMAtypeParameters();
    eamParam->latticeConstant = latticeConstant;
    eamParam->latticeType = latticeType;
    eamParam->re = re;
    eamParam->fe = fe;
    eamParam->rhoe = rhoe;
    eamParam->alpha = alpha;
    eamParam->beta = beta;
    eamParam->A = A;
    eamParam->B = B;
    eamParam->kappa = kappa;
    eamParam->lambda = lambda;
    eamParam->Fn = Fn;
    eamParam->F = F;
    eamParam->eta = eta;
    eamParam->Fe = Fe;
    eamParam->hasSplines = hasSplines;

    at_->addProperty(new EAMAtypeData(EAMtypeID, eamParam));
  }

}
