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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
#include "types/LennardJonesAdapter.hpp"
#include "utils/simError.h"

namespace OpenMD {

  bool LennardJonesAdapter::isLennardJones() {
    return at_->hasProperty(LJtypeID);
  }
  
  LJAtypeParameters* LennardJonesAdapter::getLJParam() {
    
    if (!isLennardJones()) {
      sprintf( painCave.errMsg,               
               "LennardJonesAdapter::getLJParam was passed an atomType (%s)\n"
               "\tthat does not appear to be a Lennard-Jones atom.\n",
               at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
    
    GenericData* data = at_->getPropertyByName(LJtypeID);
    if (data == NULL) {
      sprintf( painCave.errMsg, 
               "LennardJonesAdapter::getLJParam could not find Lennard-Jones\n"
               "\tparameters for atomType %s.\n", at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError(); 
    }
    
    LJAtypeData* ljData = dynamic_cast<LJAtypeData*>(data);
    if (ljData == NULL) {
      sprintf( painCave.errMsg,
               "LennardJonesAdapter::getLJParam could not convert\n"
               "\tGenericData to LJAtypeData for atom type %s\n", 
               at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();          
    }
    
    return ljData->getData();
  }
  
  RealType LennardJonesAdapter::getSigma() {    
    LJAtypeParameters* ljParam = getLJParam();
    return ljParam->sigma;
  }
  
  RealType LennardJonesAdapter::getEpsilon() {    
    LJAtypeParameters* ljParam = getLJParam();
    return ljParam->epsilon;
  }
  
  bool LennardJonesAdapter::isSoft() {    
    LJAtypeParameters* ljParam = getLJParam();
    return ljParam->isSoft;
  }
  
  void LennardJonesAdapter::makeLennardJones(RealType sigma, 
                                             RealType epsilon, bool isSoft){
    if (isLennardJones()){
      at_->removeProperty(LJtypeID);
    }

    LJAtypeParameters* ljParam = new LJAtypeParameters();
    ljParam->epsilon = epsilon;
    ljParam->sigma = sigma;
    ljParam->isSoft = isSoft;
    
    at_->addProperty(new LJAtypeData(LJtypeID, ljParam));
  }
}
