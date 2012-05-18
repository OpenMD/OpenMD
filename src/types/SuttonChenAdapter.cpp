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
 
#include "types/SuttonChenAdapter.hpp"
#include "utils/simError.h"

namespace OpenMD {

  bool SuttonChenAdapter::isSuttonChen() {
    return at_->hasProperty(SCtypeID);
  }
  
  SCAtypeParameters* SuttonChenAdapter::getSuttonChenParam() {
    
    if (!isSuttonChen()) {
      sprintf( painCave.errMsg,               
               "SuttonChenAdapter::getSuttonChenParam was passed an atomType (%s)\n"
               "\tthat does not appear to be a Sutton-Chen atom.\n",
               at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
    
    GenericData* data = at_->getPropertyByName(SCtypeID);
    if (data == NULL) {
      sprintf( painCave.errMsg, 
               "SuttonChenAdapter::getSuttonChenParam could not find Sutton-Chen\n"
               "\tparameters for atomType %s.\n", at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError(); 
    }
    
    SCAtypeData* scData = dynamic_cast<SCAtypeData*>(data);
    if (scData == NULL) {
      sprintf( painCave.errMsg,
               "SuttonChenAdapter::getSuttonChenParam could not convert\n"
               "\tGenericData to SCAtypeData for atom type %s\n", 
               at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();          
    }
    
    return scData->getData();
  }
  
  RealType SuttonChenAdapter::getC() {    
    SCAtypeParameters* scParam = getSuttonChenParam();
    return scParam->c;
  }

  RealType SuttonChenAdapter::getM() {    
    SCAtypeParameters* scParam = getSuttonChenParam();
    return scParam->m;
  }

  RealType SuttonChenAdapter::getN() {    
    SCAtypeParameters* scParam = getSuttonChenParam();
    return scParam->n;
  }

  RealType SuttonChenAdapter::getAlpha() {    
    SCAtypeParameters* scParam = getSuttonChenParam();
    return scParam->alpha;
  }

  RealType SuttonChenAdapter::getEpsilon() {    
    SCAtypeParameters* scParam = getSuttonChenParam();
    return scParam->epsilon;
  }

    
  void SuttonChenAdapter::makeSuttonChen(RealType c, RealType m, RealType n, 
                                         RealType alpha, RealType epsilon){ 

    if (isSuttonChen()){
      at_->removeProperty(SCtypeID);
    }

    SCAtypeParameters* scParam = new SCAtypeParameters();
    scParam->c = c;
    scParam->m = m;
    scParam->n = n;
    scParam->alpha = alpha;
    scParam->epsilon = epsilon;
    
    at_->addProperty(new SCAtypeData(SCtypeID, scParam));
  }
}
