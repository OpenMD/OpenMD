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
 
#include "types/GayBerneAdapter.hpp"
#include "utils/simError.h"
#include <cstdio>

namespace OpenMD {

  bool GayBerneAdapter::isGayBerne() {
    return at_->hasProperty(GBtypeID);
  }
  
  GBAtypeParameters* GayBerneAdapter::getGayBerneParam() {
    
    if (!isGayBerne()) {
      sprintf( painCave.errMsg,               
               "GayBerneAdapter::getGayBerneParam was passed an atomType (%s)\n"
               "\tthat does not appear to be a Gay-Berne atom.\n",
               at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
    
    GenericData* data = at_->getPropertyByName(GBtypeID);
    if (data == NULL) {
      sprintf( painCave.errMsg, 
               "GayBerneAdapter::getGayBerneParam could not find Gay-Berne\n"
               "\tparameters for atomType %s.\n", at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError(); 
    }
    
    GBAtypeData* gbData = dynamic_cast<GBAtypeData*>(data);
    if (gbData == NULL) {
      sprintf( painCave.errMsg,
               "GayBerneAdapter::getGayBerneParam could not convert\n"
               "\tGenericData to GBAtypeData for atom type %s\n", 
               at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();          
    }
    
    return gbData->getData();
  }
  
  RealType GayBerneAdapter::getD() {    
    GBAtypeParameters* gbParam = getGayBerneParam();
    return gbParam->GB_d;
  }

  RealType GayBerneAdapter::getL() {    
    GBAtypeParameters* gbParam = getGayBerneParam();
    return gbParam->GB_l;
  }

  RealType GayBerneAdapter::getEpsX() {    
    GBAtypeParameters* gbParam = getGayBerneParam();
    return gbParam->GB_eps_X;
  }

  RealType GayBerneAdapter::getEpsS() {    
    GBAtypeParameters* gbParam = getGayBerneParam();
    return gbParam->GB_eps_S;
  }

  RealType GayBerneAdapter::getEpsE() {    
    GBAtypeParameters* gbParam = getGayBerneParam();
    return gbParam->GB_eps_E;
  }

  RealType GayBerneAdapter::getDw() {    
    GBAtypeParameters* gbParam = getGayBerneParam();
    return gbParam->GB_dw;
  }
    
  void GayBerneAdapter::makeGayBerne(RealType d, RealType l, RealType eps_X, 
                                     RealType eps_S, RealType eps_E, 
                                     RealType dw){
    if (isGayBerne()){
      at_->removeProperty(GBtypeID);
    }

    GBAtypeParameters* gbParam = new GBAtypeParameters();
    gbParam->GB_d = d;
    gbParam->GB_l = l;
    gbParam->GB_eps_X = eps_X;
    gbParam->GB_eps_S = eps_S;
    gbParam->GB_eps_E = eps_E;
    gbParam->GB_dw = dw;
    
    at_->addProperty(new GBAtypeData(GBtypeID, gbParam));
  }
}
