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
 
#include "types/StickyAdapter.hpp"
#include "utils/simError.h"
#include <cstdio>
#include <memory>

namespace OpenMD {
  
  bool StickyAdapter::isSticky() {
    return at_->hasProperty(StickyTypeID);
  }
  
  StickyAtypeParameters StickyAdapter::getStickyParam() {
    
    if (!isSticky()) {
      sprintf( painCave.errMsg,               
               "StickyAdapter::getStickyParam was passed an atomType (%s)\n"
               "\tthat does not appear to be a Sticky atom.\n",
               at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
    
    std::shared_ptr<GenericData> data = at_->getPropertyByName(StickyTypeID);
    if (data == nullptr) {
      sprintf( painCave.errMsg, 
               "StickyAdapter::getStickyParam could not find Sticky\n"
               "\tparameters for atomType %s.\n", at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError(); 
    }
    
    std::shared_ptr<StickyAtypeData> stickyData = std::dynamic_pointer_cast<StickyAtypeData>(data);
    if (stickyData == nullptr) {
      sprintf( painCave.errMsg,
               "StickyAdapter::getStickyParam could not convert\n"
               "\tGenericData to StickyAtypeData for atom type %s\n", 
               at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();          
    }
    
    return stickyData->getData();
  }

  bool StickyAdapter::isStickyPower() {
    if (isSticky()) {
      StickyAtypeParameters stickyParam = getStickyParam();
      return stickyParam.isPower;
    } else 
      return false;
  }

  RealType StickyAdapter::getW0() {    
    StickyAtypeParameters stickyParam = getStickyParam();
    return stickyParam.w0;
  }

  RealType StickyAdapter::getV0() {    
    StickyAtypeParameters stickyParam = getStickyParam();
    return stickyParam.v0;
  }

  RealType StickyAdapter::getV0p() {    
    StickyAtypeParameters stickyParam = getStickyParam();
    return stickyParam.v0p;
  }

  RealType StickyAdapter::getRl() {    
    StickyAtypeParameters stickyParam = getStickyParam();
    return stickyParam.rl;
  }

  RealType StickyAdapter::getRu() {    
    StickyAtypeParameters stickyParam = getStickyParam();
    return stickyParam.ru;
  }

  RealType StickyAdapter::getRlp() {    
    StickyAtypeParameters stickyParam = getStickyParam();
    return stickyParam.rlp;
  }

  RealType StickyAdapter::getRup() {    
    StickyAtypeParameters stickyParam = getStickyParam();
    return stickyParam.rup;
  }

    
  void StickyAdapter::makeSticky(RealType w0, RealType v0, RealType v0p, 
                                 RealType rl, RealType ru, RealType rlp, 
                                 RealType rup, bool isPower){ 

    if (isSticky()){
      at_->removeProperty(StickyTypeID);
    }

    StickyAtypeParameters stickyParam {};
    stickyParam.w0 = w0;
    stickyParam.v0 = v0;
    stickyParam.v0p = v0p;
    stickyParam.rl = rl;
    stickyParam.ru = ru;
    stickyParam.rlp = rlp;
    stickyParam.rup = rup;
    stickyParam.isPower = isPower;
    
    at_->addProperty(std::make_shared<StickyAtypeData>(StickyTypeID, stickyParam));
  }
}
