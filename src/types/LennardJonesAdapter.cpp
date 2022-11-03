/*
 * Copyright (c) 2004-2021 The University of Notre Dame. All Rights Reserved.
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

#include "types/LennardJonesAdapter.hpp"

#include <cstdio>
#include <memory>

#include "utils/simError.h"

namespace OpenMD {

  bool LennardJonesAdapter::isLennardJones() {
    return at_->hasProperty(LJtypeID);
  }

  LJAtypeParameters LennardJonesAdapter::getLJParam() {
    if (!isLennardJones()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
              "LennardJonesAdapter::getLJParam was passed an atomType (%s)\n"
              "\tthat does not appear to be a Lennard-Jones atom.\n",
              at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    std::shared_ptr<GenericData> data = at_->getPropertyByName(LJtypeID);
    if (data == nullptr) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
              "LennardJonesAdapter::getLJParam could not find Lennard-Jones\n"
              "\tparameters for atomType %s.\n",
              at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    std::shared_ptr<LJAtypeData> ljData =
        std::dynamic_pointer_cast<LJAtypeData>(data);
    if (ljData == nullptr) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
              "LennardJonesAdapter::getLJParam could not convert\n"
              "\tGenericData to LJAtypeData for atom type %s\n",
              at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    return ljData->getData();
  }

  RealType LennardJonesAdapter::getSigma() {
    LJAtypeParameters ljParam = getLJParam();
    return ljParam.sigma;
  }

  RealType LennardJonesAdapter::getEpsilon() {
    LJAtypeParameters ljParam = getLJParam();
    return ljParam.epsilon;
  }

  bool LennardJonesAdapter::isSoft() {
    LJAtypeParameters ljParam = getLJParam();
    return ljParam.isSoft;
  }

  void LennardJonesAdapter::makeLennardJones(RealType sigma, RealType epsilon,
                                             bool isSoft) {
    if (isLennardJones()) { at_->removeProperty(LJtypeID); }

    LJAtypeParameters ljParam {};
    ljParam.epsilon = epsilon;
    ljParam.sigma   = sigma;
    ljParam.isSoft  = isSoft;

    at_->addProperty(std::make_shared<LJAtypeData>(LJtypeID, ljParam));
  }
}  // namespace OpenMD
