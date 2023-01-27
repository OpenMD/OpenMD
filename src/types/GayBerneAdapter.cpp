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

#include "types/GayBerneAdapter.hpp"

#include <cstdio>
#include <memory>

#include "utils/simError.h"

namespace OpenMD {

  bool GayBerneAdapter::isGayBerne() { return at_->hasProperty(GBtypeID); }

  GBAtypeParameters GayBerneAdapter::getGayBerneParam() {
    if (!isGayBerne()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "GayBerneAdapter::getGayBerneParam was passed an atomType (%s)\n"
               "\tthat does not appear to be a Gay-Berne atom.\n",
               at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    std::shared_ptr<GenericData> data = at_->getPropertyByName(GBtypeID);
    if (data == nullptr) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "GayBerneAdapter::getGayBerneParam could not find Gay-Berne\n"
               "\tparameters for atomType %s.\n",
               at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    std::shared_ptr<GBAtypeData> gbData =
        std::dynamic_pointer_cast<GBAtypeData>(data);
    if (gbData == nullptr) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "GayBerneAdapter::getGayBerneParam could not convert\n"
               "\tGenericData to GBAtypeData for atom type %s\n",
               at_->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    return gbData->getData();
  }

  RealType GayBerneAdapter::getD() {
    GBAtypeParameters gbParam = getGayBerneParam();
    return gbParam.GB_d;
  }

  RealType GayBerneAdapter::getL() {
    GBAtypeParameters gbParam = getGayBerneParam();
    return gbParam.GB_l;
  }

  RealType GayBerneAdapter::getEpsX() {
    GBAtypeParameters gbParam = getGayBerneParam();
    return gbParam.GB_eps_X;
  }

  RealType GayBerneAdapter::getEpsS() {
    GBAtypeParameters gbParam = getGayBerneParam();
    return gbParam.GB_eps_S;
  }

  RealType GayBerneAdapter::getEpsE() {
    GBAtypeParameters gbParam = getGayBerneParam();
    return gbParam.GB_eps_E;
  }

  RealType GayBerneAdapter::getDw() {
    GBAtypeParameters gbParam = getGayBerneParam();
    return gbParam.GB_dw;
  }

  void GayBerneAdapter::makeGayBerne(RealType d, RealType l, RealType eps_X,
                                     RealType eps_S, RealType eps_E,
                                     RealType dw) {
    if (isGayBerne()) { at_->removeProperty(GBtypeID); }

    GBAtypeParameters gbParam {};
    gbParam.GB_d     = d;
    gbParam.GB_l     = l;
    gbParam.GB_eps_X = eps_X;
    gbParam.GB_eps_S = eps_S;
    gbParam.GB_eps_E = eps_E;
    gbParam.GB_dw    = dw;

    at_->addProperty(std::make_shared<GBAtypeData>(GBtypeID, gbParam));
  }
}  // namespace OpenMD
