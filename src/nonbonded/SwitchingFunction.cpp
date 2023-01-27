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

#include "nonbonded/SwitchingFunction.hpp"

#include <cmath>
#include <cstdio>
#include <iostream>

#include "utils/simError.h"

using namespace std;
namespace OpenMD {

  SwitchingFunction::SwitchingFunction() :
      functionType_(cubic), haveSpline_(false), isCubic_(true), np_(150) {
    switchSpline_ = std::make_shared<CubicSpline>();
  }

  void SwitchingFunction::setSwitchType(SwitchingFunctionType sft) {
    if ((sft == fifth_order_poly) || (sft == cubic)) {
      if (haveSpline_) {
        switchSpline_ = std::make_shared<CubicSpline>();
        setSwitch(rin_, rout_);
      }
    } else {
      snprintf(
          painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
          "SwitchingFunction::setSwitchType was given unknown function type\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
  }

  void SwitchingFunction::setSwitch(RealType rinner, RealType router) {
    if (router < rinner) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "SwitchingFunction::setSwitch was given rinner (%lf) which was\n"
               "\tlarger than router (%lf).\n",
               rinner, router);
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
    if (router < 0.0) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "SwitchingFunction::setSwitch was given router (%lf) which was\n"
               "\tless than zero.\n",
               router);
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
    if (rinner < 0.0) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "SwitchingFunction::setSwitch was given rinner (%lf) which was\n"
               "\tless than zero.\n",
               router);
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    rin_   = rinner;
    rout_  = router;
    rin2_  = rin_ * rin_;
    rout2_ = rout_ * rout_;

    if ((router - rinner) < 1.0e-8) {
      // no reason to set up spline if the switching region is tiny
      return;
    }

    // setup r -> sw lookup spline
    if (functionType_ == fifth_order_poly) {
      isCubic_    = false;
      RealType c0 = 1.0;
      RealType c3 = -10.0;
      RealType c4 = 15.0;
      RealType c5 = -6.0;

      RealType dx, r, yval, rval, rval2, rval3, rval4, rval5;
      RealType rvaldi, rvaldi2, rvaldi3, rvaldi4, rvaldi5;

      dx = (rout_ - rin_) / RealType(np_ - 1);

      for (int i = 0; i < np_; i++) {
        r       = rin_ + RealType(i) * dx;
        rval    = (r - rin_);
        rval2   = rval * rval;
        rval3   = rval2 * rval;
        rval4   = rval2 * rval2;
        rval5   = rval3 * rval2;
        rvaldi  = 1.0 / (rout_ - rin_);
        rvaldi2 = rvaldi * rvaldi;
        rvaldi3 = rvaldi2 * rvaldi;
        rvaldi4 = rvaldi2 * rvaldi2;
        rvaldi5 = rvaldi3 * rvaldi2;
        yval    = c0 + c3 * rval3 * rvaldi3 + c4 * rval4 * rvaldi4 +
               c5 * rval5 * rvaldi5;
        switchSpline_->addPoint(r, yval);
      }
    } else {
      // cubic splines only need 2 points to do a cubic switching function...
      isCubic_ = true;
      switchSpline_->addPoint(rin_, 1.0);
      switchSpline_->addPoint(rout_, 0.0);
    }
    haveSpline_ = true;
    return;
  }

  bool SwitchingFunction::getSwitch(const RealType& r2, RealType& sw,
                                    RealType& dswdr, RealType& r) {
    sw    = 1.0;
    dswdr = 0.0;

    bool in_switching_region = false;

    if (r2 > rin2_) {
      if (r2 > rout2_) {
        sw = 0.0;
      } else {
        in_switching_region = true;
        r                   = sqrt(r2);
        switchSpline_->getValueAndDerivativeAt(r, sw, dswdr);
      }
    }
    return in_switching_region;
  }
}  // namespace OpenMD
