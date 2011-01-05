/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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
 * [4]  Vardeman & Gezelter, in progress (2009).                        
 */

#include <cmath>
#include <iostream>

#include "nonbonded/SwitchingFunction.hpp"
#include "utils/simError.h"

using namespace std;
namespace OpenMD {

  SwitchingFunction::SwitchingFunction() : np_(150), haveSpline_(false), 
                                           isCubic_(true), functionType_(cubic) {
    switchSpline_ = new CubicSpline();
  }

  void SwitchingFunction::setSwitchType(SwitchingFunctionType sft) {
    if ((sft == fifth_order_poly) || (sft == cubic)) {
      if (haveSpline_) {
        delete switchSpline_;
        switchSpline_ = new CubicSpline();
        setSwitch(rin_, rout_);
      } 
    } else {
      sprintf( painCave.errMsg,
               "SwitchingFunction::setSwitchType was given unknown function type\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();          
    }
  }

  void SwitchingFunction::setSwitch(RealType rinner, RealType router) {
    if (router < rinner) {
          sprintf( painCave.errMsg,
                   "SwitchingFunction::setSwitch was given rinner (%lf) which was\n" 
                   "\tlarger than router (%lf).\n", rinner, router);
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();          
    }
    if (router < 0.0) {
          sprintf( painCave.errMsg,
                   "SwitchingFunction::setSwitch was given router (%lf) which was\n" 
                   "\tless than zero.\n", router);
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();           
    }   
    if (rinner < 0.0) {
          sprintf( painCave.errMsg,
                   "SwitchingFunction::setSwitch was given rinner (%lf) which was\n"
                   "\tless than zero.\n", router);
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();           
    }

    rin_ = rinner;
    rout_ = router;
    rin2_ = rin_ * rin_;
    rout2_ = rout_ * rout_;

    if ((router - rinner) < 1.0e-8) {
      // no reason to set up spline if the switching region is tiny
      return;
    }
    
    // setup r -> sw lookup spline
    if (functionType_ == fifth_order_poly) {
      isCubic_ = false;
      RealType c0 = 1.0;
      RealType c3 = -10.0;
      RealType c4 = 15.0;
      RealType c5 = -6.0;
 
      RealType dx, r, yval, rval, rval2, rval3, rval4, rval5;
      RealType rvaldi, rvaldi2, rvaldi3, rvaldi4, rvaldi5;

      dx = (rout_ - rin_) / RealType(np_-1);

      for (int i = 0; i < np_; i++) {
        r = rin_ + RealType(i)*dx;
        rval = ( r - rin_ );
        rval2 = rval*rval;
        rval3 = rval2*rval;
        rval4 = rval2*rval2;
        rval5 = rval3*rval2;
        rvaldi = 1.0 / ( rout_ - rin_ );
        rvaldi2 = rvaldi*rvaldi;
        rvaldi3 = rvaldi2*rvaldi;
        rvaldi4 = rvaldi2*rvaldi2;
        rvaldi5 = rvaldi3*rvaldi2;
        yval= c0 + c3*rval3*rvaldi3 + c4*rval4*rvaldi4 + c5*rval5*rvaldi5;
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

  bool SwitchingFunction::getSwitch(const RealType &r2, RealType &sw, RealType &dswdr, 
                                    RealType &r) {

    sw = 1.0;
    dswdr = 0.0;

    bool in_switching_region = false;

    if (r2 > rin2_) {
      if (r2 > rout2_) {
        sw = 0.0;
      } else {
        in_switching_region = true;
        r = sqrt(r2);
        pair<RealType, RealType> result = switchSpline_->getValueAndDerivativeAt(r);
        sw = result.first;
        dswdr = result.second;
      }
    }
    return in_switching_region;
  }
}
