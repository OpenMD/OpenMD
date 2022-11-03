/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2001, 2002, 2003 Nicolas Di C�sar�

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#include "optimization/LineSearch.hpp"

#include <cstdio>

#include "optimization/Constraint.hpp"
#include "optimization/Problem.hpp"
#include "utils/simError.h"

namespace QuantLib {

  RealType LineSearch::update(DynamicVector<RealType>& params,
                              const DynamicVector<RealType>& direction,
                              RealType beta, const Constraint& constraint) {
    RealType diff                     = beta;
    DynamicVector<RealType> newParams = params + diff * direction;
    bool valid                        = constraint.test(newParams);
    int icount                        = 0;
    while (!valid) {
      if (icount > 200) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "can't update linesearch\n");
          painCave.isFatal  = 1;
          painCave.severity = OPENMD_ERROR;
          simError();
      }
      diff *= 0.5;
      icount++;
      newParams = params + diff * direction;
      valid     = constraint.test(newParams);
    }
    params += diff * direction;
    return diff;
  }

}  // namespace QuantLib
