/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2001, 2002, 2003 Sadruddin Rejeb

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

#include "optimization/Constraint.hpp"
#include "utils/simError.h"
#include <cstdio>

namespace QuantLib {

    Constraint::Constraint(Constraint::Impl* impl)
    : impl_(impl) {}

    RealType Constraint::update(DynamicVector<RealType>& params,
                                const DynamicVector<RealType>& direction,
                                RealType beta) {

        RealType diff=beta;
        DynamicVector<RealType> newParams = params + diff*direction;
        bool valid = test(newParams);
        int icount = 0;
        while (!valid) {
            if (icount > 200) {
                sprintf(painCave.errMsg,  "can't update parameter vector\n");
                painCave.isFatal = 1;
                painCave.severity = OPENMD_ERROR;
                simError();
            }
            diff *= 0.5;
            icount ++;
            newParams = params + diff*direction;
            valid = test(newParams);
        }
        params += diff*direction;
        return diff;
    }

}
