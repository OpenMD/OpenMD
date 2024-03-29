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

#include "optimization/Armijo.hpp"

#include "optimization/Method.hpp"
#include "optimization/Problem.hpp"

namespace QuantLib {

  RealType ArmijoLineSearch::operator()(Problem& P, EndCriteria::Type& ecType,
                                        const EndCriteria& endCriteria,
                                        const RealType t_ini) {
    // OptimizationMethod& method = P.method();
    Constraint& constraint = P.constraint();
    succeed_               = true;
    bool maxIter           = false;
    RealType qtold, t = t_ini;
    size_t loopNumber = 0;

    RealType q0  = P.functionValue();
    RealType qp0 = P.gradientNormValue();

    qt_ = q0;
    qpt_ =
        (gradient_.empty()) ? qp0 : -P.DotProduct(gradient_, searchDirection_);

    // Initialize gradient
    gradient_ = DynamicVector<RealType>(P.currentValue().size());
    // Compute new point
    xtd_ = P.currentValue();
    t    = update(xtd_, searchDirection_, t, constraint);
    // Compute function value at the new point
    qt_ = P.value(xtd_);

    // Enter in the loop if the criterion is not satisfied
    if ((qt_ - q0) > -alpha_ * t * qpt_) {
      do {
        loopNumber++;
        // Decrease step
        t *= beta_;
        // Store old value of the function
        qtold = qt_;
        // New point value
        xtd_ = P.currentValue();
        t    = update(xtd_, searchDirection_, t, constraint);

        // Compute function value at the new point
        qt_ = P.value(xtd_);
        P.gradient(gradient_, xtd_);
        // and it squared norm
        maxIter = endCriteria.checkMaxIterations(loopNumber, ecType);
      } while ((((qt_ - q0) > (-alpha_ * t * qpt_)) ||
                ((qtold - q0) <= (-alpha_ * t * qpt_ / beta_))) &&
               (!maxIter));
    }

    if (maxIter) succeed_ = false;

    // Compute new gradient
    P.gradient(gradient_, xtd_);
    // and it squared norm
    // qpt_ = P.computeGradientNormValue(gradient_);
    qpt_ = P.DotProduct(gradient_, gradient_);

    // Return new step value
    return t;
  }

}  // namespace QuantLib
