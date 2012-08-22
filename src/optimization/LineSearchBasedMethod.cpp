/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006 Ferdinando Ametrano
 Copyright (C) 2009 Frédéric Degraeve

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

#include "optimization/LineSearchBasedMethod.hpp"
#include "optimization/Problem.hpp"
#include "optimization/LineSearch.hpp"
#include "optimization/Armijo.hpp"
#include "utils/NumericConstant.hpp"

using namespace OpenMD;
namespace QuantLib {

    LineSearchBasedMethod::LineSearchBasedMethod(LineSearch* lineSearch)
    : lineSearch_(lineSearch) {
        if (!lineSearch_)
            lineSearch_ = new ArmijoLineSearch();
    }

    EndCriteria::Type
    LineSearchBasedMethod::minimize(Problem& P,
                                    const EndCriteria& endCriteria) {
        // Initializations
        RealType ftol = endCriteria.functionEpsilon();
        size_t maxStationaryStateIterations_
            = endCriteria.maxStationaryStateIterations();
        EndCriteria::Type ecType = EndCriteria::None;   // reset end criteria
        P.reset();                                      // reset problem
        DynamicVector<RealType> x_ = P.currentValue(); // store the starting point
        size_t iterationNumber_ = 0;
        // dimension line search
        lineSearch_->searchDirection() = DynamicVector<RealType>(x_.size());
        bool done = false;

        // function and squared norm of gradient values;
        RealType fnew, fold, gold2;
        RealType fdiff;
        // classical initial value for line-search step
        RealType t = 1.0;
        // Set gradient g at the size of the optimization problem
        // search direction
        size_t sz = lineSearch_->searchDirection().size();
        DynamicVector<RealType> prevGradient(sz), d(sz), sddiff(sz), direction(sz);
        // Initialize objective function, gradient prevGradient and
        // search direction
        P.setFunctionValue(P.valueAndGradient(prevGradient, x_));
        P.setGradientNormValue(P.DotProduct(prevGradient, prevGradient));
        lineSearch_->searchDirection() = -prevGradient;

        bool first_time = true;
        // Loop over iterations
        do {
            // Linesearch
            if (!first_time)
                prevGradient = lineSearch_->lastGradient();
            t = (*lineSearch_)(P, ecType, endCriteria, t);
            // don't throw: it can fail just because maxIterations exceeded
            //QL_REQUIRE(lineSearch_->succeed(), "line-search failed!");
            if (lineSearch_->succeed())
            {
                // Updates

                // New point
                x_ = lineSearch_->lastX();
                // New function value
                fold = P.functionValue();
                P.setFunctionValue(lineSearch_->lastFunctionValue());
                // New gradient and search direction vectors

                // orthogonalization coef
                gold2 = P.gradientNormValue();
                P.setGradientNormValue(lineSearch_->lastGradientNorm2());

                // conjugate gradient search direction
                direction = getUpdatedDirection(P, gold2, prevGradient);
                sddiff = direction - lineSearch_->searchDirection();
                lineSearch_->searchDirection() = direction;
                // Now compute accuracy and check end criteria
                // Numerical Recipes exit strategy on fx (see NR in C++, p.423)
                fnew = P.functionValue();
                fdiff = 2.0*std::fabs(fnew-fold) /
                    (std::fabs(fnew) + std::fabs(fold) + NumericConstant::epsilon);
                if (fdiff < ftol ||
                    endCriteria.checkMaxIterations(iterationNumber_, ecType)) {
                    endCriteria.checkStationaryFunctionValue(0.0, 0.0,
                        maxStationaryStateIterations_, ecType);
                    endCriteria.checkMaxIterations(iterationNumber_, ecType);
                    return ecType;
                }
                P.setCurrentValue(x_);      // update problem current value
                ++iterationNumber_;         // Increase iteration number
                first_time = false;
            } else {
                done = true;
            }
        } while (!done);
        P.setCurrentValue(x_);
        return ecType;
    }

}
