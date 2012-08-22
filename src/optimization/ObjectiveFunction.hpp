/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
  Copyright (C) 2001, 2002, 2003 Nicolas Di Césaré

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

/*! \file ObjectiveFunction.hpp
  \brief Optimization objective function class
*/

#ifndef optimization_objectivefunction_h
#define optimization_objectivefunction_h
#include "config.h"
#include "math/DynamicVector.hpp"

using namespace OpenMD;
namespace QuantLib {
    
    //!  Objective function abstract class for optimization problem
    class ObjectiveFunction {
    public:
        virtual ~ObjectiveFunction() {}
        //! method to overload to compute the objective function value in x
        virtual RealType value(const DynamicVector<RealType>& x)  = 0;
        
        //! method to overload to compute grad_f, the first derivative of
        //  the objective function with respect to x
        virtual void gradient(DynamicVector<RealType>& grad, const DynamicVector<RealType>& x) {
            RealType eps = finiteDifferenceEpsilon(), fp, fm;
            DynamicVector<RealType> xx(x);
            for (size_t i=0; i<x.size(); i++) {
                xx[i] += eps;
                fp = value(xx);
                xx[i] -= 2.0*eps;
                fm = value(xx);
                grad[i] = 0.5*(fp - fm)/eps;
                xx[i] = x[i];
            }
        }

        //! method to overload to compute grad_f, the first derivative
        //  of the objective function with respect to x and also the
        //  objective function
        virtual RealType valueAndGradient(DynamicVector<RealType>& grad,
                                          const DynamicVector<RealType>& x) {
            gradient(grad, x);
            return value(x);
        }

        //! Default epsilon for finite difference method :
        virtual RealType finiteDifferenceEpsilon() const { return 1e-8; }
    };

    class ParametersTransformation {
    public:
        virtual ~ParametersTransformation() {}
        virtual DynamicVector<RealType> direct(const DynamicVector<RealType>& x) const = 0;
        virtual DynamicVector<RealType> inverse(const DynamicVector<RealType>& x) const = 0;
    };
}

#endif
