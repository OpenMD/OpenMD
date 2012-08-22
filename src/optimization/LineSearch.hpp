/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
  Copyright (C) 2006 Ferdinando Ametrano
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

/*! \file linesearch.hpp
  \brief Line search abstract class
*/

#ifndef quantlib_optimization_line_search_h_
#define quantlib_optimization_line_search_h_

#include "math/DynamicVector.hpp"
#include "optimization/EndCriteria.hpp"

using namespace OpenMD;
namespace QuantLib {

    class Problem;
    class Constraint;
    class EndCriteria;

    //! Base class for line search
    class LineSearch {
    public:
        //! Default constructor
        LineSearch(RealType = 0.0)
            : qt_(0.0), qpt_(0.0), succeed_(true) {}
        //! Destructor
        virtual ~LineSearch() {}

        //! return last x value
        const DynamicVector<RealType>& lastX() { return xtd_; }
        //! return last objective function value
        RealType lastFunctionValue() { return qt_; }
        //! return last gradient
        const DynamicVector<RealType>& lastGradient() { return gradient_; }
        //! return square norm of last gradient
        RealType lastGradientNorm2() { return qpt_;}

        bool succeed() { return succeed_; }

        //! Perform line search
        virtual RealType operator()(Problem& P, // Optimization problem
                                    EndCriteria::Type& ecType,
                                    const EndCriteria&,
                                    const RealType t_ini) = 0;  // initial value of line-search step
        RealType update(DynamicVector<RealType>& params,
                        const DynamicVector<RealType>& direction,
                        RealType beta,
                        const Constraint& constraint);

        //! current value of the search direction
        const DynamicVector<RealType>& searchDirection() const { return searchDirection_; }
        DynamicVector<RealType>& searchDirection() { return searchDirection_; }
    protected:
        //! current values of the search direction
        DynamicVector<RealType> searchDirection_;
        //! new x and its gradient
        DynamicVector<RealType> xtd_, gradient_;
        //! objective function value and gradient norm corresponding to xtd_
        RealType qt_, qpt_;
        //! flag to know if linesearch succeed
        bool succeed_;
        
    };
}

#endif
