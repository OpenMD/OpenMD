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

/*! \file linesearchbasedmethod.hpp
    \brief Abstract optimization method class
*/

#ifndef quantlib_line_search_based_optimization_method_h
#define quantlib_line_search_based_optimization_method_h

#include "optimization/Method.hpp"
#include "math/DynamicVector.hpp"

using namespace OpenMD;
namespace QuantLib {

    class LineSearch;

    class LineSearchBasedMethod : public OptimizationMethod {
      public:
        LineSearchBasedMethod(LineSearch* lSearch = NULL);
        virtual ~LineSearchBasedMethod() {}

        virtual EndCriteria::Type minimize(Problem& P,
                                           const EndCriteria& endCriteria);
      protected:
        //! computes the new search direction
        virtual DynamicVector<RealType>
        getUpdatedDirection(const Problem &P,
                            RealType gold2,
                            const DynamicVector<RealType>& gradient) = 0;
        //! line search
       LineSearch* lineSearch_;
    };

}

#endif
