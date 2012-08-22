/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006, 2007 Ferdinando Ametrano
 Copyright (C) 2007 Marco Bianchetti
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

/*! \file endcriteria.hpp
    \brief Optimization criteria class
*/

#ifndef quantlib_optimization_criteria_hpp
#define quantlib_optimization_criteria_hpp

#include "config.h"
#include <ostream>

namespace QuantLib {

    //! Criteria to end optimization process:
    /*! - maximum number of iterations AND minimum number of iterations around stationary point
        - x (independent variable) stationary point
        - y=f(x) (dependent variable) stationary point
        - stationary gradient
    */
    class EndCriteria {
      public:
        enum Type {None,
                   MaxIterations,
                   StationaryPoint,
                   StationaryFunctionValue,
                   StationaryFunctionAccuracy,
                   ZeroGradientNorm,
                   Unknown};

        //! Initialization constructor
        EndCriteria(size_t maxIterations,
                    size_t maxStationaryStateIterations,
                    RealType rootEpsilon,
                    RealType functionEpsilon,
                    RealType gradientNormEpsilon);

        // Inspectors
        size_t maxIterations() const;
        size_t maxStationaryStateIterations() const;
        RealType rootEpsilon() const;
        RealType functionEpsilon() const;
        RealType gradientNormEpsilon() const;

        /*! Test if the number of iterations is not too big 
            and if a minimum point is not reached */
        bool operator()(const size_t iteration,
                        size_t& statState,
                        const bool positiveOptimization,
                        const RealType fold,
                        const RealType normgold,
                        const RealType fnew,
                        const RealType normgnew,
                        EndCriteria::Type& ecType) const;

        /*! Test if the number of iteration is below MaxIterations */
        bool checkMaxIterations(const size_t iteration,
                                EndCriteria::Type& ecType) const;
        /*! Test if the root variation is below rootEpsilon */
        bool checkStationaryPoint(const RealType xOld,
                                  const RealType xNew,
                                  size_t& statStateIterations,
                                  EndCriteria::Type& ecType) const;
        /*! Test if the function variation is below functionEpsilon */
        bool checkStationaryFunctionValue(const RealType fxOld,
                                          const RealType fxNew,
                                          size_t& statStateIterations,
                                          EndCriteria::Type& ecType) const;
        /*! Test if the function value is below functionEpsilon */
        bool checkStationaryFunctionAccuracy(const RealType f,
                                             const bool positiveOptimization,
                                             EndCriteria::Type& ecType) const;
        /*! Test if the gradient norm value is below gradientNormEpsilon */
        bool checkZeroGradientNorm(const RealType gNorm,
                                   EndCriteria::Type& ecType) const;
        
    protected:
        //! Maximum number of iterations
        size_t maxIterations_;
        //! Maximun number of iterations in stationary state
        mutable size_t maxStationaryStateIterations_;
        //! root, function and gradient epsilons
        RealType rootEpsilon_, functionEpsilon_, gradientNormEpsilon_;

    };

    std::ostream& operator<<(std::ostream& out, EndCriteria::Type ecType);

}

#endif
