/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
  Copyright (C) 2007 Ferdinando Ametrano
  Copyright (C) 2007 François du Vignaud
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

/*! \file problem.hpp
  \brief Abstract optimization problem class
*/

#ifndef quantlib_optimization_problem_h
#define quantlib_optimization_problem_h

#include "optimization/Method.hpp"
#include "optimization/ObjectiveFunction.hpp"
#include "optimization/StatusFunction.hpp"

namespace QuantLib {
    
    class Constraint;
    //! Constrained optimization problem
    class Problem {
    public:
        //! default constructor
        Problem(ObjectiveFunction& objectiveFunction,
                Constraint& constraint,
                OpenMD::StatusFunction& statFunc,
                const DynamicVector<RealType>& initialValue = DynamicVector<RealType>())
            : objectiveFunction_(objectiveFunction), constraint_(constraint),
              currentValue_(initialValue), statusFunction_(statFunc) {}
        
        /*! \warning it does not reset the current minumum to any initial value
         */
        void reset();
        
        //! call objective function computation and increment evaluation counter
        RealType value(const DynamicVector<RealType>& x);
        
        //! call objective function gradient computation and increment
        //  evaluation counter
        void gradient(DynamicVector<RealType>& grad_f,
                      const DynamicVector<RealType>& x);

        //! call objective function computation and it gradient
        RealType valueAndGradient(DynamicVector<RealType>& grad_f,
                                  const DynamicVector<RealType>& x);

        //! Constraint
        Constraint& constraint() const { return constraint_; }

        //! Objective function
        ObjectiveFunction& objectiveFunction() const { return objectiveFunction_; }

        void setCurrentValue(const DynamicVector<RealType>& currentValue) {
            currentValue_=currentValue;
            statusFunction_.writeStatus(currentValue);
        }

        //! current value of the local minimum
        const DynamicVector<RealType>& currentValue() const { return currentValue_; }

        void setFunctionValue(RealType functionValue) {
            functionValue_=functionValue;
        }

        //! value of objective function
        RealType functionValue() const { return functionValue_; }

        void setGradientNormValue(RealType squaredNorm) {
            squaredNorm_=squaredNorm;
        }
        //! value of objective function gradient norm
        RealType gradientNormValue() const { return squaredNorm_; }

        //! number of evaluation of objective function
        int functionEvaluation() const { return functionEvaluation_; }

        //! number of evaluation of objective function gradient
        int gradientEvaluation() const { return gradientEvaluation_; }

        RealType DotProduct(DynamicVector<RealType>& v1, DynamicVector<RealType>& v2);
        RealType computeGradientNormValue(DynamicVector<RealType>& grad_f);
        

    protected:
        //! Unconstrained objective function
        ObjectiveFunction& objectiveFunction_;
        //! Constraint
        Constraint& constraint_;
        //! current value of the local minimum
        DynamicVector<RealType> currentValue_;
        //! function and gradient norm values at the curentValue_ (i.e. the last step)
        RealType functionValue_, squaredNorm_;
        //! number of evaluation of objective function and its gradient
        int functionEvaluation_, gradientEvaluation_;
        //! status function
        StatusFunction& statusFunction_;

    };
    
    // inline definitions
    inline RealType Problem::value(const DynamicVector<RealType>& x) {
        ++functionEvaluation_;
        return objectiveFunction_.value(x);
    }
    
    inline void Problem::gradient(DynamicVector<RealType>& grad_f,
                                  const DynamicVector<RealType>& x) {
        ++gradientEvaluation_;
        objectiveFunction_.gradient(grad_f, x);
    }
    
    inline RealType Problem::valueAndGradient(DynamicVector<RealType>& grad_f,
                                              const DynamicVector<RealType>& x) {
        ++functionEvaluation_;
        ++gradientEvaluation_;
        return objectiveFunction_.valueAndGradient(grad_f, x);
    }
    
    inline void Problem::reset() {
        functionEvaluation_ = gradientEvaluation_ = 0;
        functionValue_ = squaredNorm_ = 0;
    }
    
}

#endif
