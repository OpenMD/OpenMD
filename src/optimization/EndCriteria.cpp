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

#include "optimization/EndCriteria.hpp"
#include "utils/simError.h"
#include <cmath>
#include <cstdio>


namespace QuantLib {
    
    EndCriteria::EndCriteria(size_t maxIterations,
                             size_t maxStationaryStateIterations,
                             RealType rootEpsilon,
                             RealType functionEpsilon,
                             RealType gradientNormEpsilon)
        : maxIterations_(maxIterations),
          maxStationaryStateIterations_(maxStationaryStateIterations),
          rootEpsilon_(rootEpsilon),
          functionEpsilon_(functionEpsilon),
          gradientNormEpsilon_(gradientNormEpsilon) {
        
        
        // replaced the QL_REQUIRE macro with OpenMD's simError calls
        if (maxStationaryStateIterations_ <= 1) {
            sprintf(painCave.errMsg,
                    "maxStationaryStateIterations_ ( %lu ) "
                    "must be greater than one\n",
                    (unsigned long)maxStationaryStateIterations_);
            painCave.isFatal = 1;
            painCave.severity = OPENMD_ERROR;
            simError();
        }
        if (maxStationaryStateIterations_ > maxIterations_) {
            sprintf(painCave.errMsg,
                    "maxStationaryStateIterations_ ( %lu ) "
                    "must be less than maxIterations_ ( %lu )\n",
                    (unsigned long)maxStationaryStateIterations_, 
                    (unsigned long)maxIterations_);
            painCave.isFatal = 1;
            painCave.severity = OPENMD_ERROR;
            simError();
        }
        
    }
    
    bool EndCriteria::checkMaxIterations(const size_t iteration,
                                         EndCriteria::Type& ecType) const{
        if (iteration < maxIterations_)
            return false;
        ecType = MaxIterations;
        return true;
    }
    
    bool EndCriteria::checkStationaryPoint(const RealType xOld,
                                           const RealType xNew,
                                           size_t& statStateIterations,
                                           EndCriteria::Type& ecType) const {
        if (std::fabs(xNew-xOld) >= rootEpsilon_) {
            statStateIterations = 0;
            return false;
        }
        ++statStateIterations;
        if (statStateIterations <= maxStationaryStateIterations_)
            return false;
        ecType = StationaryPoint;
        return true;
    }
    
    bool EndCriteria::checkStationaryFunctionValue(
                                                   const RealType fxOld,
                                                   const RealType fxNew,
                                                   size_t& statStateIterations,
                                                   EndCriteria::Type& ecType) const {
        if (std::fabs(fxNew-fxOld) >= functionEpsilon_) {
            statStateIterations = 0;
            return false;
        }
        ++statStateIterations;
        if (statStateIterations <= maxStationaryStateIterations_)
            return false;
        ecType = StationaryFunctionValue;
        return true;
    }
    
    bool EndCriteria::checkStationaryFunctionAccuracy(
                                                      const RealType f,
                                                      const bool positiveOptimization,
                                                      EndCriteria::Type& ecType) const {
        if (!positiveOptimization)
            return false;
        if (f >= functionEpsilon_)
            return false;
        ecType = StationaryFunctionAccuracy;
        return true;
    }
    
    //bool EndCriteria::checkZerGradientNormValue(
    //                                        const RealType gNormOld,
    //                                        const RealType gNormNew,
    //                                        EndCriteria::Type& ecType) const {
    //    if (std::fabs(gNormNew-gNormOld) >= gradientNormEpsilon_)
    //        return false;
    //    ecType = StationaryGradient;
    //    return true;
    //}
    
    bool EndCriteria::checkZeroGradientNorm(const RealType gradientNorm,
                                            EndCriteria::Type& ecType) const {
        if (gradientNorm >= gradientNormEpsilon_)
            return false;
        ecType = ZeroGradientNorm;
        return true;
    }
    
    bool EndCriteria::operator()(const size_t iteration,
                                 size_t& statStateIterations,
                                 const bool positiveOptimization,
                                 const RealType fold,
                                 const RealType, //normgold,
                                 const RealType fnew,
                                 const RealType normgnew,
                                 EndCriteria::Type& ecType) const {
        return
            checkMaxIterations(iteration, ecType) ||
            checkStationaryFunctionValue(fold, fnew, statStateIterations, ecType) ||
            checkStationaryFunctionAccuracy(fnew, positiveOptimization, ecType) ||
            checkZeroGradientNorm(normgnew, ecType);
    }
    
    // Inspectors
    size_t EndCriteria::maxIterations() const {
        return maxIterations_;
    }
    
    size_t EndCriteria::maxStationaryStateIterations() const {
        return maxStationaryStateIterations_;
    }
    
    RealType EndCriteria::rootEpsilon() const {
        return rootEpsilon_;
    }
    
    RealType EndCriteria::functionEpsilon() const {
        return functionEpsilon_;
    }
    
    RealType EndCriteria::gradientNormEpsilon() const {
        return gradientNormEpsilon_;
    }
    
    std::ostream& operator<<(std::ostream& out, EndCriteria::Type ec) {
        switch (ec) {
        case QuantLib::EndCriteria::None:
            return out << "None";
        case QuantLib::EndCriteria::MaxIterations:
            return out << "MaxIterations";
        case QuantLib::EndCriteria::StationaryPoint:
            return out << "StationaryPoint";
        case QuantLib::EndCriteria::StationaryFunctionValue:
            return out << "StationaryFunctionValue";
        case QuantLib::EndCriteria::StationaryFunctionAccuracy:
            return out << "StationaryFunctionAccuracy";
        case QuantLib::EndCriteria::ZeroGradientNorm:
            return out << "ZeroGradientNorm";
        case QuantLib::EndCriteria::Unknown:
            return out << "Unknown";
        default:
            sprintf(painCave.errMsg, "unknown EndCriteria::Type ( %d )\n",
                    int(ec));
            painCave.isFatal = 1;
            painCave.severity = OPENMD_ERROR;
            simError();
        }
        return out;
    }

}
