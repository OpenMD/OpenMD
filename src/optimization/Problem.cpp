#include "optimization/Problem.hpp"

#include <config.h>

namespace QuantLib {
  RealType Problem::DotProduct(DynamicVector<RealType>& v1,
                               DynamicVector<RealType>& v2) {
    RealType dp = dot(v1, v2);
    return dp;
  }

  RealType Problem::computeGradientNormValue(DynamicVector<RealType>& grad_f) {
    RealType dot = grad_f.lengthSquare();
    return dot;
  }
}  // namespace QuantLib
