#include <stdio.h>
#include <cmath>

#include "math/RealSphericalHarmonic.hpp"

using namespace oopse;

RealSphericalHarmonic::RealSphericalHarmonic() {
}

double RealSphericalHarmonic::getValueAt(double costheta, double phi) {
  
  double p, phase;
  
  // associated Legendre polynomial
  p = LegendreP(L,M,costheta);
 
  if (functionType == RSH_SIN) {
    phase = sin((double)M * phi);
  } else {
    phase = cos((double)M * phi);
  }
  
  return coefficient*p*phase;
  
}

//---------------------------------------------------------------------------//
//
// double LegendreP (int l, int m, double x);
//
// Computes the value of the associated Legendre polynomial P_lm (x)
// of order l at a given point.
//
// Input:
//   l  = degree of the polynomial  >= 0
//   m  = parameter satisfying 0 <= m <= l,
//   x  = point in which the computation is performed, range -1 <= x <= 1.
// Returns:
//   value of the polynomial in x
//
//---------------------------------------------------------------------------//
double RealSphericalHarmonic::LegendreP (int l, int m, double x) {
  // check parameters
  if (m < 0 || m > l || fabs(x) > 1.0) {
    printf("LegendreP got a bad argument: l = %d\tm = %d\tx = %lf\n", l, m, x);
    return NAN;
  }
  
  double pmm = 1.0;
  if (m > 0) {
    double h = sqrt((1.0-x)*(1.0+x)),
      f = 1.0;
    for (int i = 1; i <= m; i++) {
      pmm *= -f * h;
      f += 2.0;
    }
  }
  if (l == m)
    return pmm;
  else {
    double pmmp1 = x * (2 * m + 1) * pmm;
    if (l == (m+1))
      return pmmp1;
    else {
      double pll = 0.0;
      for (int ll = m+2; ll <= l; ll++) {
        pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
        pmm = pmmp1;
        pmmp1 = pll;
      }
      return pll;
    }
  }
}

