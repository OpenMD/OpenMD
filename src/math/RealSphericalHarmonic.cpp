/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * This software is provided "AS IS," without a warranty of any
 * kind. All express or implied conditions, representations and
 * warranties, including any implied warranty of merchantability,
 * fitness for a particular purpose or non-infringement, are hereby
 * excluded.  The University of Notre Dame and its licensors shall not
 * be liable for any damages suffered by licensee as a result of
 * using, modifying or distributing the software or its
 * derivatives. In no event will the University of Notre Dame or its
 * licensors be liable for any lost revenue, profit or data, or for
 * direct, indirect, special, consequential, incidental or punitive
 * damages, however caused and regardless of the theory of liability,
 * arising out of the use of or inability to use software, even if the
 * University of Notre Dame has been advised of the possibility of
 * such damages.
 */
 
#include <stdio.h>
#include <cmath>
#include <limits>
#include "math/RealSphericalHarmonic.hpp"

using namespace oopse;

RealSphericalHarmonic::RealSphericalHarmonic() {
}

RealType RealSphericalHarmonic::getValueAt(RealType costheta, RealType phi) {
  
  RealType p, phase;
  
  // associated Legendre polynomial
  p = LegendreP(L,M,costheta);
 
  if (functionType == RSH_SIN) {
    phase = sin((RealType)M * phi);
  } else {
    phase = cos((RealType)M * phi);
  }
  
  return coefficient*p*phase;
  
}

//---------------------------------------------------------------------------//
//
// RealType LegendreP (int l, int m, RealType x);
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
RealType RealSphericalHarmonic::LegendreP (int l, int m, RealType x) {
  // check parameters
  if (m < 0 || m > l || fabs(x) > 1.0) {
    printf("LegendreP got a bad argument: l = %d\tm = %d\tx = %lf\n", l, m, x);
//    return NAN;
	return std::numeric_limits <RealType>:: quiet_NaN();
  }
  
  RealType pmm = 1.0;
  if (m > 0) {
    RealType h = sqrt((1.0-x)*(1.0+x)),
      f = 1.0;
    for (int i = 1; i <= m; i++) {
      pmm *= -f * h;
      f += 2.0;
    }
  }
  if (l == m)
    return pmm;
  else {
    RealType pmmp1 = x * (2 * m + 1) * pmm;
    if (l == (m+1))
      return pmmp1;
    else {
      RealType pll = 0.0;
      for (int ll = m+2; ll <= l; ll++) {
        pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
        pmm = pmmp1;
        pmmp1 = pll;
      }
      return pll;
    }
  }
}

