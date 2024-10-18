/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the following paper when you publish your work:
 *
 * [1] Drisko et al., J. Open Source Softw. 9, 7004 (2024).
 *
 * Good starting points for code and simulation methodology are:
 *
 * [2] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [3] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [4] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [7] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 * [9] Drisko & Gezelter, J. Chem. Theory Comput. 20, 4986-4997 (2024).
 */

#include "math/RealSphericalHarmonic.hpp"

#include <cmath>
#include <cstdio>
#include <limits>

using namespace OpenMD;

RealSphericalHarmonic::RealSphericalHarmonic() {}

RealType RealSphericalHarmonic::getValueAt(RealType costheta, RealType phi) {
  RealType p, phase;

  // associated Legendre polynomial
  p = LegendreP(L, M, costheta);

  if (functionType == RSH_SIN) {
    phase = sin((RealType)M * phi);
  } else {
    phase = cos((RealType)M * phi);
  }

  return coefficient * p * phase;
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
RealType RealSphericalHarmonic::LegendreP(int l, int m, RealType x) {
  // check parameters
  if (m < 0 || m > l || fabs(x) > 1.0) {
    printf("LegendreP got a bad argument: l = %d\tm = %d\tx = %lf\n", l, m, x);
    //    return NAN;
    return std::numeric_limits<RealType>::quiet_NaN();
  }

  RealType pmm = 1.0;
  if (m > 0) {
    RealType h = sqrt((1.0 - x) * (1.0 + x)), f = 1.0;
    for (int i = 1; i <= m; i++) {
      pmm *= -f * h;
      f += 2.0;
    }
  }
  if (l == m)
    return pmm;
  else {
    RealType pmmp1 = x * (2 * m + 1) * pmm;
    if (l == (m + 1))
      return pmmp1;
    else {
      RealType pll = 0.0;
      for (int ll = m + 2; ll <= l; ll++) {
        pll   = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
        pmm   = pmmp1;
        pmmp1 = pll;
      }
      return pll;
    }
  }
}
