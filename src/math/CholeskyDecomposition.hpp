/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
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
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

#include "math/Vector.hpp"

#ifndef MATH_CHOLESKYDECOMPOSITION_HPP
#define MATH_CHOLESKYDECOMPOSITION_HPP

using namespace std;
namespace OpenMD {

  template<class MatrixType>
  void CholeskyDecomposition(MatrixType& A, MatrixType& L) {
    unsigned int n = A.getNRow();
    assert(n == A.getNCol() && n == L.getNRow() && n == L.getNCol());

    bool isspd(true);
    RealType eps =
        A.diagonals().abs().max() * (numeric_limits<RealType>::epsilon()) / 100;

    for (unsigned int j = 0; j < n; j++) {
      RealType d(0.0);
      for (unsigned int k = 0; k < j; k++) {
        RealType s(0.0);

        for (unsigned int i = 0; i < k; i++) {
          s += L(k, i) * L(j, i);
        }

        // if L(k,k) != 0
        if (std::abs(L(k, k)) > eps) {
          s = (A(j, k) - s) / L(k, k);
        } else {
          s     = (A(j, k) - s);
          isspd = false;
        }
        L(j, k) = s;
        d       = d + s * s;

        // this is approximately doing: isspd = isspd && ( A(k,j) == A(j,k) )
        isspd = isspd && (abs(A(k, j) - A(j, k)) < eps);
      }
      d       = A(j, j) - d;
      isspd   = isspd && (d > eps);
      L(j, j) = sqrt(d > 0.0 ? d : 0.0);
      for (unsigned int k = j + 1; k < n; k++) {
        L(j, k) = 0.0;
      }
    }
  }
}  // namespace OpenMD

#endif
