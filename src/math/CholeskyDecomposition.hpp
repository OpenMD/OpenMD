/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
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
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *                                                                      
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */

#include "math/Vector.hpp"

#ifndef MATH_CHOLESKYDECOMPOSITION_HPP
#define MATH_CHOLESKYDECOMPOSITION_HPP

using namespace std;
namespace OpenMD {

  template<class MatrixType>
  void CholeskyDecomposition(MatrixType& A, MatrixType& L) {

    int n = A.getNRow();
    assert(n == A.getNCol() && n == L.getNRow() && n == L.getNCol());

    bool isspd(true);  
    RealType eps = A.diagonals().abs().max()  *
      (numeric_limits<RealType>::epsilon())/100;

  
    for(int j = 0; j < n; j++) {
      RealType d(0.0);
      for (int k = 0; k < j; k++) {
        RealType s(0.0);
      
        for (int i = 0; i < k; i++) {
          s += L(k,i) * L(j,i);
        }
      
        // if L(k,k) != 0
        if (std::abs(L(k,k)) > eps) {
          s = (A(j,k) - s) / L(k,k);
        } else {
          s = (A(j,k) -s);
          isspd = false;
        }
        L(j,k) = s;
        d = d + s*s;
      
        // this is approximately doing: isspd = isspd && ( A(k,j) == A(j,k) )
        isspd = isspd && (abs(A(k,j) - A(j,k)) < eps );
      }
      d = A(j,j) - d;
      isspd = isspd && (d > eps);
      L(j,j) = sqrt(d > 0.0 ? d : 0.0);
      for (int k = j+1; k < n; k++)  {
        L(j,k) = 0.0;
      }
    }
  }
}

#endif
