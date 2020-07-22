/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */
 
#ifndef MATH_REALSYMMETRICTRIDIAGONAL_HPP
#define MATH_REALSYMMETRICTRIDIAGONAL_HPP

#include "config.h"
#include "math/DynamicRectMatrix.hpp"

#include <algorithm>
// for min(), max() below
#include <cmath>
// for abs() below

namespace OpenMD {

  /** 
      
      Computes eigenvalues and eigenvectors of a real (non-complex)
      symmetric tridiagonal matrix by the QL method. 
  **/
  
  template<typename Real>
  class RealSymmetricTridiagonal {

    /** Row and column dimension (square matrix).  */
    int n;

    DynamicVector<Real> d;
    DynamicVector<Real> e;

    /** Array for internal storage of eigenvectors. */
    DynamicRectMatrix<Real> V;
    
    // Symmetric tridiagonal QL algorithm.
    
    void tql2 () {
      
      //  This is derived from the Algol procedures tql2, by
      //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
      //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
      //  Fortran subroutine in EISPACK.
      
      for (int i = 1; i < n; i++) {
        e(i-1) = e(i);
      }
      e(n-1) = 0.0;
      
      Real f = 0.0;
      Real tst1 = 0.0;
      Real eps = pow(2.0,-52.0);
      for (int l = 0; l < n; l++) {
        
        // Find small subdiagonal element
        
        tst1 = max(tst1,abs(d(l)) + abs(e(l)));
        int m = l;
        
        // Original while-loop from Java code
        while (m < n) {
          if (abs(e(m)) <= eps*tst1) {
            break;
          }
          m++;
        }
        
        
        // If m == l, d(l) is an eigenvalue,
        // otherwise, iterate.
        
        if (m > l) {
          int iter = 0;
          do {
            iter = iter + 1;  // (Could check iteration count here.)
            
            // Compute implicit shift
            
            Real g = d(l);
            Real p = (d(l+1) - g) / (2.0 * e(l));
            Real r = hypot(p,1.0);
            if (p < 0) {
              r = -r;
            }
            d(l) = e(l) / (p + r);
            d(l+1) = e(l) * (p + r);
            Real dl1 = d(l+1);
            Real h = g - d(l);
            for (int i = l+2; i < n; i++) {
              d(i) -= h;
            }
            f = f + h;
            
            // Implicit QL transformation.
            
            p = d(m);
            Real c = 1.0;
            Real c2 = c;
            Real c3 = c;
            Real el1 = e(l+1);
            Real s = 0.0;
            Real s2 = 0.0;
            for (int i = m-1; i >= l; i--) {
              c3 = c2;
              c2 = c;
              s2 = s;
              g = c * e(i);
              h = c * p;
              r = hypot(p,e(i));
              e(i+1) = s * r;
              s = e(i) / r;
              c = p / r;
              p = c * d(i) - s * g;
              d(i+1) = h + s * (c * g + s * d(i));
              
              // Accumulate transformation.
              
              for (int k = 0; k < n; k++) {
                h = V(k,i+1);
                V(k,i+1) = s * V(k,i) + c * h;
                V(k,i) = c * V(k,i) - s * h;
              }
            }
            p = -s * s2 * c3 * el1 * e(l) / dl1;
            e(l) = s * p;
            d(l) = c * p;
            
            // Check for convergence.
            
          } while (abs(e(l)) > eps*tst1);
        }
        d(l) = d(l) + f;
        e(l) = 0.0;
      }
      
      // Sort eigenvalues and corresponding vectors.
      
      for (int i = 0; i < n-1; i++) {
        int k = i;
        Real p = d(i);
        for (int j = i+1; j < n; j++) {
          if (d(j) < p) {
            k = j;
            p = d(j);
          }
        }
        if (k != i) {
          d(k) = d(i);
          d(i) = p;
          for (int j = 0; j < n; j++) {
            p = V(j,i);
            V(j,i) = V(j,k);
            V(j,k) = p;
          }
        }
      }
    }

  public:
    /** Construct the eigenvalue decomposition

        @param diagonals the diagonal elements of the input matrix.
        @param subdiagonals the subdiagonal elements of the input matrix in its
                 last n-1 positions.  subdiagonals[0] is arbitrary.
    */    
    RealSymmetricTridiagonal(const DynamicVector<Real> &diagonals, const DynamicVector<Real> &subdiagonals) {

      n = diagonals.size();
      V = DynamicRectMatrix<Real>(n,n,0.0);
      d = diagonals;
      e = subdiagonals;

      for (int i = 0; i < n; i++) {
        V(i,i) = 1.0;
      }
      // Diagonalize.
      tql2();
    }

    /** Return the eigenvector matrix
        @return     V
    */    
    void getEigenvectors (DynamicRectMatrix<Real> &V_) {
      V_ = V;
      return;
    }
    
    /** Return the real parts of the eigenvalues
        @return     real(diag(D))
    */    
    void getEigenvalues (DynamicVector<Real> &d_) {
      d_ = d;
      return ;
    }
  };
}
#endif //MATH_REALSYMMETRICTRIDIAGONAL_HPP
