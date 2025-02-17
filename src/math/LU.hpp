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

/*=========================================================================

  Program:   Visualization Toolkit
  Module:    Excerpted from vtkMath.cxx

  Copyright (c) 1993-2015 Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  * Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.

  * Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

  * Neither name of Ken Martin, Will Schroeder, or Bill Lorensen nor the names
    of any contributors may be used to endorse or promote products derived
    from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

  =========================================================================*/
#ifndef MATH_LU_HPP
#define MATH_LU_HPP

#include <iostream>
#include <limits>

namespace OpenMD {

  constexpr double SMALL_NUMBER        = 1.0e-12;
  constexpr int MAX_SCRATCH_ARRAY_SIZE = 10;

  /**
   * Invert input square matrix A into matrix AI.
   * @param A input square matrix
   * @param AI output square matrix
   * @return true if inverse is computed, otherwise return false
   * @note A is modified during the inversion
   */
  template<class MatrixType>
  bool invertMatrix(MatrixType& A, MatrixType& AI) {
    using Real = typename MatrixType::ElemType;

    if (A.getNRow() != A.getNCol() || A.getNRow() != AI.getNRow() ||
        A.getNCol() != AI.getNCol()) {
      return false;
    }

    int size = A.getNRow();

    // Check on allocation of working vectors
    //
    int iScratch[MAX_SCRATCH_ARRAY_SIZE];
    int* index = (size <= MAX_SCRATCH_ARRAY_SIZE ? iScratch : new int[size]);
    Real dScratch[MAX_SCRATCH_ARRAY_SIZE];
    Real* column = (size <= MAX_SCRATCH_ARRAY_SIZE ? dScratch : new Real[size]);

    bool retVal = invertMatrix(A, AI, size, index, column);

    if (size > MAX_SCRATCH_ARRAY_SIZE) {
      delete[] index;
      delete[] column;
    }
    return retVal;
  }

  /**
   * Invert input square matrix A into matrix AI (Thread safe versions).
   * @param A input square matrix
   * @param AI output square matrix
   * @param size size of the matrix and temporary arrays
   * @param tmp1Size temporary array
   * @param tmp2Size temporary array
   * @return true if inverse is computed, otherwise return false
   * @note A is modified during the inversion.
   */
  template<class MatrixType>
  bool invertMatrix(MatrixType& A, MatrixType& AI, unsigned int size,
                    int* tmp1Size,
                    typename MatrixType::ElemPoinerType tmp2Size) {
    if (A.getNRow() != A.getNCol() || A.getNRow() != AI.getNRow() ||
        A.getNCol() != AI.getNCol() || A.getNRow() != size) {
      return false;
    }

    unsigned int i, j;

    //
    // Factor matrix; then begin solving for inverse one column at a time.
    // Note: tmp1Size returned value is used later, tmp2Size is just working
    // memory whose values are not used in LUSolveLinearSystem
    //
    if (LUFactorLinearSystem(A, tmp1Size, size, tmp2Size) == 0) {
      return false;
    }

    for (j = 0; j < size; ++j) {
      for (i = 0; i < size; ++i) {
        tmp2Size[i] = 0.0;
      }
      tmp2Size[j] = 1.0;

      LUSolveLinearSystem(A, tmp1Size, tmp2Size, size);

      for (i = 0; i < size; i++) {
        AI(i, j) = tmp2Size[i];
      }
    }

    return true;
  }

  /**
   * Factor linear equations Ax = b using LU decompostion A = LU where L is
   * lower triangular matrix and U is upper triangular matrix.
   * @param A input square matrix
   * @param index pivot indices
   * @param size size of the matrix and temporary arrays
   * @param tmpSize temporary array
   * @return true if inverse is computed, otherwise return false
   * @note A is modified during the inversion.
   */
  template<class MatrixType>
  int LUFactorLinearSystem(MatrixType& A, int* index, int size,
                           typename MatrixType::ElemPoinerType tmpSize) {
    using Real = typename MatrixType::ElemType;

    int i, j, k;
    int maxI = 0;
    Real largest, temp1, temp2, sum;

    //
    // Loop over rows to get implicit scaling information
    //
    for (i = 0; i < size; ++i) {
      for (largest = 0.0, j = 0; j < size; ++j) {
        if ((temp2 = std::abs(A(i, j))) > largest) { largest = temp2; }
      }

      if (largest == 0.0) {
        std::cerr << "Unable to factor linear system";
        return 0;
      }
      tmpSize[i] = 1.0 / largest;
    }
    //
    // Loop over all columns using Crout's method
    //
    for (j = 0; j < size; ++j) {
      for (i = 0; i < j; ++i) {
        sum = A(i, j);
        for (k = 0; k < i; ++k) {
          sum -= A(i, k) * A(k, j);
        }
        A(i, j) = sum;
      }
      //
      // Begin search for largest pivot element
      //
      for (largest = 0.0, i = j; i < size; ++i) {
        sum = A(i, j);
        for (k = 0; k < j; ++k) {
          sum -= A(i, k) * A(k, j);
        }
        A(i, j) = sum;

        if ((temp1 = tmpSize[i] * std::abs(sum)) >= largest) {
          largest = temp1;
          maxI    = i;
        }
      }
      //
      // Check for row interchange
      //
      if (j != maxI) {
        for (k = 0; k < size; ++k) {
          std::swap(A(maxI, k), A(j, k));
        }
        tmpSize[maxI] = tmpSize[j];
      }
      //
      // Divide by pivot element and perform elimination
      //
      index[j] = maxI;

      if (std::abs(A(j, j)) <= SMALL_NUMBER) {
        std::cerr << "Unable to factor linear system";
        return false;
      }

      if (j != (size - 1)) {
        temp1 = 1.0 / A(j, j);
        for (i = j + 1; i < size; ++i) {
          A(i, j) *= temp1;
        }
      }
    }

    return 1;
  }

  /**
   * Solve linear equations Ax = b using LU decompostion A = LU where L is
   * lower triangular matrix and U is upper triangular matrix.
   * @param A input square matrix
   * @param index pivot indices
   * @param x vector
   * @param size size of the matrix and temporary arrays
   * @return true if inverse is computed, otherwise return false
   * @note A=LU and index[] are generated from method LUFactorLinearSystem).
   * Also, solution vector is written directly over input load vector.
   */
  template<class MatrixType>
  void LUSolveLinearSystem(MatrixType& A, int* index,
                           typename MatrixType::ElemPoinerType x, int size) {
    using Real = typename MatrixType::ElemType;

    int i, j, ii, idx;
    Real sum;
    //
    // Proceed with forward and backsubstitution for L and U
    // matrices.  First, forward substitution.
    //
    for (ii = -1, i = 0; i < size; ++i) {
      idx    = index[i];
      sum    = x[idx];
      x[idx] = x[i];

      if (ii >= 0) {
        for (j = ii; j <= (i - 1); ++j) {
          sum -= A(i, j) * x[j];
        }
      } else if (sum != 0.0) {
        ii = i;
      }

      x[i] = sum;
    }
    //
    // Now, back substitution
    //
    for (i = size - 1; i >= 0; i--) {
      sum = x[i];
      for (j = i + 1; j < size; ++j) {
        sum -= A(i, j) * x[j];
      }
      x[i] = sum / A(i, i);
    }
  }
}  // namespace OpenMD

#endif
