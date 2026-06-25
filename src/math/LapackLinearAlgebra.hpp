/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved. (BSD 3-clause; full text omitted in this sketch.)
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite:
 * [1] Drisko et al., J. Open Source Softw. 9, 7004 (2024).
 */

/*! \file math/LapackLinearAlgebra.hpp
    \brief LAPACK-backed fast paths for the dynamic, double-precision dense
           routines, dispatched by overload.

    SKETCH / PROPOSAL -- conventions worked out, not yet build-tested.

    The generic templates in LU.hpp and CholeskyDecomposition.hpp are kept as
    the fallback for small fixed-size matrices (Mat3x3d, Mat6x6d in HydroProp,
    ...) and for non-RealType element types, where the packing + call overhead
    of going out to LAPACK would lose to the inlined VTK/Crout code. Here we add
    *non-template overloads* for exactly DynamicRectMatrix<RealType>: a
    non-template overload is preferred over a function template for an exact
    type match, so any call invertMatrix(A, AI) / CholeskyDecomposition(A, L)
    with A,L of that concrete type routes here, while every fixed-size caller
    keeps the generic template unchanged. No edits to the generic headers.

    Include this header instead of (or after) LU.hpp/CholeskyDecomposition.hpp
    at the call site (e.g. in RPYMobility.cpp) so the overloads are visible.
    With HAVE_LAPACK undefined it simply forwards to the generics.

    ASSUMPTIONS to confirm against the tree:
      - DynamicRectMatrix<RealType> stores its data contiguously and exposes a
        raw pointer (assumed getArray() below); if the accessor differs, only
        the two `a = ...getArray()` lines change.
      - Storage is row-major (OpenMD convention). Since the matrices we factor
        here (the RPY mobility M and resistance R) are symmetric, the row/column
        major mismatch with LAPACK is harmless, and for the general inverse the
        "invert the transpose, read back transposed" identity makes a plain
        dgetrf/dgetri on the row-major buffer return the correct inverse.
      - RealType is double. For a single-precision build (RealType == float),
        dispatch the s* entry points instead (a thin if-constexpr or a pair of
        overloaded extern "C" wrappers); only the BLAS/LAPACK names change.
*/

#ifndef MATH_LAPACKLINEARALGEBRA_HPP
#define MATH_LAPACKLINEARALGEBRA_HPP

#include <algorithm>
#include <vector>

#include "config.h"  // provides HAVE_LAPACK
#include "math/CholeskyDecomposition.hpp"
#include "math/DynamicRectMatrix.hpp"
#include "math/LU.hpp"

#ifdef HAVE_LAPACK

extern "C" {
// LU factor / inverse (general)
void dgetrf_(const int* m, const int* n, double* a, const int* lda, int* ipiv,
             int* info);
void dgetri_(const int* n, double* a, const int* lda, const int* ipiv,
             double* work, const int* lwork, int* info);
// Cholesky factor (symmetric positive definite)
void dpotrf_(const char* uplo, const int* n, double* a, const int* lda,
             int* info);
}

namespace OpenMD {

  /**
   * R = A^{-1} for a dynamic double matrix, via LAPACK LU.
   * Non-template overload: chosen ahead of LU.hpp's invertMatrix template for
   * this exact type. A is left intact; the result is written into AI. (Feeding
   * the row-major buffer to column-major LAPACK factors A^T, whose inverse read
   * back row-major is exactly A^{-1}, so no explicit transpose is needed.)
   */
  inline bool invertMatrix(DynamicRectMatrix<RealType>& A,
                           DynamicRectMatrix<RealType>& AI) {
    const int n = static_cast<int>(A.getNRow());
    if (n != static_cast<int>(A.getNCol()) ||
        n != static_cast<int>(AI.getNRow()) ||
        n != static_cast<int>(AI.getNCol()))
      return false;
    if (n == 0) return true;

    // factor/invert in AI's buffer, leaving A untouched
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        AI(i, j) = A(i, j);

    RealType* a = AI.getArrayPointer();  // contiguous row-major
    std::vector<int> ipiv(n);
    int info = 0, lda = n;

    dgetrf_(&n, &n, a, &lda, ipiv.data(), &info);
    if (info != 0) return false;  // singular (info>0) or bad arg (info<0)

    // workspace query, then inversion
    int lwork = -1;
    RealType wq = 0.0;
    dgetri_(&n, a, &lda, ipiv.data(), &wq, &lwork, &info);
    lwork = std::max(1, static_cast<int>(wq));
    std::vector<RealType> work(static_cast<std::size_t>(lwork));
    dgetri_(&n, a, &lda, ipiv.data(), work.data(), &lwork, &info);
    return info == 0;
  }

  /**
   * L L^T = A for a dynamic, symmetric, double matrix, via LAPACK dpotrf.
   * Returns the lower-triangular factor with the strict upper triangle zeroed,
   * matching the generic CholeskyDecomposition contract (so a downstream S*Z
   * may be a full matrix-vector product).
   *
   * Convention: our buffer is row-major; LAPACK reads it column-major, i.e. as
   * A^T. A is symmetric here (A == A^T), so requesting the 'U' factor from
   * LAPACK (column-major upper) yields, read back row-major, the lower factor
   * we want.
   */
  inline void CholeskyDecomposition(DynamicRectMatrix<RealType>& A,
                                    DynamicRectMatrix<RealType>& L) {
    const int n = static_cast<int>(A.getNRow());
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        L(i, j) = A(i, j);
    if (n == 0) return;

    RealType* l = L.getArrayPointer();  // contiguous row-major
    int info = 0, lda = n;
    const char uplo = 'U';  // row-major lower factor <-> column-major upper
    dpotrf_(&uplo, &n, l, &lda, &info);

    if (info != 0) {
      // not positive definite: zero the factor so the caller's pivot check
      // (S(k,k) > eps) trips deterministically instead of reading partial data.
      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
          L(i, j) = 0.0;
      return;
    }
    // dpotrf leaves the opposite triangle untouched; zero it to match the
    // generic routine, which returns a clean lower-triangular factor.
    for (int i = 0; i < n; ++i)
      for (int j = i + 1; j < n; ++j)
        L(i, j) = 0.0;
  }

}  // namespace OpenMD

#endif  // HAVE_LAPACK
#endif  // MATH_LAPACKLINEARALGEBRA_HPP
