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
 */

#include "hydrodynamics/RPYMobility.hpp"

#include <cmath>
#include <limits>

#include "hydrodynamics/RPYBlocks.hpp"
#include "math/DynamicVector.hpp"
// Pulls in LU.hpp + CholeskyDecomposition.hpp, and (under HAVE_LAPACK) the
// DynamicRectMatrix<RealType> overloads that dispatch to LAPACK. Without
// LAPACK this is exactly the generic templates.
#include "math/LapackLinearAlgebra.hpp"

namespace OpenMD {

  RPYMobility::RPYMobility(const std::vector<RealType>& radii,
                           RealType viscosity) :
      N_(radii.size()), eta_(viscosity), a_(radii) {
    std::size_t n3 = 3 * N_;
    M_             = DynamicRectMatrix<RealType>(n3, n3);
    R_             = DynamicRectMatrix<RealType>(n3, n3);
    S_             = DynamicRectMatrix<RealType>(n3, n3);
  }

  bool RPYMobility::update(const std::vector<Vector3d>& positions) {
    // --- assemble the 3N x 3N translational mobility from tt blocks ---
    for (std::size_t i = 0; i < N_; ++i) {
      for (std::size_t j = 0; j < N_; ++j) {
        Mat3x3d block =
            (i == j) ? RPY::muTT_self(a_[i], eta_)
                     : RPY::muTT(positions[i] - positions[j], a_[i], a_[j],
                                 eta_);
        // place the 3x3 block at (3i, 3j); confirm setSubMatrix arg order is
        // (beginRow, beginCol, block) in DynamicRectMatrix.
        M_.setSubMatrix(3 * i, 3 * j, block);
      }
    }

    // --- R = M^{-1} ---
    // invertMatrix overwrites its first argument with the LU factorization,
    // so invert a copy and keep M_ intact (HydroProp does the same).
    DynamicRectMatrix<RealType> Mwork(M_);
    if (!invertMatrix(Mwork, R_)) {
      spd_ = false;
      return false;
    }

    // --- S = chol(R), with S S^T = R (R is read, not modified) ---
    CholeskyDecomposition(R_, S_);

    // RPYC keeps M (and hence R = M^{-1}) SPD by construction, so a
    // non-positive Cholesky pivot is not a routine event: it signals an
    // upstream problem (NaN or degenerate coordinates, unphysical radii or
    // viscosity). Flag it here; the force modifier promotes it to a hard error
    // with molecule/step context rather than letting clamped (FDT-violating)
    // noise leak into the run.
    spd_          = true;
    RealType dmax = R_.diagonals().abs().max();
    RealType eps  = std::sqrt(dmax) * std::numeric_limits<RealType>::epsilon();
    for (std::size_t k = 0; k < 3 * N_; ++k)
      if (!(S_(k, k) > eps)) spd_ = false;  // false for 0, negative, and NaN

    return spd_;
  }

  std::vector<Vector3d> RPYMobility::effectiveAmbient(
      const std::vector<Vector3d>& positions,
      const std::vector<Vector3d>& ambient, const Mat3x3d& E) const {
    std::vector<Vector3d> vEff(ambient);
    for (std::size_t i = 0; i < N_; ++i) {
      for (std::size_t j = 0; j < N_; ++j) {
        if (i == j) continue;  // mu_td_ii : E = 0 (single sphere, linear flow)
        vEff[i] +=
            RPY::muTD_dot_E(positions[i] - positions[j], a_[i], a_[j], E);
      }
    }
    return vEff;
  }

  std::vector<Vector3d> RPYMobility::dragForce(
      const std::vector<Vector3d>& vEff,
      const std::vector<Vector3d>& vel) const {
    // f = R ( vEff - vel ), as one 3N matrix-vector product. Flatten the
    // per-bead velocity mismatch into a 3N vector, apply R, and unpack. When
    // the DynamicRectMatrix gemv is specialized onto BLAS (see LU/Cholesky
    // LAPACK path) this is a single dgemv.
    DynamicVector<RealType> dv(3 * N_);
    for (std::size_t i = 0; i < N_; ++i)
      for (int p = 0; p < 3; ++p)
        dv[3 * i + p] = vEff[i][p] - vel[i][p];

    DynamicVector<RealType> fv = R_ * dv;  // confirm operator* (mat * vec)

    std::vector<Vector3d> f(N_, V3Zero);
    for (std::size_t i = 0; i < N_; ++i)
      for (int p = 0; p < 3; ++p)
        f[i][p] = fv[3 * i + p];
    return f;
  }

  std::vector<Vector3d> RPYMobility::randomForce(const std::vector<RealType>& Z,
                                                 RealType kT,
                                                 RealType dt) const {
    // F^R = sqrt(2 kT / dt) * S * Z, with S the lower-triangular Cholesky
    // factor (S S^T = R). CholeskyDecomposition explicitly zeroes the upper
    // triangle of S, so a full matrix-vector product is exact; a triangular
    // multiply would save ~2x flops but this is O(N^2) against the O(N^3)
    // factorization, so it is not worth a special path.
    DynamicVector<RealType> z(3 * N_);
    for (std::size_t k = 0; k < 3 * N_; ++k)
      z[k] = Z[k];

    DynamicVector<RealType> fv = S_ * z;
    RealType scale             = std::sqrt(2.0 * kT / dt);

    std::vector<Vector3d> f(N_, V3Zero);
    for (std::size_t i = 0; i < N_; ++i)
      for (int p = 0; p < 3; ++p)
        f[i][p] = scale * fv[3 * i + p];
    return f;
  }

  std::vector<Vector3d> RPYMobility::solveImplicitVelocity(
      const std::vector<RealType>& c,
      const std::vector<Vector3d>& rhs) const {
    // Build A = I + diag(c) R (c_i scales the three rows of bead i) and solve
    // A x = rhs. A has eigenvalues 1 + c*lambda(R) >= 1, so it is nonsingular
    // for any positive-semidefinite R -- the stiff (overlapping) modes simply
    // come back strongly damped (x -> rhs's flow part) instead of diverging.
    // A is non-symmetric when bead masses differ, so this uses the LU path
    // (invertMatrix), which dispatches to LAPACK when available. One extra
    // O(N^3) solve per step; if it ever matters, factor A once (dgetrf/dgetrs)
    // rather than forming the explicit inverse.
    const std::size_t n3 = 3 * N_;
    DynamicRectMatrix<RealType> A(n3, n3);
    for (std::size_t i = 0; i < N_; ++i) {
      for (int p = 0; p < 3; ++p) {
        std::size_t row = 3 * i + p;
        for (std::size_t col = 0; col < n3; ++col)
          A(row, col) = c[i] * R_(row, col);
        A(row, row) += 1.0;  // + I
      }
    }
 
    DynamicRectMatrix<RealType> Ainv(n3, n3);
    invertMatrix(A, Ainv);  // A is overwritten; inverse returned in Ainv
 
    DynamicVector<RealType> b(n3);
    for (std::size_t i = 0; i < N_; ++i)
      for (int p = 0; p < 3; ++p)
        b[3 * i + p] = rhs[i][p];
 
    DynamicVector<RealType> x = Ainv * b;
 
    std::vector<Vector3d> velStep(N_, V3Zero);
    for (std::size_t i = 0; i < N_; ++i)
      for (int p = 0; p < 3; ++p)
        velStep[i][p] = x[3 * i + p];
    return velStep;
  }

}  // namespace OpenMD
