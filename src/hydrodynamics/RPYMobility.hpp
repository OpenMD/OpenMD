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

/*! \file hydrodynamics/RPYMobility.hpp
    \brief Per-molecule translation-only RPY mobility / resistance assembly.

    Builds, for one molecule of N spherical beads, the 3N x 3N grand
    translational mobility M (tt blocks only), inverts it to the grand
    resistance R = M^{-1}, and Cholesky-factors R for fluctuation-dissipation
    consistent random forces. The dipolar (td) rows enter only through the
    effective ambient velocity

        v_inf_eff_i = v_inf_i + sum_{j!=i} mu_td_ij : E ,

    so the hydrodynamic drag in the no-projection (J = 1) inertial scheme is

        f_hydro = R ( v_inf_eff - v ),

    and the random force is  sqrt(2 kT / dt) * S * Z,  S S^T = R, Z ~ N(0,1).

    This is the orientation-less (unstructured spherical atom) case: no
    rotational degrees of freedom, hence no tr/rt/rr/rd blocks. Rotational
    hydrodynamic coupling, which would feed the translational resistance
    through the full 6N inversion, is therefore neglected here by construction.
*/

#ifndef HYDRODYNAMICS_RPYMOBILITY_HPP
#define HYDRODYNAMICS_RPYMOBILITY_HPP

#include <vector>

#include "math/DynamicRectMatrix.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/Vector3.hpp"

namespace OpenMD {

  class RPYMobility {
  public:
    RPYMobility(const std::vector<RealType>& radii, RealType viscosity);

    //! (Re)build M, R = M^{-1}, and S = chol(R) from current bead positions.
    //! Positions must be the unwrapped intramolecular coordinates. Returns
    //! false if R is not positive definite (Cholesky hit a non-positive pivot).
    bool update(const std::vector<Vector3d>& positions);

    //! v_inf_eff_i = ambient_i + sum_{j!=i} mu_td_ij : E   (E symmetric traceless)
    std::vector<Vector3d> effectiveAmbient(
        const std::vector<Vector3d>& positions,
        const std::vector<Vector3d>& ambient, const Mat3x3d& E) const;

    //! Hydrodynamic drag: f_i = sum_j R_ij ( vEff_j - vel_j ).
    std::vector<Vector3d> dragForce(const std::vector<Vector3d>& vEff,
                                    const std::vector<Vector3d>& vel) const;

    //! FDT random force: sqrt(2 kT / dt) * S * Z, with Z a 3N vector of
    //! independent standard normal draws (row 3i+a is bead i, component a).
    std::vector<Vector3d> randomForce(const std::vector<RealType>& Z,
                                      RealType kT, RealType dt) const;

    //! Direct (non-iterative) solve of the implicit friction velocity update
    //!   (I + diag(c) R) x = rhs ,
    //! with c a per-bead coefficient (length N, applied to all three components
    //! of each bead). Unconditionally stable for any R >= 0, replacing the
    //! Picard iteration whose spectral radius c*lambda_max(R) exceeds 1 for
    //! stiff, overlapping-bead configurations.
    std::vector<Vector3d> solveImplicitVelocity(
        const std::vector<RealType>& c,
        const std::vector<Vector3d>& rhs) const;

    std::size_t size() const { return N_; }
    bool isPositiveDefinite() const { return spd_; }
    const DynamicRectMatrix<RealType>& mobility() const { return M_; }
    const DynamicRectMatrix<RealType>& resistance() const { return R_; }
    const DynamicRectMatrix<RealType>& cholesky() const { return S_; }

  private:
    std::size_t N_;
    RealType eta_;
    std::vector<RealType> a_;
    bool spd_ {false};
    DynamicRectMatrix<RealType> M_, R_, S_;
  };
}  // namespace OpenMD

#endif  // HYDRODYNAMICS_RPYMOBILITY_HPP
