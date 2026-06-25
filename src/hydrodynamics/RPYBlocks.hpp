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

/*! \file hydrodynamics/RPYBlocks.hpp
    \brief Pairwise Rotne-Prager-Yamakawa mobility blocks for a bead model.

    The six mobility blocks needed to embed a flexible atomic bead model in a
    linearized background flow:

      tt, tr, rt, rr   force/torque   -> translation/rotation   (carry 1/eta)
      td, rd           ambient strain -> translation/rotation   (eta- AND
                                                                  pi-INDEPENDENT)

    The grand velocity mobility M (the 6N x 6N assembly of tt/tr/rt/rr) is
    inverted to the grand resistance R = M^{-1}; the drag on the beads is

      f_hydro = R ( v_inf_eff - v ),   v_inf_eff = v_inf + G : E_inf,

    with v_inf the affine flow at the bead (VelocityField::getVelocity) and
    (G:E_inf)_i = sum_{j!=i}[mu_td_ij:E ; mu_rd_ij:E] the dipolar disturbance.

    Three regimes per pair, with boundaries:

      far field       R > a_i + a_j
      partial overlap |a_i - a_j| < R <= a_i + a_j
      full immersion  R <= |a_i - a_j|

    VISCOSITY: tt/tr/rt/rr carry 1/eta; td/rd are eta- AND
    pi-independent, because the 1/(pi eta) of the stresslet propagator
    cancels the (pi eta) of a sphere's stresslet response to
    strain. The (5/6) a_j^3 / R^2 structure of the far-field td block
    is exactly that cancellation.

    VERIFIED against ZWMS (Zuk, Wajnryb, Mizerski, Szymczak, J. Fluid
    Mech. 741, R5, 2014), Eqs. 3.1-3.9, in all three regimes, and
    continuity across both regime boundaries confirmed
    numerically. The full-immersion rt/tr prefactor (see rtPrefactor)
    follows Eq. 3.5, R/(8 pi eta a_i^3). The td/rd tensors are defined
    only up to their symmetric, traceless part, which is harmless here
    since they are always contracted with the symmetric traceless E.
*/

#ifndef HYDRODYNAMICS_RPYBLOCKS_HPP
#define HYDRODYNAMICS_RPYBLOCKS_HPP

#include <cmath>

#include "math/SquareMatrix3.hpp"
#include "math/Vector3.hpp"
#include "utils/Constants.hpp"

namespace OpenMD::RPY {

  enum Regime { FarField, PartialOverlap, FullImmersion };

  //! Classify a pair by centre-to-centre distance R and radii a_i, a_j.
  inline Regime classify(RealType R, RealType ai, RealType aj) {
    if (R > ai + aj) return FarField;
    if (R > std::fabs(ai - aj)) return PartialOverlap;
    return FullImmersion;
  }

  //! cross-product (skew) matrix S(v) such that S(v) * w = v x w
  inline Mat3x3d skew(const Vector3d& v) {
    Mat3x3d S(0.0);
    S(0, 1) = -v.z();
    S(0, 2) = v.y();
    S(1, 0) = v.z();
    S(1, 2) = -v.x();
    S(2, 0) = -v.y();
    S(2, 1) = v.x();
    return S;
  }

  //! Heaviside step used by the full-immersion rt/td branches: 1 if x>0 else 0.
  inline RealType heavi(RealType x) { return (x > 0.0) ? 1.0 : 0.0; }

  // ---- self blocks (i == j) --------------------------------------------- //

  inline Mat3x3d muTT_self(RealType a, RealType eta) {
    return (1.0 / (6.0 * Constants::PI * eta * a)) *
           SquareMatrix3<RealType>::identity();
  }
  inline Mat3x3d muRR_self(RealType a, RealType eta) {
    return (1.0 / (8.0 * Constants::PI * eta * a * a * a)) *
           SquareMatrix3<RealType>::identity();
  }
  // mu_tr_self = mu_rt_self = 0; mu_td_self : E = mu_rd_self : E = 0.

  // ---- off-diagonal blocks (rij = r_i - r_j, rhat = rij/|rij|) ----------- //

  //! translation-translation, carries 1/eta
  inline Mat3x3d muTT(const Vector3d& rij, RealType ai, RealType aj,
                      RealType eta) {
    RealType R   = rij.length();
    Vector3d rh  = rij / R;
    RealType R2  = R * R, R3 = R2 * R;
    Mat3x3d I    = SquareMatrix3<RealType>::identity();
    Mat3x3d op   = outProduct(rh, rh);
    RealType eye = 0.0, hat = 0.0;

    switch (classify(R, ai, aj)) {
    case FarField: {
      RealType c    = 1.0 / (8.0 * Constants::PI * eta * R);
      RealType asum = (ai * ai + aj * aj) / R2;
      eye           = c * (1.0 + asum / 3.0);
      hat           = c * (1.0 - asum);
      break;
    }
    case PartialOverlap: {
      RealType d2   = (ai - aj) * (ai - aj);
      RealType invc = 1.0 / (6.0 * Constants::PI * eta * ai * aj);
      eye = ((16.0 * R3 * (ai + aj) - (d2 + 3.0 * R2) * (d2 + 3.0 * R2)) /
             (32.0 * R3)) *
            invc;
      hat = (3.0 * (d2 - R2) * (d2 - R2) / (32.0 * R3)) * invc;
      break;
    }
    case FullImmersion: {
      RealType asup = std::max(ai, aj);
      eye           = 1.0 / (6.0 * Constants::PI * eta * asup);
      hat           = 0.0;
      break;
    }
    }
    return eye * I + hat * op;
  }

  //! rotation-rotation, carries 1/eta
  inline Mat3x3d muRR(const Vector3d& rij, RealType ai, RealType aj,
                      RealType eta) {
    RealType R   = rij.length();
    Vector3d rh  = rij / R;
    RealType R2  = R * R, R3 = R2 * R, R4 = R2 * R2, R6 = R3 * R3;
    Mat3x3d I    = SquareMatrix3<RealType>::identity();
    Mat3x3d op   = outProduct(rh, rh);
    RealType eye = 0.0, hat = 0.0;

    switch (classify(R, ai, aj)) {
    case FarField: {
      RealType c = 1.0 / (16.0 * Constants::PI * eta * R3);
      eye        = -c;
      hat        = 3.0 * c;
      break;
    }
    case PartialOverlap: {
      RealType d2   = (ai - aj) * (ai - aj);
      RealType d4   = d2 * d2;
      RealType calA = (5.0 * R6 - 27.0 * R4 * (ai * ai + aj * aj) +
                       32.0 * R3 * (ai * ai * ai + aj * aj * aj) -
                       9.0 * R2 * (ai * ai - aj * aj) * (ai * ai - aj * aj) -
                       d4 * (ai * ai + 4.0 * ai * aj + aj * aj)) /
                      (64.0 * R3);
      RealType calB = (3.0 * (d2 - R2) * (d2 - R2) *
                       (ai * ai + 4.0 * ai * aj + aj * aj - R2)) /
                      (64.0 * R3);
      RealType invc =
          1.0 / (8.0 * Constants::PI * eta * ai * ai * ai * aj * aj * aj);
      eye = calA * invc;
      hat = calB * invc;
      break;
    }
    case FullImmersion: {
      RealType asup = std::max(ai, aj);
      eye           = 1.0 / (8.0 * Constants::PI * eta * asup * asup * asup);
      hat           = 0.0;
      break;
    }
    }
    return eye * I + hat * op;
  }

  //! rt prefactor pf such that mu_rt = -pf * skew(rhat) (Omega_i from F_j).
  inline RealType rtPrefactor(RealType R, RealType ai, RealType aj,
                              RealType eta) {
    RealType R2 = R * R;
    switch (classify(R, ai, aj)) {
    case FarField:
      return 1.0 / (8.0 * Constants::PI * eta * R2);
    case PartialOverlap: {
      RealType num = (ai - aj + R) * (ai - aj + R) *
                     (aj * aj + 2.0 * aj * (ai + R) - 3.0 * (ai - R) * (ai - R));
      return (num / (8.0 * R2)) /
             (16.0 * Constants::PI * eta * ai * ai * ai * aj);
    }
    case FullImmersion:
      // ZWMS Eq. 3.5: theta(a_i - a_j) * R / zeta_i^rr, zeta_i^rr = 8 pi eta a_i^3.
      return heavi(ai - aj) * R / (8.0 * Constants::PI * eta * ai * ai * ai);
    }
    return 0.0;
  }

  //! rotation-translation: Omega_i from F_j. carries 1/eta.
  inline Mat3x3d muRT(const Vector3d& rij, RealType ai, RealType aj,
                      RealType eta) {
    RealType R  = rij.length();
    Vector3d rh = rij / R;
    return (-rtPrefactor(R, ai, aj, eta)) * skew(rh);
  }

  //! translation-rotation: U_i from T_j. Uses the i<->j swapped prefactor
  //! (equal to muRT only in the far field). carries 1/eta.
  inline Mat3x3d muTR(const Vector3d& rij, RealType ai, RealType aj,
                      RealType eta) {
    RealType R  = rij.length();
    Vector3d rh = rij / R;
    return (-rtPrefactor(R, aj, ai, eta)) * skew(rh);
  }

  // ---- dipolar blocks, contracted with the (uniform) strain E ----------- //
  // eta- AND pi-INDEPENDENT. E must be symmetric and traceless.

  //! mu_td_ij : E -> translational disturbance velocity at bead i
  //! = hat1 (E rhat) + hat3 (rhat . E . rhat) rhat
  inline Vector3d muTD_dot_E(const Vector3d& rij, RealType ai, RealType aj,
                             const Mat3x3d& E) {
    RealType R   = rij.length();
    Vector3d rh  = rij / R;
    RealType R2  = R * R, R4 = R2 * R2, R5 = R4 * R, R6 = R4 * R2;
    RealType ai2 = ai * ai, aj2 = aj * aj;
    RealType hat1 = 0.0, hat3 = 0.0;

    switch (classify(R, ai, aj)) {
    case FarField:
      hat1 = (5.0 * aj / 6.0) * (-2.0 * (5.0 * ai2 * aj2 + 3.0 * aj2 * aj2)) /
             (5.0 * R4);
      hat3 = (5.0 * aj / 6.0) * aj2 * (5.0 * ai2 + 3.0 * aj2 - 3.0 * R2) / R4;
      break;
    case PartialOverlap: {
      RealType d  = ai - aj;
      RealType dm = aj - ai;
      RealType dm5 =
          dm * dm * dm * dm * dm;  // (a_j - a_i)^5
      RealType calC = (10.0 * R6 - 24.0 * R5 * ai - 15.0 * R4 * (aj2 - ai2) +
                       dm5 * (ai + 5.0 * aj)) /
                      (40.0 * ai * aj * R4);
      RealType calD = ((d * d - R2) * (d * d - R2) *
                       (d * (ai + 5.0 * aj) - R2)) /
                      (16.0 * ai * aj * R4);
      hat1 = (5.0 * aj / 6.0) * calC;
      hat3 = (5.0 * aj / 6.0) * calD;
      break;
    }
    case FullImmersion:
      hat1 = -heavi(aj - ai) * R;
      hat3 = 0.0;
      break;
    }

    Vector3d Erh = E * rh;
    RealType rEr = dot(rh, Erh);
    return hat1 * Erh + (hat3 * rEr) * rh;
  }

  //! mu_rd_ij : E -> angular-velocity disturbance at bead i
  //! = p [ (E rhat) x rhat ]
  inline Vector3d muRD_dot_E(const Vector3d& rij, RealType ai, RealType aj,
                             const Mat3x3d& E) {
    RealType R  = rij.length();
    Vector3d rh = rij / R;
    RealType R2 = R * R, R3 = R2 * R;
    RealType p  = 0.0;

    switch (classify(R, ai, aj)) {
    case FarField:
      p = -2.5 * (aj / R) * (aj / R) * (aj / R);
      break;
    case PartialOverlap: {
      RealType d2   = (ai - aj) * (ai - aj);
      RealType calB = (3.0 * (d2 - R2) * (d2 - R2) *
                       (ai * ai + 4.0 * ai * aj + aj * aj - R2)) /
                      (64.0 * R3);
      p = -5.0 * calB / (3.0 * ai * ai * ai);
      break;
    }
    case FullImmersion:
      p = 0.0;
      break;
    }

    Vector3d Erh = E * rh;
    return p * cross(Erh, rh);
  }

}  // namespace OpenMD::RPY

#endif  // HYDRODYNAMICS_RPYBLOCKS_HPP
