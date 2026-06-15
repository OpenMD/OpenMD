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

#include "perturbations/VelocityField.hpp"

#include "utils/simError.h"

namespace OpenMD {

  VelocityField::VelocityField(SimInfo* info) : info_(info) {
    vfParams_ = info_->getSimParams()->getVelocityFieldParameters();
    initialize();
  }

  void VelocityField::initialize() {
    if (initialized_) return;
    initialized_ = true;

    if (vfParams_ == NULL || !vfParams_->getUseVelocityField()) {
      doVelocityField_ = false;
      return;
    }

    // ----- which parts were requested? -----
    bool haveSR   = vfParams_->haveStrainRate();
    bool haveSD1  = vfParams_->haveStrainDirection1();
    bool haveSD2  = vfParams_->haveStrainDirection2();
    bool haveVort = vfParams_->haveVorticity();

    bool anyStrain  = haveSR || haveSD1 || haveSD2;
    bool fullStrain = haveSR && haveSD1 && haveSD2;

    if (anyStrain && !fullStrain) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "VelocityField: an incomplete rate-of-strain block was given.\n"
               "\tThe constant-stress part requires all three of strainRate,\n"
               "\tstrainDirection1, and strainDirection2 together (set\n"
               "\tstrainDirection2 = strainDirection1 for uniaxial extension).\n");
      painCave.isFatal = 1;
      simError();
    }

    if (!fullStrain && !haveVort) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "VelocityField: nothing to do.  Specify a rate-of-strain block\n"
               "\t(strainRate, strainDirection1, strainDirection2), a vorticity\n"
               "\tvector, or both.\n");
      painCave.isFatal = 1;
      simError();
    }

    E_ = Mat3x3d(0.0);
    W_ = Mat3x3d(0.0);

    // ----- constant-stress (rate-of-strain) part -----
    // Symmetric, traceless tensor built from two directions and a rate,
    // exactly the UniformGradient construction with g -> strainRate:
    //   E_ij = sr * [ 1/2 (a_i b_j + a_j b_i) - 1/3 (a.b) delta_ij ]
    if (fullStrain) {
      RealType sr = vfParams_->getStrainRate();

      std::vector<RealType> d1 = vfParams_->getStrainDirection1();
      std::vector<RealType> d2 = vfParams_->getStrainDirection2();
      if (d1.size() != 3 || d2.size() != 3) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "VelocityField: strainDirection1 and strainDirection2 must\n"
                 "\teach have 3 components (got %zu and %zu).\n",
                 d1.size(), d2.size());
        painCave.isFatal = 1;
        simError();
      }

      Vector3d a(d1[0], d1[1], d1[2]);
      Vector3d b(d2[0], d2[1], d2[2]);
      a.normalize();
      b.normalize();
      RealType cpsi = dot(a, b);

      E_(0, 0) = sr * (a.x() * b.x() - cpsi / 3.0);
      E_(0, 1) = 0.5 * sr * (a.x() * b.y() + a.y() * b.x());
      E_(0, 2) = 0.5 * sr * (a.x() * b.z() + a.z() * b.x());
      E_(1, 0) = E_(0, 1);
      E_(1, 1) = sr * (a.y() * b.y() - cpsi / 3.0);
      E_(1, 2) = 0.5 * sr * (a.y() * b.z() + a.z() * b.y());
      E_(2, 0) = E_(0, 2);
      E_(2, 1) = E_(1, 2);
      E_(2, 2) = sr * (a.z() * b.z() - cpsi / 3.0);
    }

    // ----- constant-vorticity part -----
    // W = 1/2 [omega]_x, so that W r = 1/2 (omega x r) and curl(W r) = omega.
    if (haveVort) {
      std::vector<RealType> w = vfParams_->getVorticity();
      if (w.size() != 3) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "VelocityField: vorticity must have 3 components (got %zu).\n",
                 w.size());
        painCave.isFatal = 1;
        simError();
      }

      W_(0, 1) = -0.5 * w[2];
      W_(0, 2) = 0.5 * w[1];
      W_(1, 0) = 0.5 * w[2];
      W_(1, 2) = -0.5 * w[0];
      W_(2, 0) = -0.5 * w[1];
      W_(2, 1) = 0.5 * w[0];
    }

    // ----- optional uniform background velocity -----
    if (vfParams_->haveBackgroundVelocity()) {
      std::vector<RealType> v0 = vfParams_->getBackgroundVelocity();
      if (v0.size() != 3) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "VelocityField: backgroundVelocity must have 3 components\n"
                 "\t(got %zu).\n",
                 v0.size());
        painCave.isFatal = 1;
        simError();
      }
      v0_ = Vector3d(v0[0], v0[1], v0[2]);
    }

    // ----- assemble and canonicalize from K -----
    K_ = E_ + W_;

    // Recompute E_, W_, omega_ from K_ as the single source of truth.
    E_         = 0.5 * (K_ + K_.transpose());
    W_         = 0.5 * (K_ - K_.transpose());
    omega_.x() = K_(2, 1) - K_(1, 2);
    omega_.y() = K_(0, 2) - K_(2, 0);
    omega_.z() = K_(1, 0) - K_(0, 1);

    doVelocityField_ = true;
  }
}  // namespace OpenMD
