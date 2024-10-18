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

#include "restraints/ObjectRestraint.hpp"

namespace OpenMD {

  void ObjectRestraint::calcForce(Vector3d struc) {
    pot_   = 0.0;
    force_ = V3Zero;

    if (restType_ & rtDisplacement) {
      Vector3d del = struc - refPos_;
      RealType r   = del.length();
      Vector3d frc = -kDisp_ * del;
      RealType p   = 0.5 * kDisp_ * del.lengthSquare();

      pot_ += p;
      force_ += frc * scaleFactor_;
      if (printRest_) restInfo_[rtDisplacement] = std::make_pair(r, p);
    }

    if (restType_ & rtAbsoluteZ) {
      RealType r   = struc(2) - posZ0_;
      Vector3d frc = Vector3d(0.0, 0.0, -kAbs_ * r);
      RealType p   = 0.5 * kAbs_ * r * r;

      pot_ += p;
      force_ += frc * scaleFactor_;
      if (printRest_) restInfo_[rtAbsoluteZ] = std::make_pair(r, p);
    }
  }

  void ObjectRestraint::calcForce(Vector3d struc, RotMat3x3d A) {
    calcForce(struc);

    // rtDisplacement is 1, rtAbsolute is 2, so anything higher than 3
    // that requires orientations:
    if (restType_ > 3) {
      Vector3d tBody(0.0);

      RotMat3x3d temp = A * refA_.transpose();
      Quat4d quat     = temp.toQuaternion();

      RealType twistAngle;
      RealType swingX, swingY;

      quat.toTwistSwing(twistAngle, swingX, swingY);

      RealType p;
      Vector3d tTwist, tSwing;

      if (restType_ & rtTwist) {
        RealType dTwist = twistAngle - twist0_;
        /// RealType dVdtwist = kTwist_ * sin(dTwist);
        /// p = kTwist_ * (1.0 - cos(dTwist) );
        RealType dVdtwist = kTwist_ * dTwist;
        p                 = 0.5 * kTwist_ * dTwist * dTwist;
        pot_ += p;
        tBody -= dVdtwist * V3Z;

        if (printRest_) restInfo_[rtTwist] = std::make_pair(twistAngle, p);
      }

      if (restType_ & rtSwingX) {
        RealType dSwingX = swingX - swingX0_;
        /// RealType dVdswingX = kSwingX_ * 0.5 * sin(2.0 * dSwingX);
        /// p = 0.25 * kSwingX_ * (1.0 - cos(2.0 * dSwingX));
        RealType dVdswingX = kSwingX_ * dSwingX;
        p                  = 0.5 * kSwingX_ * dSwingX * dSwingX;
        pot_ += p;
        tBody -= dVdswingX * V3X;
        if (printRest_) restInfo_[rtSwingX] = std::make_pair(swingX, p);
      }

      if (restType_ & rtSwingY) {
        RealType dSwingY = swingY - swingY0_;
        /// RealType dVdswingY = kSwingY_ * 0.5 * sin(2.0 * dSwingY);
        /// p = 0.25 * kSwingY_ * (1.0 - cos(2.0 * dSwingY));
        RealType dVdswingY = kSwingY_ * dSwingY;
        p                  = 0.5 * kSwingX_ * dSwingY * dSwingY;
        pot_ += p;
        tBody -= dVdswingY * V3Y;
        if (printRest_) restInfo_[rtSwingY] = std::make_pair(swingY, p);
      }

      Vector3d tLab = A.transpose() * tBody;
      torque_       = tLab * scaleFactor_;
    }
  }
}  // namespace OpenMD
