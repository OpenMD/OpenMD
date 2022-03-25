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

#include "restraints/MolecularRestraint.hpp"

#include <utility>

#include "math/SVD.hpp"
#include "math/SquareMatrix3.hpp"

// using namespace JAMA;

namespace OpenMD {

  void MolecularRestraint::calcForce(std::vector<Vector3d> struc,
                                     Vector3d molCom) {
    assert(struc.size() == ref_.size());

    std::vector<Vector3d>::iterator it;

    // clear out initial values:
    pot_ = 0.0;
    for (it = forces_.begin(); it != forces_.end(); ++it)
      (*it) = 0.0;

    if (restType_ & rtDisplacement) {
      Vector3d del = molCom - refCom_;

      RealType r = del.length();
      RealType p = 0.5 * kDisp_ * r * r;

      pot_ += p;

      if (printRest_) restInfo_[rtDisplacement] = std::make_pair(r, p);

      for (it = forces_.begin(); it != forces_.end(); ++it)
        (*it) += -kDisp_ * del * scaleFactor_;
    }

    if (restType_ & rtAbsoluteZ) {
      RealType r   = molCom(2) - posZ0_;
      RealType p   = 0.5 * kAbs_ * r * r;
      Vector3d frc = Vector3d(0.0, 0.0, -kAbs_ * r);

      pot_ += p;

      if (printRest_) restInfo_[rtAbsoluteZ] = std::make_pair(r, p);

      for (it = forces_.begin(); it != forces_.end(); ++it)
        (*it) += frc * scaleFactor_;
    }

    // rtDisplacement is 1, rtAbsolute is 2, so anything higher than 3
    // that requires orientations:
    if (restType_ > 3) {
      Vector3d tBody(0.0);

      Mat3x3d R(0.0);

      for (unsigned int n = 0; n < struc.size(); n++) {
        /*
         * First migrate the center of mass:
         */
        struc[n] -= molCom;

        /*
         * correlation matrix R:
         *   R(i,j) = sum(over n): y(n,i) * x(n,j)
         *   where x(n) and y(n) are two vector sets
         */

        R += outProduct(struc[n], ref_[n]);
      }

      // SVD class uses dynamic matrices, so we must wrap the correlation
      // matrix before calling SVD and then unwrap the results into Mat3x3d
      // and Vector3d before we use them.

      DynamicRectMatrix<RealType> Rtmp(3, 3, 0.0);
      DynamicRectMatrix<RealType> vtmp(3, 3);
      DynamicVector<RealType> stmp(3);
      DynamicRectMatrix<RealType> wtmp(3, 3);

      Rtmp.setSubMatrix(0, 0, R);

      // Heavy lifting goes here:

      JAMA::SVD<RealType> svd(Rtmp);

      svd.getU(vtmp);
      svd.getSingularValues(stmp);
      svd.getV(wtmp);

      Mat3x3d v;
      Vector3d s;
      Mat3x3d w_tr;

      vtmp.getSubMatrix(0, 0, v);
      stmp.getSubVector(0, s);
      wtmp.getSubMatrix(0, 0, w_tr);

      bool is_reflection = (v.determinant() * w_tr.determinant()) < 0.0;

      if (is_reflection) {
        v(2, 0) = -v(2, 0);
        v(2, 1) = -v(2, 1);
        v(2, 2) = -v(2, 2);
      }

      RotMat3x3d Atrans = v * w_tr.transpose();
      RotMat3x3d A      = Atrans.transpose();

      Quat4d quat = A.toQuaternion();

      RealType twistAngle;
      RealType swingX, swingY;

      quat.toTwistSwing(twistAngle, swingX, swingY);

      RealType dVdtwist, dVdswingX, dVdswingY;
      RealType dTwist, dSwingX, dSwingY;
      RealType p;

      if (restType_ & rtTwist) {
        dTwist = twistAngle - twist0_;
        /// dVdtwist = kTwist_ * sin(dTwist) ;
        /// p = kTwist_ * (1.0 - cos(dTwist) ) ;
        dVdtwist = kTwist_ * dTwist;
        p        = 0.5 * kTwist_ * dTwist * dTwist;
        pot_ += p;
        tBody -= dVdtwist * V3Z;
        if (printRest_) restInfo_[rtTwist] = std::make_pair(twistAngle, p);
      }

      if (restType_ & rtSwingX) {
        dSwingX = swingX - swingX0_;
        /// dVdswingX = kSwingX_ * 2.0 * sin(2.0 * dSwingX);
        /// p = kSwingX_ * (1.0 - cos(2.0 * dSwingX));
        dVdswingX = kSwingX_ * dSwingX;
        p         = 0.5 * kSwingX_ * dSwingX * dSwingX;
        pot_ += p;
        tBody -= dVdswingX * V3X;
        if (printRest_) restInfo_[rtSwingX] = std::make_pair(swingX, p);
      }
      if (restType_ & rtSwingY) {
        dSwingY = swingY - swingY0_;
        /// dVdswingY = kSwingY_ * 2.0 * sin(2.0 * dSwingY);
        /// p = kSwingY_ * (1.0 - cos(2.0 * dSwingY));
        dVdswingY = kSwingY_ * dSwingY;
        p         = 0.5 * kSwingX_ * dSwingY * dSwingY;
        pot_ += p;
        tBody -= dVdswingY * V3Y;
        if (printRest_) restInfo_[rtSwingY] = std::make_pair(swingY, p);
      }

      RealType t2 = dot(tBody, tBody);

      Vector3d rLab, rBody, txr, fBody, fLab;

      for (unsigned int i = 0; i < struc.size(); i++) {
        rLab  = struc[i];
        rBody = A * rLab;

        txr   = cross(tBody, rBody);
        fBody = txr * t2;
        fLab  = Atrans * fBody;
        fLab *= scaleFactor_;

        forces_[i] += fLab;
      }

      // test the force vectors and see if it is the right orientation
      //       std::cout << struc.size() << std::endl << std::endl;
      //       for (int i = 0; i != struc.size(); ++i){
      //         std::cout << "H\t" << struc[i].x() << "\t" << struc[i].y() <<
      //         "\t" << struc[i].z() << "\t"; std::cout << forces_[i].x() <<
      //         "\t"
      //         << forces_[i].y()
      //         << "\t" << forces_[i].z() << std::endl;
      //       }
    }
  }
}  // namespace OpenMD
