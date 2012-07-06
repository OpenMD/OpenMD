/*
 * Copyright (c) 2009 The University of Notre Dame. All Rights Reserved.
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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */

#include "restraints/ObjectRestraint.hpp"

namespace OpenMD {

  void ObjectRestraint::calcForce(Vector3d struc) {

    pot_ = 0.0;
    if (restType_ & rtDisplacement) {
      Vector3d del = struc - refPos_;
      RealType r = del.length();
      Vector3d frc = -kDisp_ * del;
      RealType p = 0.5 * kDisp_ * del.lengthSquare();
      pot_ = p;
      force_ = frc * scaleFactor_;
      restInfo_[rtDisplacement] = std::make_pair(r,p);
    }
  }
    
  void ObjectRestraint::calcForce(Vector3d struc, RotMat3x3d A) {

    calcForce(struc);

    // rtDisplacement is 1, so anything higher than that requires orientations:
    if (restType_ > 1) {
      
      Vector3d tBody(0.0);
      
      RotMat3x3d temp = A * refA_.transpose();

      Quat4d quat = temp.toQuaternion();

      RealType twistAngle;
      Vector3d swingAxis;
      RealType swingX, swingY;
      
      quat.toSwingTwist(swingX, swingY, twistAngle);

      RealType dVdtwist, dVdswingX, dVdswingY;
      RealType dTwist, dSwingX, dSwingY;
      RealType p;
      Vector3d tTwist, tSwing;

      if (restType_ & rtTwist){
        dTwist = twistAngle - twist0_;
        dVdtwist = kTwist_ * sin(dTwist);
        p = kTwist_ * (1.0 - cos(dTwist) );
        pot_ += p;
        tBody -= dVdtwist * V3Z;
        restInfo_[rtTwist] = std::make_pair(twistAngle, p);
      }

      if (restType_ & rtSwingX){
        dSwingX = swingX - swingX0_;
        dVdswingX = kSwingX_ * 0.5 * sin(2.0 * dSwingX);
        p = 0.25 * kSwingX_ * (1.0 - cos(2.0 * dSwingX));
        pot_ += p;
        tBody -= dVdswingX * V3X;
        restInfo_[rtSwingX] = std::make_pair(swingX, p);
      }

      if (restType_ & rtSwingY){
        dSwingY = swingY - swingY0_;
        dVdswingY = kSwingY_ * 0.5 * sin(2.0 * dSwingY);
        p = 0.25 * kSwingY_ * (1.0 - cos(2.0 * dSwingY));
        pot_ += p;
        tBody -= dVdswingY * V3Y;
        restInfo_[rtSwingY] = std::make_pair(swingY, p);
      }

      Vector3d tLab = A.transpose() * tBody;      
      torque_ = tLab * scaleFactor_;      
    }
  }
}
