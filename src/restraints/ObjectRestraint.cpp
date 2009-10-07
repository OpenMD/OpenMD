/*
 * Copyright (c) 2009 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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
 */

#include "restraints/ObjectRestraint.hpp"

namespace oopse {

  void ObjectRestraint::calcForce(Vector3d struc) {

    pot_ = 0.0;
    
    if (restType_ & rtDisplacement) {
      Vector3d del = struc - refPos_;
      RealType r = del.length();  
      Vector3d frc = -kDisp_ * del;
      RealType p = 0.5 * kDisp_ * del.lengthSquare();
      pot_ += p;
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

      RealType twistAngle, swingAngle;
      Vector3d swingAxis;
      RealType tw, swingX, swingY;
      
      quat.getTwistSwingAxisAngle(twistAngle, swingAngle, swingAxis);
      quat.toSwingTwist(tw, swingX, swingY);
      
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
