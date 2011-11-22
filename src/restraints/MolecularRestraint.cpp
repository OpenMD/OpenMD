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
 
#include "restraints/MolecularRestraint.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/SVD.hpp"
#include <utility>

//using namespace JAMA;

namespace OpenMD {


  void MolecularRestraint::calcForce(std::vector<Vector3d> struc, 
                                     Vector3d molCom){

    assert(struc.size() == ref_.size());

    std::vector<Vector3d>::iterator it;

    // clear out initial values:
    pot_ = 0.0;
    for(it = forces_.begin(); it != forces_.end(); ++it)
      (*it) = 0.0;

   
    if (restType_ & rtDisplacement) {
      Vector3d del = molCom - refCom_;     
      
      RealType r = del.length();
      RealType p = 0.5 * kDisp_ * r * r;

      pot_ += p;

      restInfo_[rtDisplacement] = std::make_pair(r, p);

      for(it = forces_.begin(); it != forces_.end(); ++it)
        (*it) = -kDisp_ * del * scaleFactor_;
    }

    for(it = struc.begin(); it != struc.end(); ++it)
      (*it) -= molCom;

    // rtDisplacement = 1, so anything higher than that requires orientations:
    if (restType_ > 1) {
      Vector3d tBody(0.0);
      
      Mat3x3d R(0.0);
      
      for (int n = 0; n < struc.size(); n++){

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
      DynamicVector<RealType>     stmp(3);
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
      
      if (is_reflection){        
        v(2, 0) = -v(2, 0);
        v(2, 1) = -v(2, 1);
        v(2, 2) = -v(2, 2);
      }
      
      RotMat3x3d Atrans = v * w_tr.transpose();
      RotMat3x3d A = Atrans.transpose();

      Vector3d eularAngles = A.toEulerAngles();


      RealType twistAngle, swingAngle;
      Vector3d swingAxis;

      Quat4d quat = A.toQuaternion();  

      RealType tw, sx, sy, ttw, swingX, swingY;
      quat.toSwingTwist(swingX, swingY, twistAngle);

      RealType dVdtwist, dVdswing, dVdswingX, dVdswingY;
      RealType dTwist, dSwing, dSwingX, dSwingY;
      RealType p;

      if (restType_ & rtTwist){
        dTwist = twistAngle - twist0_;
        dVdtwist = kTwist_ * sin(dTwist) ;
        p = kTwist_ * (1.0 - cos(dTwist) ) ;
        pot_ += p;
        tBody -= dVdtwist * V3Z;
        restInfo_[rtTwist] = std::make_pair(twistAngle, p);
      }

//       if (restType_ & rtSwing){
//         dSwing = swingAngle - swing0_;
//         dVdswing = kSwing_ * 2.0 * sin(2.0 * dSwing);
//         p = kSwing_ * (1.0 - cos(2.0 * dSwing));
//         pot_ += p;
//         tBody -= dVdswing * swingAxis;
//         restInfo_[rtSwing] = std::make_pair(swingAngle, p);
//       }

      if (restType_ & rtSwingX){
        dSwingX = swingX - swingX0_;
        dVdswingX = kSwingX_ * 2.0 * sin(2.0 * dSwingX);
        p = kSwingX_ * (1.0 - cos(2.0 * dSwingX));
        pot_ += p;
        tBody -= dVdswingX * V3X;
        restInfo_[rtSwingX] = std::make_pair(swingX, p);
      }
      if (restType_ & rtSwingY){
        dSwingY = swingY - swingY0_;
        dVdswingY = kSwingY_ * 2.0 * sin(2.0 * dSwingY);
        p = kSwingY_ * (1.0 - cos(2.0 * dSwingY));
        pot_ += p;
        tBody -= dVdswingY * V3Y;
        restInfo_[rtSwingY] = std::make_pair(swingY, p);
      }

            
      RealType t2 = dot(tBody, tBody);     
      
      Vector3d rLab, rBody, txr, fBody, fLab;

      for (int i = 0; i < struc.size(); i++) {
                   
        rLab = struc[i];        
        rBody = A * rLab;
        
        txr = cross(tBody, rBody);
        fBody = txr * t2;
        fLab = Atrans * fBody;        
        fLab *= scaleFactor_;
          
        forces_[i] += fLab;
      }

      // test the force vectors and see if it is the right orientation
//       std::cout << struc.size() << std::endl << std::endl;
//       for (int i = 0; i != struc.size(); ++i){
//         std::cout << "H\t" << struc[i].x() << "\t" << struc[i].y() << "\t" << struc[i].z() << "\t";
//         std::cout << forces_[i].x() << "\t" << forces_[i].y() << "\t" << forces_[i].z() << std::endl;
//       }
    }
  }
}
