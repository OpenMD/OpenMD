/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */

#include "hydrodynamics/Ellipsoid.hpp"
#include "utils/PhysicalConstants.hpp"
#include "math/LU.hpp"

namespace OpenMD {
  
  Ellipsoid::Ellipsoid(Vector3d origin, RealType rAxial, RealType rEquatorial,
		       Mat3x3d rotMat) : origin_(origin), rAxial_(rAxial), 
					 rEquatorial_(rEquatorial), 
					 rotMat_(rotMat) {
    if (rAxial_ > rEquatorial_) {
      rMajor_ = rAxial_;
      rMinor_ = rEquatorial_;
    } else {
      rMajor_ = rEquatorial_;
      rMinor_ = rAxial_;
    }          
  }

  bool Ellipsoid::isInterior(Vector3d pos) {
    Vector3d r = pos - origin_;
    Vector3d rbody = rotMat_ * r;

    RealType xoverb = rbody[0]/rEquatorial_;
    RealType yoverb = rbody[1]/rEquatorial_;
    RealType zovera = rbody[2]/rAxial_;
    
    bool result;
    if (xoverb*xoverb + yoverb*yoverb + zovera*zovera < 1)
      result = true;
    else
      result = false;
    
    return result;    
  }
  
  std::pair<Vector3d, Vector3d> Ellipsoid::getBoundingBox() {
    
    std::pair<Vector3d, Vector3d>  boundary;
    //make a cubic box
    RealType rad  = rAxial_ > rEquatorial_ ? rAxial_ : rEquatorial_; 
    Vector3d r(rad, rad, rad);
    boundary.first = origin_ - r;
    boundary.second = origin_ + r;
    return boundary;
  }
  
  HydroProp* Ellipsoid::getHydroProp(RealType viscosity, 
				     RealType temperature) {
    
    RealType a = rAxial_;
    RealType b = rEquatorial_;
    RealType a2 = a * a;
    RealType b2 = b * b;
    
    RealType p = a / b;
    RealType S;
    if (p > 1.0) {  
      // Ellipsoid is prolate:
      S = 2.0/sqrt(a2 - b2) * log((a + sqrt(a2-b2))/b);
    } else { 
      // Ellipsoid is oblate:
      S = 2.0/sqrt(b2 - a2) * atan(sqrt(b2-a2)/a);
    }
    
    RealType pi = NumericConstant::PI;
    RealType XittA = 16.0 * pi * viscosity * (a2 - b2) /((2.0*a2-b2)*S -2.0*a);
    RealType XittB = 32.0 * pi * viscosity * (a2 - b2) /((2.0*a2-3.0*b2)*S +2.0*a);
    RealType XirrA = 32.0/3.0 * pi * viscosity *(a2 - b2) * b2 /(2.0*a -b2*S);
    RealType XirrB = 32.0/3.0 * pi * viscosity *(a2*a2 - b2*b2)/((2.0*a2-b2)*S-2.0*a);
    
    
    Mat6x6d Xi, XiCopy, D;
    
    Xi(0,0) = XittB;
    Xi(1,1) = XittB;
    Xi(2,2) = XittA;
    Xi(3,3) = XirrB;
    Xi(4,4) = XirrB;
    Xi(5,5) = XirrA;

    Xi *= PhysicalConstants::viscoConvert;    
    
    XiCopy = Xi;
    invertMatrix(XiCopy, D);
    RealType kt = PhysicalConstants::kb * temperature; // in kcal mol^-1
    D *= kt;
   
    HydroProp* hprop = new HydroProp(V3Zero, Xi, D);

    return hprop;

  }  
}
