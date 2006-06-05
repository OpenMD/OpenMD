/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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

#include "hydrodynamics/Ellipsoid.hpp"
#include "utils/OOPSEConstant.hpp"
#include "math/LU.hpp"

namespace oopse {
  
  Ellipsoid::Ellipsoid(Vector3d origin, RealType rMajor, RealType rMinor,Mat3x3d rotMat) 
    : origin_(origin), rMajor_(rMajor), rMinor_(rMinor), rotMat_(rotMat) {
    
  }
  bool Ellipsoid::isInterior(Vector3d pos) {
    Vector3d r = pos - origin_;
    Vector3d rbody = rotMat_ * r;
    RealType xovera = rbody[0]/rMajor_;
    RealType yovera = rbody[1]/rMajor_;
    RealType zoverb = rbody[2]/rMinor_;
    
    bool result;
    if (xovera*xovera + yovera*yovera + zoverb*zoverb < 1)
      result = true;
    else
      result = false;
    
    return result;    
  }
  
  std::pair<Vector3d, Vector3d> Ellipsoid::getBoundingBox() {
    
    std::pair<Vector3d, Vector3d>  boundary;
    //make a cubic box
    RealType rad  = rMajor_ > rMinor_ ? rMajor_ : rMinor_; 
    Vector3d r(rad, rad, rad);
    boundary.first = origin_ - r;
    boundary.second = origin_ + r;
    return boundary;
  }
  
  HydroProp* Ellipsoid::getHydroProp(RealType viscosity, RealType temperature) {
    
    RealType a = rMinor_;
    RealType b = rMajor_;
    RealType a2 = a * a;
    RealType b2 = b* b;
    
    RealType p = a /b;
    RealType S;
    if (p > 1.0) { //prolate
      S = 2.0/sqrt(a2 - b2) * log((a + sqrt(a2-b2))/b);
    } else { //oblate
      S = 2.0/sqrt(b2 - a2) * atan(sqrt(b2-a2)/a);
    }
    
    RealType P = 1.0/(a2 - b2) * (S - 2.0/a);
    RealType Q = 0.5/(a2-b2) * (2.0*a/b2 - S);
    
    RealType transMinor = 16.0 * NumericConstant::PI * viscosity * (a2 - b2) /((2.0*a2-b2)*S -2.0*a);
    RealType transMajor = 32.0 * NumericConstant::PI * viscosity * (a2 - b2) /((2.0*a2-3.0*b2)*S +2.0*a);
    RealType rotMinor = 32.0/3.0 * NumericConstant::PI * viscosity *(a2 - b2) * b2 /(2.0*a -b2*S);
    RealType rotMajor = 32.0/3.0 * NumericConstant::PI * viscosity *(a2*a2 - b2*b2)/((2.0*a2-b2)*S-2.0*a);
    
    
    Mat6x6d Xi, XiCopy, D;
    
    Xi(0,0) = transMajor;
    Xi(1,1) = transMajor;
    Xi(2,2) = transMinor;
    Xi(3,3) = rotMajor;
    Xi(4,4) = rotMajor;
    Xi(5,5) = rotMinor;
    
    const RealType convertConstant = 6.023; //convert poise.angstrom to amu/fs
    Xi *= convertConstant;    
    
    XiCopy = Xi;
    invertMatrix(XiCopy, D);
    RealType kt = OOPSEConstant::kB * temperature;
    D *= kt;
    Xi *= OOPSEConstant::kb * temperature;
   
    HydroProp* hprop = new HydroProp(V3Zero, Xi, D);

    return hprop;

  }  
}
