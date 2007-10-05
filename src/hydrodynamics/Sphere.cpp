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

#include "hydrodynamics/Sphere.hpp"
#include "utils/OOPSEConstant.hpp"
#include "math/LU.hpp"

namespace oopse {
  
  Sphere::Sphere(Vector3d origin, RealType radius) : origin_(origin), radius_(radius){
    
  }
  
  bool Sphere::isInterior(Vector3d pos) {
    Vector3d r = pos - origin_;
    
    bool result;
    if (r.length() < radius_)
      result = true;
    else
      result = false;
    
    return result;
  }
  
  std::pair<Vector3d, Vector3d> Sphere::getBoundingBox() {
    std::pair<Vector3d, Vector3d>  boundary;
    Vector3d r(radius_, radius_, radius_);
    boundary.first = origin_ - r;
    boundary.second = origin_ + r;
    return boundary;
  }
  
  HydroProp* Sphere::getHydroProp(RealType viscosity, RealType temperature) {
    
    RealType Xitt  = 6.0 * NumericConstant::PI * viscosity * radius_;
    RealType Xirr = 8.0 * NumericConstant::PI * viscosity * radius_ * radius_ * radius_;

    Mat6x6d Xi, XiCopy, D;

    Xi(0, 0) = Xitt;
    Xi(1, 1) = Xitt;
    Xi(2, 2) = Xitt;
    Xi(3, 3) = Xirr;
    Xi(4, 4) = Xirr;
    Xi(5, 5) = Xirr;
    
    const RealType convertConstant = 1.439326479e4; //converts Poise angstroms
                                                    // to kcal fs mol^-1 Angstrom^-1
    Xi *= convertConstant;
    XiCopy = Xi;

    invertMatrix(XiCopy, D);
    RealType kt = OOPSEConstant::kb * temperature; // in kcal mol^-1
    D *= kt;  // now in angstroms^2 fs^-1  (at least for Trans-trans)

    HydroProp* hprop = new HydroProp(V3Zero, Xi, D);
    
    return hprop;
  }
  
}
