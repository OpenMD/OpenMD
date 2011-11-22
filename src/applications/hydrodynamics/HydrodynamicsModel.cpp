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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
#include "applications/hydrodynamics/HydrodynamicsModel.hpp"
#include "hydrodynamics/Shape.hpp"
#include "hydrodynamics/Sphere.hpp"
#include "hydrodynamics/Ellipsoid.hpp"
#include "applications/hydrodynamics/CompositeShape.hpp"

namespace OpenMD {
  
  bool HydrodynamicsModel::calcHydroProps(Shape* shape, RealType viscosity, RealType temperature) {
    return false;
  }
  
  void HydrodynamicsModel::writeHydroProps(std::ostream& os) {

    Vector3d center;
    Mat6x6d Xi, D;
    
    os << sd_->getType() << "\t";
    
    //center of resistance

    center = cr_->getCOR();

    os << center[0] <<  "\t" << center[1] <<  "\t" << center[2] <<  "\t"; 
    
    //resistance tensor at center of resistance
    //translation

    Xi = cr_->getXi();

    os << Xi(0, 0) <<  "\t" << Xi(0, 1) <<  "\t" << Xi(0, 2) <<  "\t" 
       << Xi(1, 0) <<  "\t" << Xi(1, 1) <<  "\t" << Xi(1, 2) <<  "\t" 
       << Xi(2, 0) <<  "\t" << Xi(2, 1) <<  "\t" << Xi(2, 2) <<  "\t"; 
    
    //rotation-translation
    os << Xi(0, 3) <<  "\t" << Xi(0, 4) <<  "\t" << Xi(0, 5) <<  "\t" 
       << Xi(1, 3) <<  "\t" << Xi(1, 4) <<  "\t" << Xi(1, 5) <<  "\t" 
       << Xi(2, 3) <<  "\t" << Xi(2, 4) <<  "\t" << Xi(2, 5) <<  "\t";
    
    //translation-rotation
    os << Xi(3, 0) <<  "\t" << Xi(3, 1) <<  "\t" << Xi(3, 2) <<  "\t" 
       << Xi(4, 0) <<  "\t" << Xi(4, 1) <<  "\t" << Xi(4, 2) <<  "\t" 
       << Xi(5, 0) <<  "\t" << Xi(5, 1) <<  "\t" << Xi(5, 2) <<  "\t"; 
    
    //rotation
    os << Xi(3, 3) <<  "\t" << Xi(3, 4) <<  "\t" << Xi(3, 5) <<  "\t" 
       << Xi(4, 3) <<  "\t" << Xi(4, 4) <<  "\t" << Xi(4, 5) <<  "\t" 
       << Xi(5, 3) <<  "\t" << Xi(5, 4) <<  "\t" << Xi(5, 5) <<  "\t"; 
    
    
    //diffusion tensor at center of resistance
    //translation

    D = cr_->getD();

    os << D(0, 0) <<  "\t" << D(0, 1) <<  "\t" << D(0, 2) <<  "\t" 
       << D(1, 0) <<  "\t" << D(1, 1) <<  "\t" << D(1, 2) <<  "\t" 
       << D(2, 0) <<  "\t" << D(2, 1) <<  "\t" << D(2, 2) <<  "\t"; 
    
    //rotation-translation
    os << D(0, 3) <<  "\t" << D(0, 4) <<  "\t" << D(0, 5) <<  "\t" 
       << D(1, 3) <<  "\t" << D(1, 4) <<  "\t" << D(1, 5) <<  "\t" 
       << D(2, 3) <<  "\t" << D(2, 4) <<  "\t" << D(2, 5) <<  "\t"; 
    
    //translation-rotation
    os << D(3, 0) <<  "\t" << D(3, 1) <<  "\t" << D(3, 2) <<  "\t" 
       << D(4, 0) <<  "\t" << D(4, 1) <<  "\t" << D(4, 2) <<  "\t" 
       << D(5, 0) <<  "\t" << D(5, 1) <<  "\t" << D(5, 2) <<  "\t"; 
    
    //rotation
    os << D(3, 3) <<  "\t" << D(3, 4) <<  "\t" << D(3, 5) <<  "\t" 
       << D(4, 3) <<  "\t" << D(4, 4) <<  "\t" << D(4, 5) <<  "\t" 
       << D(5, 3) <<  "\t" << D(5, 4) <<  "\t" << D(5, 5) <<  "\t"; 
    
    //---------------------------------------------------------------------
    
    //center of diffusion

    center = cd_->getCOR();

    os << center[0] <<  "\t" << center[1] <<  "\t" << center[2] <<  "\t"; 
    
    //resistance tensor at center of diffusion
    //translation

    Xi = cd_->getXi();

    os << Xi(0, 0) <<  "\t" << Xi(0, 1) <<  "\t" << Xi(0, 2) <<  "\t" 
       << Xi(1, 0) <<  "\t" << Xi(1, 1) <<  "\t" << Xi(1, 2) <<  "\t" 
       << Xi(2, 0) <<  "\t" << Xi(2, 1) <<  "\t" << Xi(2, 2) <<  "\t"; 
    
    //rotation-translation
    os << Xi(0, 3) <<  "\t" << Xi(0, 4) <<  "\t" << Xi(0, 5) <<  "\t" 
       << Xi(1, 3) <<  "\t" << Xi(1, 4) <<  "\t" << Xi(1, 5) <<  "\t" 
       << Xi(2, 3) <<  "\t" << Xi(2, 4) <<  "\t" << Xi(2, 5) <<  "\t"; 
    
    //translation-rotation
    os << Xi(3, 0) <<  "\t" << Xi(3, 1) <<  "\t" << Xi(3, 2) <<  "\t" 
       << Xi(4, 0) <<  "\t" << Xi(4, 1) <<  "\t" << Xi(4, 2) <<  "\t" 
       << Xi(5, 0) <<  "\t" << Xi(5, 1) <<  "\t" << Xi(5, 2) <<  "\t"; 
    
    //rotation
    os << Xi(3, 3) <<  "\t" << Xi(3, 4) <<  "\t" << Xi(3, 5) <<  "\t" 
       << Xi(4, 3) <<  "\t" << Xi(4, 4) <<  "\t" << Xi(4, 5) <<  "\t" 
       << Xi(5, 3) <<  "\t" << Xi(5, 4) <<  "\t" << Xi(5, 5) <<  "\t"; 
    
    
    //diffusion tensor at center of diffusion
    //translation

    D = cd_->getD();

    os << D(0, 0) <<  "\t" << D(0, 1) <<  "\t" << D(0, 2) <<  "\t" 
       << D(1, 0) <<  "\t" << D(1, 1) <<  "\t" << D(1, 2) <<  "\t" 
       << D(2, 0) <<  "\t" << D(2, 1) <<  "\t" << D(2, 2) <<  "\t"; 
    
    //rotation-translation
    os << D(0, 3) <<  "\t" << D(0, 4) <<  "\t" << D(0, 5) <<  "\t" 
       << D(1, 3) <<  "\t" << D(1, 4) <<  "\t" << D(1, 5) <<  "\t" 
       << D(2, 3) <<  "\t" << D(2, 4) <<  "\t" << D(2, 5) <<  "\t"; 
    
    //translation-rotation
    os << D(3, 0) <<  "\t" << D(3, 1) <<  "\t" << D(3, 2) <<  "\t" 
       << D(4, 0) <<  "\t" << D(4, 1) <<  "\t" << D(4, 2) <<  "\t" 
       << D(5, 0) <<  "\t" << D(5, 1) <<  "\t" << D(5, 2) <<  "\t"; 
    
    //rotation
    os << D(3, 3) <<  "\t" << D(3, 4) <<  "\t" << D(3, 5) <<  "\t" 
       << D(4, 3) <<  "\t" << D(4, 4) <<  "\t" << D(4, 5) <<  "\t" 
       << D(5, 3) <<  "\t" << D(5, 4) <<  "\t" << D(5, 5) <<  "\n";    
    
  }
  
}
