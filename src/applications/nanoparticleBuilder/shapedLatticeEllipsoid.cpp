/*
 * Copyright (c) 2006 The University of Notre Dame. All Rights Reserved.
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
 * [4]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011). 
 */


#include "lattice/shapedLattice.hpp"
#include "shapedLatticeEllipsoid.hpp"
#include <cmath>

using namespace std;
namespace OpenMD {
  
  shapedLatticeEllipsoid::shapedLatticeEllipsoid(RealType latticeConstant,
                                                 std::string latticeType, 
                                                 RealType rAxial, 
                                                 RealType rEquatorial) : 
    shapedLattice(latticeConstant, latticeType){
    
    rAxial_= rAxial;
    rEquatorial_= rEquatorial;

    Vector3d dimension;
    dimension[0]=2.0*rEquatorial_;
    dimension[1]=2.0*rEquatorial_;
    dimension[2]=2.0*rAxial_;
    setGridDimension(dimension);
    Vector3d origin;
    origin[0] = latticeConstant / 2.0;
    origin[1] = latticeConstant / 2.0;
    origin[2] = latticeConstant / 2.0;
    setOrigin(origin);
  }

  /**
   * Determines whether a point lies within an ellipsoid centered at
   * origin (0,0,0)
   *
   */ 
  bool shapedLatticeEllipsoid::isInterior(Vector3d pos) {
    RealType xoverb = pos[0]/rEquatorial_;
    RealType yoverb = pos[1]/rEquatorial_;
    RealType zovera = pos[2]/rAxial_;
    
    bool result;
    if (xoverb*xoverb + yoverb*yoverb + zovera*zovera < 1)
      result = true;
    else
      result = false;
    
    return result;    
  }    
}
