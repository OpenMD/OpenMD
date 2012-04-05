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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [4]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011). *
 *  Created by Kelsey M. Stocker on 2/9/12.
 *  @author  Kelsey M. Stocker 
 *  @version $Id: shapedLatticeRod.cpp 1665 2011-11-22 20:38:56Z gezelter $
 *
 */

#include "lattice/shapedLattice.hpp"
#include "shapedLatticePentRod.hpp"
#include "shapedLatticeRod.hpp"
#include <cmath>

using namespace std;
namespace OpenMD {
  
  shapedLatticePentRod::shapedLatticePentRod(RealType latticeConstant,
				     std::string latticeType, 
				     RealType radius, RealType length) : shapedLattice(latticeConstant, latticeType){
    
    rodRadius_= radius;
    rodLength_= length;
    double x, y, z, new_x, new_y, new_z, new_x72, new_z72, left_newx, right_newx, left_newx72, right_newx72;
    int z_int;					
    Vector3d dimension;
    dimension[0] = 2.0*length;
    dimension[1] = 2.0*length;
    dimension[2] = 2.0*length;
    cerr << "using dimension = " << dimension << "\n";
    setGridDimension(dimension);
    cerr << "done!\n";
    Vector3d origin;
    origin[0] = 0;
    origin[1] = 0;
    origin[2] = 2.04;
    setOrigin(origin);
  }


  /**
   * Creates a wedge for pentagonal nanorods
   *
   */ 

  bool shapedLatticePentRod::isInterior(Vector3d point){

    RealType x, y, z, new_x, new_y, left_newx, right_newx;
    int z_int;

    bool isIT=false;

    x = point[0];
    y = point[1];
    z = point[2];

    z_int = z/2.04;

    //Rotate by 45 degrees around z-axis so length of rod lies along y axis
    new_x = (sqrt(2)/2)*(x - y);
    new_y = (sqrt(2)/2)*(x + y);
 
    left_newx = (z - 1.44)*(0.577350269/0.816496581);
    right_newx = (z + 1.44)*(-0.577350269/0.816496581);

    //Make center spine of nanorod down new_y axis
    //This is now done directly in nanorod_pentBuilder.cpp
    /*if ( (new_x == 0) && (z == 0) ) {

      if ( abs(new_y) <= 0.5*rodLength_ + 0.5773502692*rodRadius_ ) {

	isIT = true;

      }  
      }*/

    //Make one wedge
    if ( (z < 0) && (z >= -0.816496581*rodRadius_ - 1.44) ) {
      
      if ( abs(new_y) <= 1.44*(z/2.04) + 0.5*rodLength_ + 0.5773502692*rodRadius_ ) { 
	
	if ( (new_x >= left_newx) && (new_x <= right_newx) ) {
	  
	  isIT=true;
	  
	}	  
      }
    }	   
    return isIT;
  }
}
