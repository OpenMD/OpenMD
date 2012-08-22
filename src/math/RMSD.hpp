/*
 *  *******************************************************************
 *
 *  rmsd.h 
 *  (c) 2005 Bosco K Ho
 * 
 *  Implementation of the Kabsch algorithm to find the RMSD, and 
 *  the least-squares rotation matrix for a superposition between 
 *  two sets of vectors.
 *
 *  This implementation is completely self-contained. No other dependencies.
 *
 *  **************************************************************************
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published
 *  by the Free Software Foundation; either version 2.1 of the License, or (at
 *  your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,  but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details. 
 *  
 *  You should have received a copy of the GNU Lesser General Public License 
 *  along with this program; if not, write to the Free Software Foundation, 
 *  Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 *  **************************************************************************
 * 
 */

#ifndef MATH_RMSD_HPP
#define MATH_RMSD_HPP

#include <vector>

#include "config.h"
#include "math/SquareMatrix3.hpp"

namespace OpenMD {

  class RMSD{
  public:
    RMSD();
    RMSD(std::vector<Vector3d> ref) {
      set_reference_structure(ref);
    }
    virtual ~RMSD() { };

    void set_reference_structure(std::vector<Vector3d> ref) {
      ref_ = ref;
      ref_com = V3Zero;
      
      for (unsigned int n = 0; n < ref_.size(); n++) {
        ref_com += ref_[n];
      }
      ref_com /= (RealType)ref.size();      
    }
    
    /*
     * calculate_rmsd()
     *
     *   given a vector of Vector3 coordinates, constructs
     *    - mov_com: the centre of mass of the mov list
     *    - mov_to_ref: vector between the com of mov and ref
     *    - U: the rotation matrix for least-squares, usage of
     *
     *   returns
     *    - rmsd: measures similarity between the vectors
     */
    RealType calculate_rmsd(std::vector<Vector3d> mov,
                            Vector3d mov_com,
                            Vector3d mov_to_ref);
    
    /* 
     * optimal_superposition()
     *
     *   Returns best-fit rotation matrix 
     */
    RotMat3x3d optimal_superposition(std::vector<Vector3d> mov,
                            Vector3d mov_com,
                            Vector3d mov_to_ref);

    
  protected:
    std::vector<Vector3d> ref_;
    Vector3d ref_com;

  };
}
    

#endif
