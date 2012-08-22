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
 
#ifndef LATTICE_BASELATTICE_HPP
#define LATTICE_BASELATTICE_HPP

#include <vector>
#include "math/Vector3.hpp"

namespace OpenMD {

  class Lattice{
  protected:
    Lattice(){
      
      setOrigin(V3Zero);
    }
    
  public:

    //virtual destructor of Lattice
    virtual ~Lattice() {}
    
    int getNumSitesPerCell() {return nCellSites;}

    void getLatticePointsPos(std::vector<Vector3d>& latticePos, int nx, int ny, int nz);

    std::vector<Vector3d> getLatticePointsOrt() {return cellSitesOrt;}
    
    //get lattice constant of unit cell
    virtual  std::vector<RealType> getLatticeConstant() =0;

    //set lattice constant of unit cell
    virtual void setLatticeConstant(const  std::vector<RealType>& lc)=0;

    //get origin of unit cell
    Vector3d getOrigin( ) {return origin;} 

    //set origin of unit cell
    void setOrigin(const Vector3d& newOrigin){
      this->origin = newOrigin;
    }

  protected:
    virtual void update() =0;
    
    int nCellSites;
    Vector3d origin;    
    std::vector<Vector3d> cellSitesPos;
    std::vector<Vector3d> cellSitesOrt;
    Vector3d cellLen;
  };


}

#endif
