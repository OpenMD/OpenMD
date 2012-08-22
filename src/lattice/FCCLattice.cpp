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
 
#include "lattice/FCCLattice.hpp"

namespace OpenMD {

  FCCLattice::FCCLattice() : CubicLattice(){
    nCellSites = 4;
    cellSitesPos.resize(nCellSites);
    cellSitesOrt.resize(nCellSites);
    update();

  }

  void FCCLattice::update(){

    RealType cellLenOver2;
    RealType oneOverRoot3;

    cellLenOver2  = 0.5 * latticeParam;
    oneOverRoot3 = 1.0 / sqrt(3.0);

    // Molecule 1
    cellSitesPos[0][0] = 0.0; 
    cellSitesPos[0][1] = 0.0;
    cellSitesPos[0][2] = 0.0;
  
    cellSitesOrt[0][0] = oneOverRoot3;
    cellSitesOrt[0][1] = oneOverRoot3;
    cellSitesOrt[0][2] = oneOverRoot3;

    // Molecule 2  
    cellSitesPos[1][0]   = 0.0;
    cellSitesPos[1][1]   = cellLenOver2;
    cellSitesPos[1][2]   = cellLenOver2;

    cellSitesOrt[1][0] = -oneOverRoot3;
    cellSitesOrt[1][1] = oneOverRoot3;
    cellSitesOrt[1][2] = -oneOverRoot3;
   
    // Molecule 3
    cellSitesPos[2][0]   = cellLenOver2;
    cellSitesPos[2][1]   = cellLenOver2;
    cellSitesPos[2][2]   = 0.0;

    cellSitesOrt[2][0] = oneOverRoot3;
    cellSitesOrt[2][1] = -oneOverRoot3;
    cellSitesOrt[2][2] = -oneOverRoot3;

    // Molecule 4

    cellSitesPos[3][0]   = cellLenOver2;
    cellSitesPos[3][1]   = 0.0;
    cellSitesPos[3][2]   = cellLenOver2;

    cellSitesOrt[3][0] = -oneOverRoot3;
    cellSitesOrt[3][1] = oneOverRoot3;
    cellSitesOrt[3][2] = oneOverRoot3;
  }

}

