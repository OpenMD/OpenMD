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
 
#include "lattice/FCCLattice.hpp"

namespace oopse {

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

