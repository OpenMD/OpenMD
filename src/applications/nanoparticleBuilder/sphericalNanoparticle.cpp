/*
 *  sphericalNanoparticle.cpp
 *  OOPSE-2.0
 *
 *  Created by Charles F. Vardeman II on 9/28/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
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

#include "sphericalNanoparticle.hpp"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <fstream>

#include "config.h"


#include "lattice/LatticeFactory.hpp"
#include "utils/MoLocator.hpp"
#include "lattice/Lattice.hpp"
#include "brains/Register.hpp"
#include "brains/SimInfo.hpp"
#include "brains/SimCreator.hpp"
#include "io/DumpWriter.hpp"
#include "math/Vector3.hpp"
#include "math/SquareMatrix3.hpp"
#include "utils/StringUtils.hpp"

using namespace std;
using namespace oopse;

sphericalNanoparticle::sphericalNanoparticle(double radius, int latticeConstant){
  
  particleRadius_ = radius;
  particleDiameter_ = 2.0 * radius;
  latticeConstant_ = latticeConstant;
  
  // Number of Unit Cells in Length first
  ny_ = (int)(particleDiameter/latticeConstant);
  // Number of unit cells in Width
  nx_ = (int)(particleDiameter/latticeConstant);
  nz_ = (int)(particleDiameter/latticeConstant);
}

int sphericalNanoparticle::getNMol(Lattice& simpleLat){
  int numMol;
  int numMolPerCell
    std::vector<Vector3d> latticePos;
  
  numMol = 0;
  numMolPerCell = simpleLat->getNumSitesPerCell();
  for(int i = -nx_; i < nx_; i++) {     
    for(int j = -ny_; j < ny_; j++) {       
      for(int k = -nz_; k < nz_; k++) {
        simpleLat->getLatticePointsPos(latticePos, i, j, k);
        
        for(int l = 0; l < numMolPerCell; l++) {
          rx = latticePos[l][0];
          ry = latticePos[l][1];
          rz = latticePos[l][2];
          dist = sqrt(rx*rx + ry*ry + rz*rz)
            if (dist <= particleRadius){
              numMol++;
            }
        }
      }
    }
  }
  
  return numMol;
}

int sphericalNanoparticle::getNMol(Lattice &simpleLat, int nComponents, double &shellRadius,int &numSites){
  int numMol;
  int numMolPerCell
    std::vector<Vector3d> latticePos;
  
  numMol = 0;
  numMolPerCell = simpleLat->getNumSitesPerCell();
  for(int i = -nx_; i < nx_; i++) {     
    for(int j = -ny_; j < ny_; j++) {       
      for(int k = -nz_; k < nz_; k++) {
        simpleLat->getLatticePointsPos(latticePos, i, j, k);
        
        for(int l = 0; l < numMolPerCell; l++) {
          rx = latticePos[l][0];
          ry = latticePos[l][1];
          rz = latticePos[l][2];
          dist = sqrt(rx*rx + ry*ry + rz*rz)
            if (dist <= particleRadius){
              for (int l = 0; l< nComponents -1;l++){
                if (dist <= shellRadius[l]{
                  numMol++;
                }
            }
        }
      }
    }
  }
  
  return numMol;
}




void sphericalNanoparticle::buildSphericalNanoparticle(Molecule* ){
  
  
  
  
  
}

