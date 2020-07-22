/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

#include "applications/hydrodynamics/RoughShell.hpp"
#include "applications/hydrodynamics/ShapeBuilder.hpp"
#include "brains/SimInfo.hpp"

namespace OpenMD {

  RoughShell::RoughShell(StuntDouble* sd, SimInfo* info) :
    ApproximationModel(sd, info) {

    shape_ = ShapeBuilder::createShape(sd);

    Globals* simParams = info->getSimParams();
    if (simParams->haveBeadSize()) {
      sigma_ = simParams->getBeadSize();
    }
  }

  struct BeadLattice {
    Vector3d origin;
    RealType radius;
    bool interior;
  };

  struct ExteriorFunctor : public std::unary_function<BeadLattice, bool>{

    bool operator() (const BeadLattice& bead) {
      return !bead.interior;
    }

  };

  struct InteriorFunctor  : public std::unary_function<BeadLattice, bool>{

    bool operator() (const BeadLattice& bead) {
      return bead.interior;
    }

  };

   //OpenMD files: functions
  //ShapeBuilder.cpp: currShape = new Sphere(atom->getPos(), lja.getSigma()/2.0);
  //Sphere.cpp: Sphere::Sphere(Vector3d origin, RealType radius) : origin_(origin), radius_(radius)
  //std::pair<Vector3d, Vector3d> Sphere::getBoundingBox()
  bool RoughShell::createBeads(std::vector<BeadParam>& beads) {
    std::pair<Vector3d, Vector3d> boxBoundary = shape_->getBoundingBox();
    RealType firstMin = std::min(std::min(floor(boxBoundary.first[0]),
                                          floor(boxBoundary.first[1])),
                                 floor(boxBoundary.first[2]));
    RealType secondMax = std::max(std::max(ceil(boxBoundary.second[0]),
                                           ceil(boxBoundary.second[1])),
                                  ceil(boxBoundary.second[2]));   //std::max is a binary function, i.e., it accepts two values at a time
    int len = secondMax - firstMin;   //floor() and ceil() functions return an integer number

    //number of lattices after shifting positions to integer values (to be symmetric under reflection)
    int numLattices = ceil(len/sigma_) + 2;   //+2: one extra lattice to the left (-1) and one extra lattice to the right (+1)
    Grid3D<BeadLattice>  grid(numLattices, numLattices, numLattices);   //Grid3d.hpp

    //fill beads
    for (int i = 0; i < numLattices; ++i) {
      for (int j = 0; j < numLattices; ++j) {
        for (int k = 0; k < numLattices; ++k) {
          BeadLattice& currentBead = grid(i, j, k);
          currentBead.origin = Vector3d((i-1)*sigma_ + floor(boxBoundary.first[0]),
                                        (j-1)*sigma_ + floor(boxBoundary.first[1]),
                                        (k-1)*sigma_ + floor(boxBoundary.first[2]));
          currentBead.radius = sigma_ / 2.0;
          currentBead.interior = shape_->isInterior(grid(i, j, k).origin);   //return True if bead is inside the given radius (bool Sphere::isInterior(Vector3d pos))
        }
      }
    }

    //remove embedded beads
    for (int i = 0; i < numLattices; ++i) {
      for (int j = 0; j < numLattices; ++j) {
        for (int k = 0; k < numLattices; ++k) {
          std::vector<BeadLattice> neighborCells = grid.getAllNeighbors(i, j, k);
          //if one of its neighbor cells (beads) is exterior, current cell (bead) is
          //on the surface; loop over beads' center (lattice point)

          if (grid(i, j, k).interior){

            bool allNeighBorIsInterior = true;
            for (std::vector<BeadLattice>::iterator l = neighborCells.begin();
                 l != neighborCells.end(); ++l) {
              if (!l->interior) {
                allNeighBorIsInterior = false;
                break;
              }
            }

            if (allNeighBorIsInterior)
              continue;   //if allNeighBorIsInterior = true, skip the remaining code below, i.e., bead is not on the surface

            BeadParam surfaceBead;   //if allNeighBorIsInterior=false; bead is on the surface; grab lattice point
            surfaceBead.atomName = "H";
            surfaceBead.pos = grid(i, j, k).origin;   //loop over the i, j, k positions (lattice point) of the grid (loop is above)
            surfaceBead.mass = 1.0079;
            surfaceBead.radius = grid(i, j, k).radius;
            beads.push_back(surfaceBead);

          }
        }
      }
    }
    return true;
  }
}
