/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
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

#include "hydrodynamics/RoughShell.hpp"

#include "brains/SimInfo.hpp"
#include "hydrodynamics/ShapeBuilder.hpp"

namespace OpenMD {

  struct BeadLattice {
    Vector3d origin;
    RealType radius;
    bool interior;
  };

  struct ExteriorFunctor {
    bool operator()(const BeadLattice& bead) { return !bead.interior; }
  };

  struct InteriorFunctor {
    bool operator()(const BeadLattice& bead) { return bead.interior; }
  };

  std::size_t RoughShell::assignElements() {
    std::pair<Vector3d, Vector3d> boxBoundary = shape_->getBoundingBox();
    RealType firstMin                         = std::min(
        std::min(floor(boxBoundary.first[0]), floor(boxBoundary.first[1])),
        floor(boxBoundary.first[2]));
    RealType secondMax = std::max(
        std::max(ceil(boxBoundary.second[0]), ceil(boxBoundary.second[1])),
        ceil(boxBoundary.second[2]));
    // std::max is a binary function, i.e., it accepts two values at a time
    // floor() and ceil() functions return an integer number
    int len = secondMax - firstMin;

    // number of lattices after shifting positions to integer values (to be
    // symmetric under reflection)

    int numLattices = ceil(len / sigma_) + 2;

    // +2: one extra lattice point to the left (-1) and one extra
    // lattice point to the right (+1)

    Grid3D<BeadLattice> grid(numLattices, numLattices, numLattices);

    // fill beads
    for (int i = 0; i < numLattices; ++i) {
      for (int j = 0; j < numLattices; ++j) {
        for (int k = 0; k < numLattices; ++k) {
          BeadLattice& currentBead = grid(i, j, k);
          currentBead.origin =
              Vector3d((i - 1) * sigma_ + floor(boxBoundary.first[0]),
                       (j - 1) * sigma_ + floor(boxBoundary.first[1]),
                       (k - 1) * sigma_ + floor(boxBoundary.first[2]));
          currentBead.radius   = sigma_ / 2.0;
          currentBead.interior = shape_->isInterior(grid(i, j, k).origin);
          // returns True if bead is inside the given shape
        }
      }
    }

    // remove embedded beads
    for (int i = 0; i < numLattices; ++i) {
      for (int j = 0; j < numLattices; ++j) {
        for (int k = 0; k < numLattices; ++k) {
          std::vector<BeadLattice> neighborCells =
              grid.getAllNeighbors(i, j, k);
          // if one of its neighbor cells (beads) is exterior, current cell
          // (bead) is on the surface; loop over beads' center (lattice point)

          if (grid(i, j, k).interior) {
            bool allNeighborsAreInterior = true;
            for (std::vector<BeadLattice>::iterator l = neighborCells.begin();
                 l != neighborCells.end(); ++l) {
              if (!l->interior) {
                allNeighborsAreInterior = false;
                break;
              }
            }

            if (allNeighborsAreInterior)
              continue;  // if allNeighborsAreInterior == true, skip
                         // the remaining code below, i.e., bead is
                         // not on the surface

            HydrodynamicsElement surfaceBead;
            // if allNeighborsAreInterior == false, bead is on the
            // surface; grab this lattice point
            surfaceBead.name = "H";
            surfaceBead.pos  = grid(i, j, k).origin;
            // loop over the i, j, k positions (lattice point) of the
            // grid (loop is above)
            surfaceBead.radius = grid(i, j, k).radius;
            elements_.push_back(surfaceBead);
          }
        }
      }
    }
    return elements_.size();
  }
}  // namespace OpenMD
