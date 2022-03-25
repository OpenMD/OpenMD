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

#include "shapedLatticePentRod.hpp"

#include <cmath>

#include "lattice/shapedLattice.hpp"
#include "shapedLatticeRod.hpp"

using namespace std;
namespace OpenMD {

  shapedLatticePentRod::shapedLatticePentRod(RealType latticeConstant,
                                             std::string latticeType,
                                             RealType radius, RealType length) :
      shapedLattice(latticeConstant, latticeType) {
    rodRadius_ = radius;
    rodLength_ = length;
    Vector3d dimension;
    dimension[0] = 2.0 * length;
    dimension[1] = 2.0 * length;
    dimension[2] = 2.0 * length;
    setGridDimension(dimension);
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

  bool shapedLatticePentRod::isInterior(Vector3d point) {
    RealType x, y, z, new_x, new_y, left_newx, right_newx;
    // int z_int;

    bool isIT = false;

    x = point[0];
    y = point[1];
    z = point[2];

    // z_int = int(z/2.04);

    // Rotate by 45 degrees around z-axis so length of rod lies along y axis
    new_x = (sqrt(2.0) / 2) * (x - y);
    new_y = (sqrt(2.0) / 2) * (x + y);

    left_newx  = (z - 1.44) * (0.577350269 / 0.816496581);
    right_newx = (z + 1.44) * (-0.577350269 / 0.816496581);

    // Make center spine of nanorod down new_y axis
    // This is now done directly in nanorod_pentBuilder.cpp
    /*if ( (new_x == 0) && (z == 0) ) {

    if ( abs(new_y) <= 0.5*rodLength_ + 0.5773502692*rodRadius_ ) {

      isIT = true;

    }
    }*/

    // Make one wedge
    if ((z < 0) && (z >= -0.816496581 * rodRadius_ - 1.44)) {
      if (abs(new_y) <=
          1.44 * (z / 2.04) + 0.5 * rodLength_ + 0.5773502692 * rodRadius_) {
        if ((new_x >= left_newx) && (new_x <= right_newx)) { isIT = true; }
      }
    }
    return isIT;
  }
}  // namespace OpenMD
