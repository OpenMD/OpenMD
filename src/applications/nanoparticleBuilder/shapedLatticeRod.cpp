/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
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
 * research, please cite the following paper when you publish your work:
 *
 * [1] Drisko et al., J. Open Source Softw. 9, 7004 (2024).
 *
 * Good starting points for code and simulation methodology are:
 *
 * [2] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [3] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [4] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [7] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 * [9] Drisko & Gezelter, J. Chem. Theory Comput. 20, 4986-4997 (2024).
 */

/*
 *  Created by Kelsey M. Stocker on 2/9/12.
 *  @author  Kelsey M. Stocker
 *  @version $Id: shapedLatticeRod.cpp 1665 2011-11-22 20:38:56Z gezelter $
 *
 */

#include "shapedLatticeRod.hpp"

#include <cmath>

#include "lattice/shapedLattice.hpp"

using namespace std;
namespace OpenMD {

  shapedLatticeRod::shapedLatticeRod(RealType latticeConstant,
                                     std::string latticeType, RealType radius,
                                     RealType length) :
      shapedLattice(latticeConstant, latticeType) {
    rodRadius_ = radius;
    rodLength_ = length;
    Vector3d dimension;
    dimension[0] = 3.0 * radius;
    dimension[1] = 3.0 * radius;
    dimension[2] = length + 3.0 * radius;
    setGridDimension(dimension);
    Vector3d origin;
    origin[0] = 0;
    origin[1] = 0;
    origin[2] = 0;
    setOrigin(origin);
  }

  /**
   * Determines whether a point lies within a spherically-capped
   * nanorod at origin (0,0,0)
   *
   */
  bool shapedLatticeRod::isInterior(Vector3d point) {
    RealType x, y, z, distance;

    distance = 0;

    x = point[0];
    y = point[1];
    z = point[2];

    if (abs(z) >= rodLength_ / 2.0) {
      RealType delta_z = abs(z) - rodLength_ / 2.0;
      distance         = sqrt((x * x) + (y * y) + (delta_z * delta_z));
    } else {
      distance = sqrt((x * x) + (y * y));
    }

    bool isIT = false;
    if (distance <= rodRadius_) { isIT = true; }
    return isIT;
  }

}  // namespace OpenMD
