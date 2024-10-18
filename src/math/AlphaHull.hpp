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

#ifndef MATH_ALPHAHULL_HPP_
#define MATH_ALPHAHULL_HPP_

#include <config.h>

#include <cassert>
#include <string>
#include <vector>

#include "math/Hull.hpp"
#include "math/Triangle.hpp"
#include "math/Vector3.hpp"

namespace OpenMD {

  /**
   * @class AlphaHull
   * @brief Compute alpha complex or alpha shape
   *
   * Builds the alpha shape (H. Edelsbrunner and P.Mucke,
   * "Three-dimensional Alpha Shapes," ACM Trans. Graph. 13, 1994) from a
   * set of atomic locations using the Qhull library,
   * http://www.qhull.org/
   *
   * For a given value of \f$\alpha\f$, the \f$\alpha\f$-shape
   * includes all of the tetrahedra in the Delaunay triangulation
   * which have an empty circumsphere with radius equal or smaller
   * than \f$\alpha\f$.  To carry out this calculation, all points are
   * lifted to 4D so that each point is \f$(x, y, z, x^2+y^2+z^2)\f$.
   * The convex hull is computed in 4D, and tetrahedral facets with
   * circumsphere radii larger than \f$\alpha\f$ are removed.
   *
   *   \param alpha the circumsphere radius to test tetrahedra for elimination
   */
  class AlphaHull : public Hull {
  public:
    AlphaHull(RealType alpha);
    virtual ~AlphaHull() {};

    void computeHull(std::vector<StuntDouble*> bodydoubles);

    /* Total area of Hull*/
    RealType getArea() { return area_; }

    /* Total Volume enclosed by Hull */
    RealType getVolume() { return volume_; }

    vector<Triangle> getMesh() { return Triangles_; }

  protected:
    int dim_;
    RealType alpha_;
    const std::string options_;

  private:
    // These variables are private so that each new hull returns
    // information about itself.
    RealType volume_;
    RealType area_;
    std::vector<Triangle> Triangles_;
  };
}  // namespace OpenMD

#endif
