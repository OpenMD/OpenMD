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

#ifndef MATH_ALPHAHULL_HPP_
#define MATH_ALPHAHULL_HPP_

#include "math/Vector3.hpp"
#include "config.h"
#include "math/Hull.hpp"
#include "math/Triangle.hpp"

#include <cassert>
#include <vector>
#include <string>

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
   * circumsphere radii larger than alpha are removed.
   *
   *   \param alpha the circumsphere radius to test tetrahedra for elimination
   */
  class AlphaHull : public Hull {
  public:
    
    AlphaHull(RealType alpha);    
    virtual ~AlphaHull(){};
    
    void computeHull( std::vector<StuntDouble*> bodydoubles );
    
    /* Total area of Hull*/
    RealType getArea(){ return area_; }
    
    /* Total Volume enclosed by Hull */
    RealType getVolume(){ return volume_; } 
    
    vector<Triangle> getMesh(){ return Triangles_; }
    
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
}
#endif
