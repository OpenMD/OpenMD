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

#include "hydrodynamics/BoundaryElementModel.hpp"

#include "hydrodynamics/Mesh.hpp"
#include "math/DynamicRectMatrix.hpp"
#include "math/LU.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/integration/StrangFixCowperTriangleQuadrature.hpp"
#include "math/integration/TriangleQuadrature.hpp"
#include "utils/Constants.hpp"
#include "utils/simError.h"

namespace OpenMD {

  BoundaryElementModel::BoundaryElementModel() : ApproximateModel() {}

  std::size_t BoundaryElementModel::assignElements() {
    if (shape_ != NULL) {
      if (shape_->isMesh()) {
        createTriangles(dynamic_cast<Mesh*>(shape_));
      } else {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "BoundaryElementModel::assignElements Error: No mesh was "
                 "given as the shape\n");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal  = 1;
        simError();
      }
      return elements_.size();
    }
    return 0;
  }

  void BoundaryElementModel::createTriangles(Mesh* m) {
    if (m != NULL) {
      std::string name             = m->getName();
      std::vector<Triangle> facets = m->getFacets();
      for (std::vector<Triangle>::iterator i = facets.begin();
           i != facets.end(); ++i) {
        HydrodynamicsElement currTri;
        currTri.name = name;
        currTri.pos  = (*i).getCentroid();
        currTri.t    = (*i);
        elements_.push_back(currTri);
      }
    }
  }

  void BoundaryElementModel::checkElement(std::size_t i) {}

  void BoundaryElementModel::writeElements(std::ostream& os) {
    std::vector<HydrodynamicsElement>::iterator iter;
    std::string name = shape_->getName();
    os << "solid"
       << " " << name << std::endl;
    for (iter = elements_.begin(); iter != elements_.end(); ++iter) {
      Triangle t  = iter->t;
      Vector3d n  = t.getUnitNormal();
      Vector3d v1 = t.vertex1();
      Vector3d v2 = t.vertex2();
      Vector3d v3 = t.vertex3();
      os << "  "
         << "facet normal"
         << " " << n[0] << " " << n[1] << " " << n[2] << std::endl;
      os << "    "
         << "outer loop" << std::endl;
      os << "      "
         << " "
         << "vertex"
         << " " << v1[0] << " " << v1[1] << " " << v1[2] << std::endl;
      os << "      "
         << " "
         << "vertex"
         << " " << v2[0] << " " << v2[1] << " " << v2[2] << std::endl;
      os << "      "
         << " "
         << "vertex"
         << " " << v3[0] << " " << v3[1] << " " << v3[2] << std::endl;
      os << "    "
         << "endloop" << std::endl;
      os << "  "
         << "endfacet" << std::endl;
    }
    os << "endsolid"
       << " " << name << std::endl;
  }

  Mat3x3d BoundaryElementModel::interactionTensor(const std::size_t i,
                                                  const std::size_t j,
                                                  const RealType viscosity) {
    Mat3x3d B;
    Mat3x3d I = SquareMatrix3<RealType>::identity();

    StrangFixCowperTriangleQuadratureRule rule(6);

    Vector3d centroid = elements_[i].pos;
    Triangle t        = elements_[j].t;

    auto Tij = [&t, &centroid, &I, &viscosity](const Vector3d& p) {
      // p are in barycentric coordinates
      Vector3d r   = t.barycentricToCartesian(p);
      Vector3d ab  = centroid - r;
      RealType abl = ab.length();
      Mat3x3d T;
      T = (I + outProduct(ab, ab) / (abl * abl));
      T /= (8.0 * Constants::PI * viscosity * abl);
      return T;
    };

    B = TriangleQuadrature<RectMatrix<RealType, 3, 3>, RealType>::Integrate(
        Tij, rule, 1.0);

    centroid = elements_[j].pos;
    t        = elements_[i].t;

    auto Tji = [&t, &centroid, &I, &viscosity](const Vector3d& p) {
      // p are in barycentric coordinates
      Vector3d r   = t.barycentricToCartesian(p);
      Vector3d ab  = centroid - r;
      RealType abl = ab.length();
      Mat3x3d T;
      T = (I + outProduct(ab, ab) / (abl * abl));
      T /= (8.0 * Constants::PI * viscosity * abl);
      return T;
    };

    B += TriangleQuadrature<RectMatrix<RealType, 3, 3>, RealType>::Integrate(
        Tji, rule, 1.0);
    B *= 0.5;
    return B;
  }
}  // namespace OpenMD
