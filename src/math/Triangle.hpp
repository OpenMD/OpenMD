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

#ifndef MATH_FACET_HPP
#define MATH_FACET_HPP

#include <config.h>

#include <vector>

#include "math/SquareMatrix3.hpp"
#include "math/Vector3.hpp"
#include "primitives/StuntDouble.hpp"

namespace OpenMD {

  /**
   * @class Triangle
   *
   * Triangle provides geometric data to OpenMD. Triangle includes
   * information about the normal, centroid and the atoms
   * that belong to this triangle.
   */
  class Triangle {
  public:
    Triangle();
    Triangle(Vector3d P1, Vector3d P2, Vector3d P3);
    virtual ~Triangle() {};

    void setUnitNormal(Vector3d normal) {
      unitnormal_     = normal;
      HaveUnitNormal_ = true;
    }

    void addVertices(Vector3d P1, Vector3d P2, Vector3d P3);

    void addVertexSD(StuntDouble* thisSD) { vertexSD_.push_back(thisSD); }

    std::vector<StuntDouble*> getVertices() { return vertexSD_; }

    void setArea(RealType area) {
      area_     = area;
      HaveArea_ = true;
    }

    Vector3d getNormal() {
      if (HaveNormal_) {
        return normal_;
      } else {
        return computeNormal();
      }
    }
    Vector3d getUnitNormal() {
      if (HaveUnitNormal_) {
        return unitnormal_;
      } else {
        return computeUnitNormal();
      }
    }

    void flipNormal() {
      std::swap(vertices_[1], vertices_[2]);

      // Compute some quantites like a,b,c
      a_ = vertices_[0] - vertices_[1];
      b_ = vertices_[0] - vertices_[2];
      c_ = vertices_[1] - vertices_[2];

      HaveUnitNormal_ = false;
      HaveNormal_     = false;
    }

    RealType getArea() {
      if (HaveArea_) {
        return area_;
      } else {
        return computeArea();
      }
    }

    RealType computeArea();
    Vector3d computeNormal();
    Vector3d computeCentroid();
    Vector3d computeUnitNormal();

    void setCentroid(Vector3d centroid) {
      centroid_     = centroid;
      HaveCentroid_ = true;
    }

    Vector3d getCentroid() {
      if (HaveCentroid_) {
        return centroid_;
      } else {
        return computeCentroid();
      }
    }

    Vector3d getFacetVelocity() { return facetVelocity_; }

    void setFacetVelocity(Vector3d facetVelocity) {
      facetVelocity_ = facetVelocity;
    }

    void setFacetMass(RealType mass) { mass_ = mass; }

    RealType getFacetMass() { return mass_; }

    Vector3d vertex1() const { return vertices_[0]; }
    Vector3d vertex2() const { return vertices_[1]; }
    Vector3d vertex3() const { return vertices_[2]; }

    RealType a() { return a_.length(); }

    RealType b() { return b_.length(); }

    RealType c() { return c_.length(); }

    RealType getHydroLength() {
      RealType a1 = a();
      RealType b1 = b();
      RealType c1 = c();
      RealType t1 = a1 + b1 + c1;
      RealType t4 = a1 + b1 - c1;

      return 32.0 * c1 / log(t1 * t1 / t4 / t4);
    }

    RealType getIncircleRadius() { return 2.0 * getArea() / (a() + b() + c()); }

    RealType getCircumcircleRadius() {
      RealType a1 = a();
      RealType b1 = b();
      RealType c1 = c();
      RealType t1 = a1 + b1 + c1;
      RealType t2 = -a1 + b1 + c1;
      RealType t3 = a1 - b1 + c1;
      RealType t4 = a1 + b1 - c1;
      return a1 * b1 * c1 / sqrt(t1 * t2 * t3 * t4);
    }

    Mat3x3d computeHydrodynamicTensor(RealType viscosity);

    Vector3d cartesionToBarycentric(Vector3d p) {
      RealType area = getArea();

      Vector3d v0 = vertices_[0] - p;
      Vector3d v1 = vertices_[1] - p;
      Vector3d v2 = vertices_[2] - p;

      RealType u = 0.5 * cross(v1, v2).length() / area;
      RealType v = 0.5 * cross(v0, v2).length() / area;
      RealType w = 0.5 * cross(v0, v1).length() / area;

      return Vector3d(u, v, w);
    }

    Vector3d barycentricToCartesian(const Vector3d& barycentric) const {
      return barycentric.x() * vertices_[0] + barycentric.y() * vertices_[1] +
             barycentric.z() * vertices_[2];
    }

  private:
    Mat3x3d hydro_tensor(const Vector3d& ri, const Vector3d& rj0,
                         const Vector3d& rj1, const Vector3d& rj2, RealType s,
                         RealType viscosity);

    /* Local Indentity of vertex atoms in pos array*/
    std::vector<StuntDouble*> vertexSD_;
    Vector3d normal_;
    Vector3d unitnormal_;
    Vector3d centroid_;
    Vector3d vertices_[3];
    RealType area_;
    RealType mass_;
    Vector3d facetVelocity_;
    // Length of triangle sides
    Vector3d a_, b_, c_;
    bool HaveArea_;
    bool HaveNormal_;
    bool HaveUnitNormal_;
    bool HaveCentroid_;

  };  // End class Triangle

}  // namespace OpenMD

#endif  // MATH_FACET_HPP
