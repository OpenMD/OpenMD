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

#include "math/Triangle.hpp"

using namespace OpenMD;

Triangle::Triangle() :
    normal_(V3Zero), centroid_(V3Zero), area_(0.0), mass_(0.0),
    facetVelocity_(V3Zero), a_(V3Zero), b_(V3Zero), c_(V3Zero),
    HaveArea_(false), HaveNormal_(false), HaveUnitNormal_(false),
    HaveCentroid_(false) {}

Triangle::Triangle(Vector3d P1, Vector3d P2, Vector3d P3) :
    mass_(0.0), facetVelocity_(V3Zero), HaveArea_(false), HaveNormal_(false),
    HaveUnitNormal_(false), HaveCentroid_(false) {
  vertices_[0] = P1;
  vertices_[1] = P2;
  vertices_[2] = P3;

  // Compute some quantites like a,b,c
  a_ = P1 - P2;
  b_ = P1 - P3;
  c_ = P2 - P3;
  computeArea();
  computeNormal();
  computeCentroid();
}

void Triangle::addVertices(Vector3d P1, Vector3d P2, Vector3d P3) {
  vertices_[0] = P1;
  vertices_[1] = P2;
  vertices_[2] = P3;

  // Compute some quantites like a,b,c
  a_ = P1 - P2;
  b_ = P1 - P3;
  c_ = P2 - P3;
}

RealType Triangle::computeArea() {
  HaveArea_ = true;
  area_     = getNormal().length() * 0.5;
  return area_;
}
// This should return the normal for our calculations.
Vector3d Triangle::computeNormal() {
  HaveNormal_ = true;
  normal_     = cross(a_, b_);
  return normal_;
}
// This should return the normal for our calculations.
Vector3d Triangle::computeUnitNormal() {
  HaveUnitNormal_ = true;
  unitnormal_     = cross(a_, b_);
  unitnormal_.normalize();
  return unitnormal_;
}

Vector3d Triangle::computeCentroid() {
  HaveCentroid_ = true;
  centroid_     = (vertices_[0] + vertices_[1] + vertices_[2]) / RealType(3.0);
  return centroid_;
}

Mat3x3d Triangle::computeHydrodynamicTensor(RealType viscosity) {
  Vector3d u0 = -a_;
  Vector3d v0 = centroid_ - vertices_[0];
  RealType s0 = 0.5 * cross(u0, v0).length();

  Vector3d u1 = -c_;
  Vector3d v1 = centroid_ - vertices_[1];
  RealType s1 = 0.5 * cross(u1, v1).length();

  Vector3d u2 = b_;
  Vector3d v2 = centroid_ - vertices_[2];
  RealType s2 = 0.5 * cross(u2, v2).length();

  Mat3x3d H;
  H = hydro_tensor(centroid_, centroid_, vertices_[1], vertices_[0], s0,
                   viscosity) +
      hydro_tensor(centroid_, centroid_, vertices_[1], vertices_[2], s1,
                   viscosity) +
      hydro_tensor(centroid_, centroid_, vertices_[2], vertices_[0], s2,
                   viscosity);
  return H.inverse();
}

Mat3x3d Triangle::hydro_tensor(const Vector3d& ri, const Vector3d& rj0,
                               const Vector3d& rj1, const Vector3d& rj2,
                               RealType s, RealType viscosity) {
  Vector3d v2  = (rj0 + rj1 + rj2) / RealType(3.0);  // sub-centroid
  Vector3d dr  = ri - v2;  // real centroid to sub-centroid
  RealType ri2 = RealType(1.0) / dr.lengthSquare();

  Mat3x3d G;
  G = (SquareMatrix3<RealType>::identity() + outProduct(dr, dr) * ri2) *
      RealType(sqrt(ri2));

  G *= 1.0 / (8.0 * 3.14159285358979 * viscosity);
  return G;
}
